#!/usr/bin/env python3
"""Compare binary128 trigonometric builds with one GCC-built timing driver.

Run from any directory. Builds, fixed input files, disassembly, commands and
raw RDTSC samples are retained in the output directory; the checkout is not
built in place. Times are TSC ticks/call, not hardware core-cycle counts.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--cpu", type=int, required=True)
    parser.add_argument("--count", type=int, default=100000)
    parser.add_argument("--repeat", type=int, default=10)
    parser.add_argument("--trials", type=int, default=10)
    parser.add_argument("--functions", nargs="+", choices=["sinq", "cosq", "sincosq", "tanq"],
                        required=True)
    parser.add_argument("--compilers", nargs="+", choices=["gcc", "clang", "icx"],
                        default=["gcc", "clang", "icx"])
    parser.add_argument("--modes", nargs="+", choices=["rndn", "rndz", "rndu", "rndd"],
                        default=["rndn", "rndz", "rndu", "rndd"])
    parser.add_argument("--metrics", nargs="+", choices=["throughput", "latency"],
                        default=["throughput", "latency"])
    parser.add_argument("--icx", default="icx")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bands", nargs="+", choices=["default", "full", "small"],
                        default=["default"])
    args = parser.parse_args()
    if min(args.count, args.repeat, args.trials) <= 0:
        parser.error("count, repeat and trials must be positive")
    if args.cpu not in os.sched_getaffinity(0):
        parser.error("CPU is not in the current affinity mask")
    root = Path(__file__).resolve().parents[3]
    out = args.output or Path(tempfile.mkdtemp(prefix="core-math-compilers-"))
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    print(f"Artifacts: {out}", flush=True)
    commands = []

    def run(command):
        command = list(map(str, command))
        commands.append(command)
        return subprocess.run(command, cwd=root, check=True, text=True,
                              stdout=subprocess.PIPE).stdout

    compilers = {name: shutil.which(args.icx if name == "icx" else name)
                 for name in args.compilers}
    if not shutil.which("gcc") or not all(compilers.values()):
        parser.error("GCC and every selected compiler must be available; use --icx PATH")
    flags = {name: ["-O3", "-march=native", "-fno-pie"] +
             (["-fp-model=precise"] if name == "icx" else ["-frounding-math"])
             for name in args.compilers}
    # icx supplies its conversion/memory runtime helpers; every executable
    # still has the same GCC driver and is linked by the same linker driver.
    linker = compilers.get("icx", "gcc")
    record = {"baseline": run(["git", "rev-parse", args.baseline]).strip(),
              "cpu": args.cpu, "cpu_info": run(["lscpu"]),
              "count": args.count, "repeat": args.repeat, "trials": args.trials,
              "flags": flags, "versions": {c: run([p, "--version"])
                                           for c, p in compilers.items()},
              "functions": {}}
    for fun in args.functions:
        folder = root / "src/binary128" / fun[:-1]
        driver = out / f"{fun}-driver.o"
        harness = root / "src/generic/support/perf_univariate.c"
        if fun == "sincosq":
            # Match the binary128 univariate GCC driver's bit dependency.
            # The stock paired driver uses software binary128 arithmetic.
            text = (root / "src/generic/support/perf_sincos.c").read_text()
            macro = r"""#define CALL_LATENCY(accu,i) do { \
  TYPE_UNDER_TEST x, y; \
  __m128i in, out, m = {1, 0}; \
  asm("" : "=x" (in) : "0" (accu)); \
  asm("" : "=x" (out) : "0" (randoms[i])); \
  in = _mm_and_si128(in, m); \
  out = _mm_or_si128(out, in); \
  asm("" : "=x" (accu) : "0" (out)); \
  p_function_under_test(accu, &x, &y); \
  accu = x; \
} while (0)"""
            text, n = re.subn(r"#define CALL_LATENCY.*?while \(0\)",
                              lambda _: macro, text, flags=re.S)
            if n != 1:
                raise RuntimeError("unrecognized sincos latency driver")
            harness = out / "perf_sincos.c"
            harness.write_text(text)
        run(["gcc", "-O3", "-march=native", "-frounding-math", "-DCORE_MATH_F128",
             "-DTYPE_UNDER_TEST=__float128", f"-Dcr_function_under_test=cr_{fun}",
             f"-Dfunction_under_test={fun}", "-I", root / "src/generic" / fun[:-1],
             "-I", root / "src/generic/support", "-c", harness, "-o", driver])
        sources = {
            "baseline": run(["git", "show", f"{args.baseline}:src/binary128/{fun[:-1]}/{fun}.c"]),
            "candidate": (folder / f"{fun}.c").read_text(),
        }
        executables = []
        for version, source in sources.items():
            path = out / f"{fun}-{version}.c"
            path.write_text(source)
            for compiler, executable in compilers.items():
                stem = out / f"{fun}-{version}-{compiler}"
                obj = stem.with_suffix(".o")
                run([executable, *flags[compiler], "-c", path, "-o", obj])
                run([linker, "-no-pie", driver, obj, "-lquadmath", "-lm", "-o", stem])
                stem.with_suffix(".asm").write_text(run(["objdump", "-drwC", "-Mintel", obj]))
                executables.append(stem)
        data = {"input_sha256": {},
                "source_sha256": {v: hashlib.sha256(s.encode()).hexdigest()
                                  for v, s in sources.items()}, "samples": {}}
        for band in args.bands:
            inputs = out / f"{fun}-{band}.inputs"
            if band == "default":
                run([executables[0], "--file", inputs, "--count", args.count, "--reference"])
            else:
                rng = random.Random(1)
                upper = 16383 if band == "full" else -9
                with inputs.open("wb") as stream:
                    for _ in range(args.count):
                        bits = rng.getrandbits(1) << 127
                        bits |= (rng.randint(-57, upper) + 16383) << 112
                        bits |= rng.getrandbits(112)
                        stream.write(bits.to_bytes(16, sys.byteorder))
            data["input_sha256"][band] = hashlib.sha256(inputs.read_bytes()).hexdigest()
            for metric in args.metrics:
                for mode in args.modes:
                    samples = {exe.name: [] for exe in executables}
                    for trial in range(args.trials):
                        offset = trial % len(executables)
                        for exe in executables[offset:] + executables[:offset]:
                            command = ["taskset", "-c", args.cpu, exe, "--file", inputs,
                                       "--count", args.count, "--repeat", args.repeat,
                                       "--rdtsc", "--" + mode]
                            if metric == "latency":
                                command.append("--latency")
                            samples[exe.name].append(float(run(command)))
                    data["samples"][band + ":" + metric + ":" + mode] = samples
                    for name, values in samples.items():
                        print(f"{name} {band} {metric} {mode}: min {min(values):.3f}, "
                              f"median {statistics.median(values):.3f}", flush=True)
        record["functions"][fun] = data
        (out / "results.json").write_text(json.dumps(record, indent=2) + "\n")
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")


if __name__ == "__main__":
    main()
