#!/usr/bin/env python3
"""Compare the working tanq.c with a commit using the repository perf harness.

Example (all modes, throughput and latency, three input distributions):
  python3 src/binary128/tan/bench_compare.py --baseline fb944b75 --cpu 4
CC and CFLAGS select the compiler and identical flags for both versions.
No correctness claim follows from this benchmark.
"""
import argparse
import os
from pathlib import Path
import random
import shlex
import statistics
import subprocess
import sys
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--baseline', required=True)
parser.add_argument('--cpu', type=int)
parser.add_argument('--count', type=int, default=20000)
parser.add_argument('--repeat', type=int, default=5)
parser.add_argument('--trials', type=int, default=3)
parser.add_argument('--bands', nargs='+', choices=['default', 'full', 'small'],
                    default=['default', 'full', 'small'])
parser.add_argument('--force-wide', action='store_true')
args = parser.parse_args()
if min(args.count, args.repeat, args.trials) <= 0:
    parser.error('count, repeat and trials must be positive')
root = Path(__file__).resolve().parents[3]
source = Path('src/binary128/tan/tanq.c')
cc = shlex.split(os.environ.get('CC', 'cc'))
flags = shlex.split(os.environ.get('CFLAGS',
    '-O3 -march=native -frounding-math -fno-finite-math-only'))
pin = ['taskset', '-c', str(args.cpu)] if args.cpu is not None else []

def run(cmd):
    return subprocess.check_output(cmd, cwd=root, text=True).strip()

with tempfile.TemporaryDirectory(prefix='tanq-bench.') as tmp:
    work = Path(tmp)
    before = work/'before.c'
    before.write_text(run(['git', 'show', f'{args.baseline}:{source}'])+'\n')
    revision = run(['git', 'rev-parse', args.baseline])
    print(f'baseline={revision} cpu={args.cpu} compiler={shlex.join(cc)}', flush=True)
    print(f'flags={shlex.join(flags)} count={args.count} repeat={args.repeat} '
          f'trials={args.trials} forced_wide={args.force_wide}', flush=True)
    executables = {}
    for name, src in [('before', before), ('after', root/source)]:
        obj, exe = work/f'{name}.o', work/name
        force = ['-DCORE_MATH_TANQ_FORCE_ACCURATE'] if args.force_wide else []
        run(cc+flags+force+['-c', str(src), '-o', str(obj)])
        run(cc+flags+['-DCORE_MATH_F128', '-DTYPE_UNDER_TEST=__float128',
            '-Dcr_function_under_test=cr_tanq', '-Dfunction_under_test=tanq',
            '-Isrc/generic/tan', 'src/generic/support/perf_univariate.c',
            str(obj), '-lquadmath', '-lm', '-o', str(exe)])
        executables[name] = str(exe)
    print('band,mode,kind,before_median,after_median,speedup,before_samples,after_samples', flush=True)
    for band in args.bands:
        inputs = work/f'{band}.dat'
        if band == 'default':
            run([executables['before'], '--file', str(inputs), '--reference',
                 '--count', str(args.count)])
        else:
            rng = random.Random(1)
            # Exact binary128 encodings, including random low significand bits.
            # Full: unbiased exponents [-57,16383]; small: [-57,-9].
            upper = 16383 if band == 'full' else -9
            with inputs.open('wb') as f:
                for _ in range(args.count):
                    bits = (rng.getrandbits(1) << 127)
                    bits |= (rng.randint(-57, upper)+16383) << 112
                    bits |= rng.getrandbits(112)
                    f.write(bits.to_bytes(16, sys.byteorder))
        for mode in ('rndn', 'rndz', 'rndu', 'rndd'):
            for kind in ('throughput', 'latency'):
                options = ['--file', str(inputs), '--count', str(args.count),
                           '--repeat', str(args.repeat), '--rdtsc', '--'+mode]
                if kind == 'latency':
                    options.append('--latency')
                samples = {'before': [], 'after': []}
                # Alternate A/B and B/A pairs to reduce ordering bias.
                for trial in range(args.trials):
                    order = ('before', 'after') if trial % 2 == 0 else ('after', 'before')
                    for name in order:
                        samples[name].append(float(run(pin+[executables[name]]+options)))
                a, b = (statistics.median(samples[x]) for x in ('before', 'after'))
                raw_a = ';'.join(map(str, samples['before']))
                raw_b = ';'.join(map(str, samples['after']))
                print(f'{band},{mode},{kind},{a:.3f},{b:.3f},{a/b:.3f},'
                      f'{raw_a},{raw_b}', flush=True)
