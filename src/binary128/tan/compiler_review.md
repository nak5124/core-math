# Compiler review: tanq, 2026-09-21

The branch preserves the existing integer results while expressing reduction
and frame operations in 64-bit limbs. The generated GCC code improves
without changing constants, approximation order, error bounds or rounding
gates. The GCC/icx gap remains; this is not a claim of a compiler defect.

## Changes and scope

The reduction product handles its first and last columns explicitly and
accumulates the low-limb sum separately from the two high halves. This
avoids detecting overflow of a full 128-bit column. The 192-bit right shift
and frame normalization use explicit 64-bit funnels instead of general
128-bit variable shifts. The per-function error analysis gives the exact
integer identities and the shift/range preconditions.

Tangent also splits `shr_round` at the 64-bit boundary, explicitly taking
the folded bit from the low limb when the count is exactly 64. Here that
rewrite improved measurements after removing an experimental multiply
change. The reciprocal/Newton divider is unchanged; trial branch-removal
rewrites did not improve its measured performance. The branch includes
Paul Zimmermann's hard-case commit `5f34efa3` before this optimization.

## Measurement conditions

Intel Core i9-14900K, CPU 4 affinity; GCC 15.2.0, Clang 21.1.8, icx
2026.1.1 (20260724). Times are **RDTSC ticks/call**, not hardware core cycles.
These compiler versions differ from Paul's GCC 16.2.0 and icx 2026.1.0.

Baseline `5f34efa3` and working source use `-O3 -march=native -fno-pie`,
plus `-frounding-math` for GCC/Clang or `-fp-model=precise` for icx.
Every executable has the same GCC-built timing driver, input file and icx
link step (`-no-pie`, including the Intel runtime helpers). Trials rotate
executable order. The paired latency driver uses the same integer bit
dependency as the univariate GCC driver, avoiding the stock paired driver's
software binary128 multiply/add overhead. It measures dependency through
the sine output.

Default inputs use the repository generator over `[-pi,pi]`: 100,000 inputs,
10 repeats, 20 trials. Reported values are medians. Baseline -> candidate:

| Function, nearest throughput | GCC | Clang | icx |
| --- | ---: | ---: | ---: |
| tanq | 187.720 -> 171.977 | 160.865 -> 151.351 | 159.416 -> 151.043 |

Across all four rounding modes and both throughput and latency, the range
of percentage reductions in median time is:

| Function | GCC | Clang | icx |
| --- | ---: | ---: | ---: |
| tanq | 7.4% to 8.4% | 4.8% to 5.9% | 3.6% to 5.3% |

Additional distributions use random signs and all 112 fraction bits, with
uniform exponents in `[-57,16383]` (full) or `[-57,-9]` (small), seed 1.
Each uses 20,000 inputs, 50 repeats and 10 trials. Percentage reductions
across all four modes and both metrics are below; negative values denote
regressions. These results are input- and machine-specific.

| Function / band | GCC | Clang | icx |
| --- | ---: | ---: | ---: |
| tanq / full | 6.9% to 8.2% | 4.5% to 5.8% | 3.4% to 4.8% |
| tanq / small | -1.3% to 6.4% | -0.9% to -0.2% | -0.9% to -0.0% |

## Reproduction

The benchmark saves source snapshots, exact commands, disassembly, input
hashes, metadata and every raw sample outside the checkout:

```sh
python3 src/generic/support/bench_trig_compilers.py \
  --baseline 5f34efa3 --functions tanq --cpu 4 --trials 20 \
  --icx /opt/intel/oneapi/compiler/2026.1/bin/icx

# Additional distributions (same compiler selection):
python3 src/generic/support/bench_trig_compilers.py \
  --baseline 5f34efa3 --functions tanq --cpu 4 \
  --bands full small --count 20000 --repeat 50 --trials 10 \
  --icx /opt/intel/oneapi/compiler/2026.1/bin/icx
```

Use an available CPU and compiler paths. This script targets Linux/x86-64
and requires GCC, the selected function compilers, quadmath, `taskset`, and
`objdump`. `--compilers gcc clang` omits icx when unavailable. Performance
comparisons use CORE-MATH before/after, not a libc correctness oracle.

## Validation of the retained patch

All checks below passed; input samples supplement the existing error
arguments and do not exhaust the binary128 domain.

* Complete regression corpora in all four rounding modes, with GCC-,
  Clang-, and icx-compiled function objects, inexact and errno enabled.
  Parsed input records: `tanq` 61,960.
* The existing special checks at the 100,000 setting for every function
  and compiler, all four modes, inexact/errno enabled. The cross-compiler
  checks use a common GCC-built MPFR reference and harness. Harness copies
  fix the random seed at 1; the tangent harness uses its seed environment
  variable. Special checks include full-range inputs, the band near zero,
  neighborhoods of multiples of pi/2, and exceptional/subnormal inputs.
* GCC ASan+UBSan special checks at the 10,000 setting in all four modes.
* UBSan-enabled old/new helper comparisons: 2,000,000 reduction products
  across both window lengths, plus 1,010,287 frames covering shifts,
  normalization, random limbs, zeros, maximal limbs and carry chains.
* Exact-rational error-bound generation/checks passed. No constants or
  error radii changed. `git diff --check` passed.
* Forced-wide checks passed on all 61,960 corpus inputs and the special
  checks at 10,000, all modes, inexact/errno enabled.
* `check_internal` passed its one million narrow-quotient checks,
  100,000 exact wide divisions, interval checks and quotient-bound checks.

For standard GCC MPFR checks:

```sh
for f in tanq; do
  for kind in worst special; do
    CC=gcc OMP_NUM_THREADS=2 \
      EXTRA_CFLAGS='-Werror -DCORE_MATH_TESTS=100000 -DCORE_MATH_CHECK_INEXACT -DCORE_MATH_SUPPORT_ERRNO' \
      ./check.sh --$kind "$f"
  done
done
```

Paul's compiler versions have not been run for this patch.

## Validation follow-up

Pipeline [1486925](https://gitlab.inria.fr/core-math/core-math/-/pipelines/1486925)
passed at `c03e11cc`. Its ordinary, errno, inexact and no-LZCNT jobs ran
the 61,960-input corpus and special checks. The GCC sanitizer job skipped
binary128 after the combined quadmath/SSE4.1 probe; its green status is
not sanitizer coverage of tangent.

Manual baseline x86-64 runs, with a GCC-built MPFR reference and harness
linked against separately compiled function objects, passed the complete
tangent corpus and special checks at 10,000 with icx 2026.1.1, and with
GCC 15.2 ASan+UBSan. All runs checked four rounding modes, inexact and
errno. These runs are not part of the pipeline.
