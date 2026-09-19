# Binary128 tangent: implementation and error argument

`tanq.c` is standalone. Its constants, Payne–Hanek reduction, residual
polynomials, and table reconstruction come from `master` at `671f2c73`
(the merged sine/cosine implementation). It needs no commit or linked
function from `chhe/sincosq`. Both components share one reduction and residual
evaluation. The quotient uses the **unrounded** component fractions.

This is a new implementation with analytic error bounds and MPFR regression
coverage. The corpus is not an exhaustive search for tangent's hardest cases.
The adaptive fallback avoids assuming that a fixed precision resolves every
rounding boundary. The argument below is not a machine-checked C proof.

## Edges and environment

* Preserve both zeros exactly, without raising any exception.
* Quiet NaNs, preserve their sign/payload, and raise invalid only for sNaNs.
* Infinities produce a quiet NaN and invalid. With errno enabled they set
  EDOM, matching `sinq`, `cosq`, and the binary128 worst-case harness. The
  root README's blanket statement that infinity leaves errno unchanged is
  inconsistent with these existing implementations and that harness.
* For `0 < |x| < 2^-57`, `0 < tan(|x|)-|x| < |x|^3/2` is below half an ulp.
  Nearest and rounding toward zero return x. Rounding away from zero adds
  one to its magnitude encoding. All these results are inexact.
* Every nonzero subnormal input raises underflow. At the largest subnormal,
  rounding away from zero returns the smallest normal, but the result
  rounded to 113 bits with an **unbounded exponent** is still subnormal.
  Consequently underflow (and optional ERANGE) is required there too.
* Non-tiny paths use integer arithmetic exclusively. The sign is resolved
  before rounding, and existing flags and the caller's mode are preserved.
  The fast and wide paths have explicit reduction guards keeping their
  quotient exponents far from either range limit. The adaptive rounder
  handles overflow and underflow without needing a global pole-separation
  theorem. No nonzero finite tangent is exact.

## Constants and inherited kernels

Run from the repository root:

```
python3 src/binary128/sin/sinq_tables.py > /tmp/tanq-tables.inc
python3 src/binary128/tan/tanq_error.py
python3 src/binary128/tan/tanq_cases.py > /tmp/tanq.wc
cmp /tmp/tanq.wc src/binary128/tan/tanq.wc
```

The first command reproduces all seven constant arrays/scalars in `tanq.c`
and checks the sine minimax error and tiered polynomial arithmetic with
exact rational inequalities. `tanq_error.py` verifies their identity with
`sinq.c` and the additional rational inequalities below. The accompanying
`../sin/sinq_error.md` explains the component errors. This implementation
uses master's 256-bit fast reconstruction, rather than the later 192-bit
reconstruction on `sincosq`; its placement operations are exact and its
omitted table limbs contribute less than 2^-256.

For the direct small-argument cosine, the shared reconstruction multiplies
the cosine correction by the saturated table value for one. Relative to the
old direct branch this adds at most four correction units, scaled by
`2^(2*et) <= 2^-16`. Its total remains below two normalized Q128 units.

## Narrow quotient and rounding gate

Write the represented normalized magnitudes as a and b in `[1/2,1)`, with
exponents carried separately, and let `d = 2^-128`. The inherited component
analysis gives `|a-A|, |b-B| < 16*d` for the exact scaled components A and B.
Since `a/b < 2`,

```
|a/b - A/B| <= (|a-A| + (a/b)*|b-B|) / B
             < 48*d / (1/2 - 16*d).
```

`divide_fast` truncates an exact normalized integer quotient, losing less
than one final Q128 unit. Normalizing the quotient divides by either one
or two, so its error is below 98 units. The gate uses a radius of **128**.
`div_digit` estimates a base-2^64 quotient digit from the high divisor limb,
then corrects it against the exact 192-bit product. The residual is less
than the divisor after each digit; two digits suffice for the fraction.

The 113-bit output has 15 guard bits. Reject a closed radius-128 window
around each midpoint (nearest mode) or grid point (directed modes). Outside
it the true magnitude and represented quotient round identically. This
also rejects ambiguous binade crossings. The final sign controls the
direction of the increment.

## Wide quotient and rounding gate

The wide reduction extracts 448 fraction bits. The missing window tail is
less than 2^-462, and extraction loses less than 2^-448. Reject `lz > 160`
**before** normalization; this replaces an implicit global separation
assumption with a checked condition. Multiplication by pi/2 and normalization
then give less than 2^-286 error in the normalized residual fraction,
including the constant and integer-product errors.

The inherited exact-rational analysis bounds the polynomial evaluation
error after multiplying by the squared residual by 2^-296. Including the
argument and product truncations and the possible one-bit normalization,
the floating sine's normalized error is less than 2^-284. In table
reconstruction, nonzero breakpoints have magnitude greater than 2^-8; the
two polynomial errors and all table, product, and placement errors remain
below 2^-284 after normalization. The `j=0` sine retains its floating form,
so cancellation does not destroy relative precision.

Applying the quotient inequality above with component errors 2^-284, plus
less than 2^-384 for exact restoring division, leaves the normalized
quotient error below **2^-270**. `tanq_error.py` checks these conservative
inequalities. The code adds/subtracts `2^114` integer units in the 384-bit
fraction and rounds both endpoints. It returns only if their rounded
encodings agree. An interval crossing the normalized fraction's binade
boundary is sent to the adaptive path instead of wrapping.

## Adaptive reduction

Use base-2^32 integers at scale `2^P`, where the requested working fraction
has K bits (initially 256) and

```
P = K + 32*ceil(max(e,0)/32) + 64.
```

Machin's identity gives `pi/2 = 8*atan(1/5) - 2*atan(1/239)`. Each arctangent
uses its alternating series. Successive integer powers are divided by
`q^2`; their error relative to the exact scaled power is below
`1/(1-1/q^2) < 2` units. Dividing by each odd denominator adds less than
one unit. Summing N terms therefore costs less than `3*N` units, and the
omitted tail after the first zero power is below two units. Combining the
radii produces an explicit integer interval `[p, p+width]` for scaled pi/2,
with `width < 32*P+64`.

The significand shifted to the wide scale is divided **exactly** by the
lower endpoint p, retaining the remainder and quotient modulo four.
The discrepancy from reduction by the true pi/2 is at most the full
quotient times `width`. This quotient (or its nearest-integer adjustment)
is below `2^(max(e,0)+2)`. The 64 extra bits make that discrepancy less than
one QK unit for every representable array size. Choosing the closest
multiple of p and truncating the residual to K fraction bits adds less
than one more unit. The represented magnitude t is at most pi/4. If its
sign is uncertain, the subsequent interval contains zero and forces more
precision.

This path recomputes pi as precision increases; it does not exhaust a
finite table of 2/pi digits. Integer precision, term counts, and allocation
are naturally limited by the C implementation's address space and stack;
no empirical fixed-precision termination cutoff is used.

## Adaptive sine/cosine and quotient

Let `u = floor(t*t*2^K)/2^K`. Its error is below one QK unit. Generate sine
and cosine terms using the usual Taylor recurrences, truncating each
product and division. Compare first to the exact series at represented u:
for a preceding term error below three units, the next error is bounded by

```
(3*(5/8) + 1)/2 + 1 < 3,
```

because `u < (pi/4)^2 < 5/8` and every denominator is at least two. When
both represented terms vanish, the remaining tails are below two units.
The sum of k term errors is below `3*k` units. Perturbing u by one unit
changes the complete sine/cosine series by less than one unit; perturbing
the angle by the reduction error changes either by less than two units.
The shared conservative radius `4*k + 16` covers all of these contributions.

If either component interval contains zero, double K and restart. Otherwise
monotonicity of positive division gives a rational enclosure from the
appropriate numerator and denominator endpoints. Odd quadrants invert the
quotient and change its sign. `round_ratio` performs exact rational binary
long division to get 113 significand bits, a rounding bit, and an exact
sticky remainder. It handles ties, subnormals and overflow explicitly.
Both endpoints must agree in their rounded encoding **and** range flags.

For a nonzero algebraic x, tan(x) is transcendental: an algebraic value
would make `exp(2*i*x)` algebraic via the tangent identity, contradicting
Hermite–Lindemann. It is therefore neither a pole nor a dyadic rounding
boundary. The enclosing intervals converge, so increasing precision
resolves the rounding decision without a hard-case completeness assumption.

## Reproducing checks

```
EXTRA_CFLAGS='-Werror -DCORE_MATH_CHECK_INEXACT -DCORE_MATH_SUPPORT_ERRNO' \
  OMP_NUM_THREADS=8 ./check.sh --worst tanq
EXTRA_CFLAGS='-Werror -DCORE_MATH_TESTS=1000000 -DCORE_MATH_CHECK_INEXACT -DCORE_MATH_SUPPORT_ERRNO' \
  CORE_MATH_TEST_SEED=1 OMP_NUM_THREADS=8 ./check.sh --special tanq
make -C src/binary128/tan check_internal CFLAGS='-O2 -frounding-math -DMPFR_WANT_FLOAT128'
src/binary128/tan/check_internal
```

`CORE_MATH_TANQ_FORCE_ACCURATE` skips the narrow path. Define that and
`CORE_MATH_TANQ_FORCE_ADAPTIVE` together to exercise adaptive refinement.
The latter is much slower for huge exponents. Setting
`CORE_MATH_TANQ_ADAPTIVE_BITS=64` starts below the target precision and
exercises repeated refinement. `check_internal` compares
100,000 narrow and 10,000 wide divisions against exact GMP quotients, checks
Machin and 3,000 paired Taylor intervals against MPFR, checks 400 signed
rational rounding/range cases, and samples 10,000 quotient error bounds.
These tests validate the implementation; sampling does not replace the
analytic inequalities above.

## Validation recorded on 2026-09-19

On an Intel Core i9-14900K, GCC 15.2 and Clang 21.1:

* The 4,959 generated regression inputs passed all four modes, normally and
  with the wide path forced, with inexact and errno enabled.
* All 196,978 input records from the combined existing `sinq.wc` and
  `cosq.wc` corpora passed as tangent inputs in all four modes, with inexact
  and errno enabled. `WORST_SYMMETRIC` also exercises negative inputs.
* The one-million special-test setting passed all four modes, plus its
  subnormal, exception, NaN, signed-zero, and environment-preservation checks.
  A final reproducible run used seed 1 and four OpenMP threads.
* Forced wide special checks passed at the 10,000 setting. Forced adaptive
  checks passed at 600, both at the default precision and starting at 64 bits
  to exercise repeated refinement. Inexact and errno were enabled.
* The independent internal arithmetic checks passed. Constant regeneration
  and exact-rational bound checks passed, and the corpus regenerated exactly.
* GCC address/undefined-behavior sanitizers passed forced-adaptive special
  checks at 60, including refinement from 64 bits; Clang undefined-behavior
  sanitizer passed normal special
  checks at 10,000. Clang needed GCC's quadmath include directory, supplied
  using `-idirafter "$(gcc -print-file-name=include)"` in CFLAGS.
* A `-march=core2` build passed special checks at 10,000 in all modes,
  exercising the conversion helpers without SSE4.1. Clang also compiled
  the forced-adaptive configuration with `-Wall -Wextra -Werror`.
* The CI selector schedules both tangent checks, and the CI dry build passed
  with `-march=x86-64-v2 -Wdouble-promotion`. As with the existing pipeline,
  this ISA flag is needed by the repository's quadmath capability probe.

The corpus contains regression inputs, including cancellation cases, rather
than a completed tangent-specific hard-to-round search. These sample counts
are not an exhaustive binary128 verification.

For an initial throughput baseline:

```
CORE_MATH_PERF_MODE=rdtsc PERF_COUNT=100000 PERF_REPEAT=100 \
  CORE_MATH_LAUNCHER='taskset -c 4' ./perf.sh tanq
```

Using the harness's default nearest-mode `[-pi,pi]` distribution and build
flags, the minimum of 20 trials was 404.494 TSC ticks/call for `cr_tanq`
(median 412.499), versus 1227.496 for glibc 2.43 (median 1228.844). This is
an initial libc comparison; no previous CORE-MATH `tanq` implementation
exists on the base branch. It does not measure latency, other input bands,
or directed-mode performance, and libc is not the correctness oracle.
