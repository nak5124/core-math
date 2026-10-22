# Binary128 tangent: implementation and error argument

`tanq.c` is standalone. Its constants, Payne–Hanek reduction, residual
polynomials, and table reconstruction come from `master` at `671f2c73`
(the merged sine/cosine implementation). It needs no commit or linked
function from `chhe/sincosq`. Both components share one reduction and residual
evaluation. The quotient uses the **unrounded** component fractions.

This is a new implementation with analytic error bounds and MPFR regression
coverage. The corpus is not an exhaustive search for tangent's hardest cases.
An input that neither gated path decides is printed and the process exits,
following `binary64/pow/pow.c`; there is no arbitrary-precision fallback.
The argument below is not a machine-checked C proof.

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
  Both paths have explicit reduction guards keeping their quotient
  exponents far from either range limit, so neither can overflow or
  underflow. An input that neither path decides is reported, not rounded
  (see below). No nonzero finite tangent is exact.

## Constants and inherited kernels

Run from the repository root:

```
python3 src/binary128/sin/sinq_tables.py > /tmp/tanq-tables.inc
python3 src/binary128/tan/tanq_error.py
python3 src/binary128/tan/tanq_cases.py > /tmp/tanq-generated.wc
```

The first command reproduces all seven constant arrays/scalars in `tanq.c`
and checks the sine minimax error and tiered polynomial arithmetic with
exact rational inequalities. `tanq_error.py` verifies their identity with
`sinq.c` and the additional rational inequalities below. The accompanying
`../sin/sinq_error.md` explains the component errors.

`tanq_cases.py` reproduces the first 4,961 records, not the whole current
corpus: `5f34efa3` appended 57,001 bacsel cases and `d8c2b9ea` appended
coverage inputs. Keep those additions when checking or regenerating the file.
The full corpus has 61,983 records.

The narrow path departs from master's reconstruction in two places, neither
of which changes a bound used below. First, `reduce` normalizes the product
with pi/2 by a one-bit select: both factors are normalized, so the high word
is at least 2^126 and `lzt <= 1`, as the inherited analysis already states.
The selected word equals what the general 256-bit shift produced.

Second, recombination uses a 192-bit frame `hi*2^-64 + lo*2^-192` instead of
256 bits, and forms sin A and cos A once for the breakpoint-plus-residual
angle A; the quadrant only chooses numerator and denominator and flips the
sign. With `-57 <= et <= -7` (direct arguments give `-56 <= et <= -8`):

* the 192-bit table word truncates the nearest 384-bit entry: less than one
  frame unit `2^-192`, including the entry's own rounding;
* the sine term `p*2^(et-128)` is placed by a left shift of `et+64`, in
  `[7,57]`. It is exact, and `p < 2^128` keeps it inside the frame;
* the cosine term rides `2^(2*et)`: two truncating right shifts of `-et`,
  each less than 64 bits and each losing less than one frame unit.

The frame therefore errs by less than three units of `2^-192`. With the
normalization shift `L <= 7` that is `3*2^-57` output units, replacing the
`2^-121` table row of the inherited budget (3.6 + 4.3 + 4 + 24/128 + 1, plus
this row), which remains below **14**. `tanq_error.py` checks that sum. The
floating `j=0` sine is unchanged, below **15.8**. All variable shifts are
narrower than 64 bits, which is the point of the frame.

For the direct small-argument cosine, the shared reconstruction multiplies
the cosine correction by the saturated table value for one. Relative to the
old direct branch this adds at most four correction units, scaled by
`2^(2*et) <= 2^-16`. Its total remains below two normalized Q128 units.

## Narrow quotient and rounding gate

Write the represented normalized magnitudes as a and b in `[1/2,1)`, with
exponents carried separately, `q = a/b`, and `d = 2^-128`. For the exact
scaled components A and B, `|a-A| < ea*d` and `|b-B| < eb*d`. Table
recombination gives less than 14 and the floating sine less than 15.8; at
most one of the two components is a floating sine. In general

```
|a/b - A/B| <= (ea + q*eb)*d / B,        B >= b - eb*d.
```

The two unfavourable extremes, `q` near 2 and `b` near 1/2, cannot both be
charged at full weight to the *normalized* fraction:

* `a < b`. The fraction is q itself, in `[1/2,1)`. Since `a >= 1/2`,
  `b = a/q >= 1/(2q)`, and the error is below
  `(ea + q*eb)*2q/(1 - 2q*eb*d)`. This increases with q; at `q = 1` it is
  `2*(ea+eb)`, less than 59.6 plus a negligible term.
* `a >= b`. Now `q < 2` and `b >= 1/2`, but the fraction is `q/2`, halving
  the error: below `(ea + 2*eb)` up to the same negligible factor, 45.6.

### Reciprocal division

`divide_fast` does not develop quotient digits. With `b = bh*2^64 + bl`, one
limb division gives `x0 = floor((2^128-1)/bh) - 2^64`, that is
`(2^64+x0)*bh = 2^128-1-rho` with `0 <= rho < bh`. Its dividend's top limb
`~bh` is below `bh` because b is normalized, so the quotient fits one limb.
For `Y0 = 1 + x0/2^64`, the residual `eps0 = 1 - Y0*b/2^128` satisfies

```
2^128*eps0 = 1 + rho - (2^64+x0)*bl/2^64,   in (-2^65, 2^64].
```

The code forms `t = b + floor(x0*b/2^64)` modulo 2^128 (both partial products
are exact) and reads `e = -t` as a signed integer: it is the ceiling of
`2^128*eps0`, below 2^66 in magnitude. One Newton step,
`X = x0*2^64 + e + floor(x0*e/2^64)`, approximates
`2^128*(Y0*(1+eps0) - 1) = 2^128*(2^128/b*(1-eps0^2) - 1)`. Then
`r + mhi_approx(r, X)` approximates `r*2^128/b` for the remainder `r < b`
left after the optional subtraction. In units of the quotient before the
final shift:

| Source | Interval |
| --- | ---: |
| Newton truncation `q*eps0^2`, with `eps0^2 < 2^-126` | `(-4, 0]` |
| Ceiling of e times `1+x0/2^64 < 2`, less one floor, scaled by `r/2^128 < 1` | `(-1, 2)` |
| `mhi_approx`, relative to the real product | `(-3, 0]` |

The quotient errs by more than -8 and less than 2. For `a >= b` the halving
shift truncates less than one further unit: the fraction error is in
`(-5, 1)`.

Three edges are closed explicitly. (1) The exact quotient is below 2^128 but
can exceed `2^128-2` when `r = b-1`; an overflowing sum saturates to
`2^128-1`, which only reduces the error. (2) X is an unsigned fraction of
the reciprocal, but for b within a few units of 2^128 its true value drops
below the error above and the computed X wraps. That can happen only when
`bh = 2^64-1`, where the true X is below `2^64+2`; a set top bit then
identifies the wrap, and the signed high product `mhi(r,X) - r` cancels the
leading r exactly. (3) For `a < b` the exact fraction exceeds 2^127; a
computed value below it is clamped to 2^127, again reducing the error, so
the frame reaching `round_fast` is always normalized.
(4) X does not reach 2^128 at the other end either. For `b = 2^127+k` with
`1 <= k < 2^127`, the true X is `2^128 - 2^129*k/(2^127+k) <= 2^128-2k`,
and the computed one is less than two above it. For `b = 2^127` every step is exact: `x0 = 2^64-1`,
`e = 2^63`, and `X = 2^128-1`.

The total narrow error is therefore below `59.6 + 8 <` **68** units, checked
with exact rationals by `tanq_error.py`. The gate radius is **72**. Each
unit of radius costs about 0.05 TSC ticks per call through the fallback, so
no effort is spent on shaving the constant further.

The 113-bit output has 15 guard bits. Reject a closed radius-72 window
around each midpoint (nearest mode) or grid point (directed modes). Outside
it the true magnitude and represented quotient round identically. This
also rejects ambiguous binade crossings. The final sign controls the
direction of the increment.

## Wide quotient and rounding gate

The wide reduction extracts 448 fraction bits. The missing window tail is
less than 2^-462, and extraction loses less than 2^-448. Reject `lz > 160`
**before** normalization; this replaces an implicit global separation
assumption with a checked condition. The wide path works on five
little-endian 64-bit limbs, that is 320 fraction bits. The residual, pi/2
and every table entry are the top 320 bits of the 384-bit constants shared
with `sinq`, so each costs less than one unit of `2^-320` beyond the entry's
own rounding. `mul_hi_320` keeps the product columns from bit 192 up; the
three dropped columns sum to less than `2^-62` of a unit, so each product is
at most one unit short and never over. Multiplication by pi/2 and
normalization then give less than 2^-286 error in the normalized residual
fraction, including the constant and integer-product errors.

The tiered polynomial evaluation of `../sin/sinq_tables.py`, with its top
tier at 320 instead of 384 bits, still errs by less than 2^-282 before and
2^-296 after multiplying by the squared residual: the 320-bit tier
contributes under `6*2^-320` per step, far below the 256-bit tier's
`2^-253` slack riding `u^2 < 2^-29`. Including the argument and product
truncations and the possible one-bit normalization, the floating sine's
normalized error is less than 2^-284. In table reconstruction, nonzero
breakpoints have magnitude greater than 2^-8; the two polynomial errors and
all table, product, and placement errors remain below 2^-284 after
normalization. The `j=0` sine retains its floating form, so cancellation
does not destroy relative precision.

Applying the quotient inequality above with component errors 2^-284, plus
less than 2^-320 for exact radix-2^64 division, leaves the normalized
quotient error below **2^-270**: the budget is set by the 448-bit
reduction and its `lz <= 160` guard, not by the limb count, which is why
the sixth limb could be dropped without loosening the bound.
`tanq_error.py` checks these conservative inequalities. The code
adds/subtracts `2^50` integer units in the 320-bit fraction and rounds both
endpoints. It returns only if their rounded encodings agree. An interval
crossing the normalized fraction's binade boundary is reported instead of
wrapping.

The wide divider generates five exact base-2^64 digits. Put `B = 2^64`, let
the five-limb divisor be `D = B^4*h + l` with `h >= B/2` and `l < B^4`, and
let `R < D` be the remainder before a zero digit is appended. The estimate
`qhat = min(B-1, floor(floor(R/B^3)/h))` is at least the exact digit
`floor(R*B/D)`, and `qhat*D - R*B <= qhat*l < B^5 <= 2*D`, so it exceeds
that digit by at most two. Clamping can only reduce the excess.
The code compares the complete six-limb product against `R*B` and, when
necessary, decrements the estimate and subtracts D from that product.
The resulting remainder is exact and in `[0,D)`. The product accumulation
fits in 128 bits: `(B-1)^2 + (B-1) < B^2`. Subtractions keep a 128-bit
borrow operand so that `UINT64_MAX+1` cannot lose its carry. Five such
steps, followed by the same optional one-bit normalization, give the
truncated 320-bit quotient that single-bit restoring steps would.

## Undecided inputs

`accurate` reports failure instead of a result when wide reduction refuses
more than 160 leading zeros, when the quotient's radius-`2^50` enclosure
leaves its normalized binade, or when its two endpoints round differently.
`cr_tanq` then prints

```
Unexpected worst-case found.
Please report to core-math@inria.fr:
Worst-case of tanq found: x=<sign>0x1.<28 hex digits>p<exponent>
```

and calls `exit(1)`, the convention of `binary64/pow/pow.c`. The input is
always normal here, and the printed form is the exact hexadecimal notation
of `tanq.wc`, so it can be appended to the corpus directly.

No input is known to reach this message: the worst-case corpus and the
special checks, forced through the wide path, are all decided there. That
is sampled evidence. The wide gate rejects only quotients whose normalized
fraction lies within `2^-270` of a rounding boundary, and the reduction
guard only arguments within `2^-160` of a multiple of pi/2, in units of
pi/2. Whether any binary128 input meets either condition is not settled
by this document.

## Reproducing checks

```
EXTRA_CFLAGS='-Werror -DCORE_MATH_CHECK_INEXACT -DCORE_MATH_SUPPORT_ERRNO' \
  OMP_NUM_THREADS=8 ./check.sh --worst tanq
EXTRA_CFLAGS='-Werror -DCORE_MATH_TESTS=1000000 -DCORE_MATH_CHECK_INEXACT -DCORE_MATH_SUPPORT_ERRNO' \
  CORE_MATH_TEST_SEED=1 OMP_NUM_THREADS=8 ./check.sh --special tanq
make -C src/binary128/tan check_internal CFLAGS='-O2 -frounding-math -DMPFR_WANT_FLOAT128'
src/binary128/tan/check_internal
```

To force the wide path, edit a scratch copy only: replace
`if (fast_pair(m, e, f, ex, &flip)) {` by
`if (0 && fast_pair(m, e, f, ex, &flip)) {` in `cr_tanq`, as the sine and
cosine notes do. `bench_compare.py --force-wide` applies the same edit to
its temporary copies. Special values and the tiny-input shortcut still run
first. This tests the wide arithmetic, not the narrow gate's rejection rate.

`check_internal` checks 321 leading-zero frames and 4,173 shifts against GMP,
including zero, limb boundaries and shifts at or beyond 320 bits, and
1,000,000 wide products against their one-unit bound. It compares
1,000,000 narrow quotients with their documented error interval and 100,000
exact wide quotients against GMP. The narrow operands include divisors at
and beside both ends of the binade and extreme limbs; the wide ones include
long carry/borrow chains. It reports the extreme narrow errors observed, and
samples 10,000 quotient error bounds against MPFR. These tests validate the implementation; sampling does not replace the
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
* Forced wide special checks passed at the 10,000 setting, with inexact
  and errno enabled.
* The independent internal arithmetic checks passed. Constant regeneration
  and exact-rational bound checks passed, and the corpus regenerated exactly.
* Clang undefined-behavior sanitizer passed normal special checks at
  10,000. Clang needed GCC's quadmath include directory, supplied
  using `-idirafter "$(gcc -print-file-name=include)"` in CFLAGS.
* A `-march=core2` build passed special checks at 10,000 in all modes,
  exercising the conversion helpers without SSE4.1.
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


## Removal of the adaptive fallback (2026-10-13)

The adaptive fixed-point interval tier (Machin pi/2, Taylor enclosures and
exact rational rounding) and its compile-time test macros were removed in
favour of the report-and-exit convention above, following review. The
narrow and wide paths, their gates and every constant are unchanged.
Validation of the resulting source (GCC 15.2, `-march=x86-64-v2`, inexact
and errno enabled, `-Werror`, all four rounding modes):

* `check.sh --worst`: 61,980 corpus records passed, normally and with the
  wide path forced in a scratch copy.
* `check.sh --special`: the deterministic suite plus 1,000,000 random inputs
  per mode (seed 1) passed; forced wide, 100,000 per mode (seed 2) passed.
* GCC address/undefined-behavior sanitizers passed the full corpus and
  100,000 special inputs per mode (seed 4) on the normal build.
* `check_internal`: 385 leading-zero frames, 10,010 shifts, 1,000,000
  narrow quotients (floor error -5 to 1), 100,000 exact wide quotients and
  10,000 quotient-bound samples passed; largest narrow error 4.0 of 72.
* `tanq_error.py` passed. Clang 21.1 and GCC `-march=core2` compiled the
  source with `-Wall -Wextra -Werror`, Clang also with the wide path forced.
* A scratch copy whose `accurate` always fails printed the three-line
  message for every tested normal input; each printed value parsed back to
  the identical binary128 encoding.

## Optimizations versus fb944b75

The initial implementation was committed as `fb944b75`. The follow-up keeps
the constants, the residual polynomials, the wide path's arithmetic and
every rounding rule, and changes:

1. wide division: six exact radix-2^64 digits instead of 383/384 single-bit
   restoring steps;
2. narrow recombination: the 192-bit frame, with sin A and cos A formed once;
3. narrow division: one limb division and a signed Newton step instead of
   quotient digits, under the error interval derived above;
4. the gate radius: 72 instead of 128, from the sharper quotient inequality.
   The fallback rate falls from 257/32768 to 145/32768;
5. limb divisions use the x86-64 `div` instruction directly when available.
   Their quotients provably fit one limb; the C expression is kept elsewhere.

On the same i9-14900K with GCC 15.2, CPU 4 affinity and identical input
files, the default-distribution comparison used 100,000 inputs, 20 repeats,
and five trials per version, alternating A/B and B/A pairs. Both versions
used `-O3 -march=native -frounding-math -fno-finite-math-only`. Values below
are median TSC ticks/call, not frequency-dependent core-cycle counts.

| Mode | Throughput before | Throughput after | Speedup | Latency before | Latency after |
| --- | ---: | ---: | ---: | ---: | ---: |
| rndn | 405.1 | 186.9 | 2.17x | 407.8 | 202.7 |
| rndz | 420.9 | 193.5 | 2.18x | 429.5 | 202.7 |
| rndu | 428.2 | 193.6 | 2.21x | 427.5 | 202.4 |
| rndd | 421.3 | 193.4 | 2.18x | 418.2 | 201.9 |

Nearest-mode throughput after each step, same inputs and session: exact
narrow digits with wide radix division 243.1; 192-bit frame 220.5; gate 64
with the digit divider 207.9; branch-free digit correction 198.3; direct
`div` 194.2; reciprocal division with gate 72, 187.4. With the gate disabled
the narrow path alone measured about 6.5 ticks less, which is what the
fallback still costs per call on average.

A separate three-trial run with 20,000 inputs and five repeats covered
random binary128 significands with exponents in `[-57,16383]` and
`[-57,-9]`, with both signs. Throughput speedups were 2.06-2.29x for the
full exponent band and 2.71-2.73x for small arguments; latency speedups were
1.98-2.03x and 2.50-2.60x. Forcing the wide path on the default distribution
costs about 3,150 ticks/call, against 21,057-21,988 before radix division.

Measured and not adopted: replacing both series by a single tangent series,
`(sin A + cos A*tan t)/(cos A - sin A*tan t)`. A timing-only stub that drops
the cosine series and two recombination products while lengthening the
remaining polynomial by four products gained about 8 ticks. That does not
justify new coefficients and a new error analysis at present.

Reproduce the matched baseline/candidate measurements with:

```
python3 src/binary128/tan/bench_compare.py --baseline fb944b75 --cpu 4
python3 src/binary128/tan/bench_compare.py --baseline fb944b75 --cpu 4 \
  --bands default --count 100000 --repeat 20 --trials 5
python3 src/binary128/tan/bench_compare.py --baseline fb944b75 --cpu 4 \
  --bands default --count 10000 --repeat 2 --trials 3 --force-wide
```

The script builds the committed baseline and working source with identical
flags, uses the repository's throughput/latency harness, and prints all raw
samples alongside medians. It is a performance comparison, not an oracle.

Validation after these changes (GCC 15.2, `-march=x86-64-v2`, inexact and
errno enabled throughout, all four modes unless stated):

* 1,000,000 narrow quotients stayed within the documented interval against
  GMP (observed error relative to the floor: -5 to 1), and 100,000 wide
  quotients were exact. The internal quotient-bound checks passed; the
  largest sampled narrow error was 4.4 of 72 units.
* The 4,959-input regression corpus passed normally and with the wide path
  forced.
* All 196,978 inherited sine/cosine input records passed as tangent inputs.
* Special checks passed at the 50-million setting (seed 3) and at one
  million (seed 1) normally, and at 100,000 (seed 2) with the wide path
  forced.
* GCC address/undefined-behavior sanitizers passed normal special checks at
  100,000 (seed 4).
* Clang 21.1 compiled the normal and forced-wide configurations with
  `-Wall -Wextra -Werror`; GCC compiled `-march=core2` likewise. These two were compile checks only in this round.
* The exact-rational bound check and `git diff --check` passed.

As with the initial implementation, the MPFR samples are regression evidence
and the hard-case corpus is not exhaustive.

## Five-limb accurate path (2026-10-22)

Following review, the accurate path was narrowed from three 128-bit
limbs (384 bits) to five 64-bit limbs (320 bits), the width `log2q` and
`log10q` use. The constants stay the 384-bit arrays shared with `sinq`;
`top_320` reads their top five limbs. The adders are plain carry chains on
the limbs, the product is column-wise, the divider produces five digits,
and the rounding gate keeps its radius of `2^-270`, now `2^50` units. The
rational bookkeeping above was redone for this width and `tanq_error.py`
checks it; the quotient bound is unchanged because the reduction, not the
limb count, dominates the budget. Validation is recorded in the merge
request.

## Equivalent limb arithmetic (compiler review, 2026-09-21)

The follow-up compiler changes preserve every intermediate integer result:

* `shr_round` distinguishes shifts below 64, equal to 64, and above 64.
  The boundary case folds bit 63 of the low limb; all evaluated 64-bit
  shifts are in `[0,63]`. The result is the same shifted word with the
  adjacent discarded bit ORed into its low bit.
* The 192-bit shift and frame normalization explicitly funnel 64-bit limbs.
  The existing contracts (`0 < shift < 64`, nonzero frame top limb) ensure
  valid shifts. The split right funnel handles zero normalization shift.
* In each reduction-product column, with `B=2^64`, `A=low*current`, and
  `C=high*previous+carry`, the output is `(A+C) mod B` and the next carry is
  `floor(A/B)+floor(C/B)+floor(((A mod B)+(C mod B))/B)`.
  This exactly reproduces the original column even when `A+C >= B^2`.
  The low-limb sum is below `2*B`; `high < 2^49` bounds the incoming carry
  below `B+2^50` and `C` below `2^114`, so the 128-bit additions cannot
  overflow. Explicit first/last columns preserve all `W+2` product limbs
  and the zero sentinel for both the six- and ten-limb windows.

The reciprocal/Newton divider, clamping rules and error radii are
unchanged. See [the compiler review](compiler_review.md)
for measurements, verification coverage and reproduction commands.
