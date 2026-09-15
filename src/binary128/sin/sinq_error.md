# Binary128 sine/cosine fast-path error

This note applies to the identical fixed-point pipelines in `sinq.c` and
`../cos/cosq.c`. It bounds the **unrounded fast result**; it does not prove
that the 384-bit fallback resolves every binary128 rounding boundary.
`python3 sinq_tables.py > tables.inc` checks the rational inequalities below
and the minimax approximation bound. The generated constants must match both
C files. MPFR samples are regression checks, not a proof over the input domain.

## Units and elementary operations

Put `d = 2^-128`, `U = 0.00003765`, and `R = 1/6 + 0.000001`.
For a normalized residual, write `T = a*2^et`, where
`a = t1*d` is in `[1/2,1)` and `et <= -7`. Both `T^2` and the Q128 value
`u` are below `U`. The slack between `(pi/512)^2` and `U` exceeds the
reduction errors below. Direct evaluation has `et <= -8` and no reduction
error. Throughout the following kernel bounds, errors are in units of `d`.

* `wmul` is exact. Integer addition/subtraction and placement into the
  256-bit frame are exact when their results are in range.
* `mhi_approx(X,Y)` omits the low-by-low product and the fractional parts of
  the two cross products. Consequently its error relative to the **real**
  product `X*Y/2^128` is in `[0,3)`, downward. Relative to the integer high
  half it is at most two units downward. These are different bounds.
* `shr_round` discards low bits and ORs the next bit into the retained low
  bit. Its error is less than one unit; it is not round-to-nearest.
* Truncating a positive normalized product to its top word loses less than
  one unit. Normalizing `sin_frac` can shift left by at most one bit because
  its value is greater than `0.49*2^128`.
* The generated pi/2 constant is nearest at scale `2^127`. The sine/cosine
  table is nearest at scale `2^384`, except for saturated `cos(0)`. Its top
  128 and 256 bits have absolute errors less than `d` and `2^-256`,
  respectively. Taylor coefficient top words differ from their exact
  factorial reciprocals by less than `d`; the coefficient `1/2` is exact.

All intermediate polynomial partial sums are positive and below one. This
follows from the alternating coefficient magnitudes on `[0,U]`; their
separations are much larger than the error allowances below. There is no
unsigned wrap in these evaluations.

## Argument reduction: less than 3.6 residual units

The six-limb window omits less than `2^-206` from `x*2/pi`, and extracting
the 192-bit fraction loses less than `2^-192`. The omitted leading limbs
contribute multiples of four and do not affect the quadrant. Thus the
retained signed residual differs from the exact residual at the selected
breakpoint by less than `(1 + 2^-14)*2^-192`. Even if this perturbation changes
which neighboring breakpoint is selected, the same bound holds relative to
the **selected** breakpoint, and the true residual is at most that much
outside `[-1/256,1/256]`.

For accepted reductions, `lz <= 56`; multiplying the normalized residual by
pi/2 gives `lzt <= 1`, and `et = 1-lz-lzt`. In units of `2^(et-128)`, the
four sources of residual error are bounded by:

| Source | Bound |
| --- | ---: |
| Window tail and 192-bit extraction | `(1+2^-14)*2^(lz+lzt-65) <= (1+2^-14)/256` |
| Dropping residual bits below the normalized 128-bit word | `(pi/2)*2^(lzt-1) < 11/7` |
| Rounding the pi/2 constant | `< 1` |
| Truncating the normalized product | `< 1` |

Their sum is less than 3.6. The retained residual has magnitude at least
`2^-57`, so the extraction error cannot change its sign. Quadrant symmetry
therefore preserves the sign of the result, including the `j=0` sine case.

## Sine kernel: less than 4.3 residual units

Here the target is `sin(T)`, with `T` the represented residual. Reduction
error is added later. `squares` gives an error in `u` of less than
`eu*d`, where `eu = 1 + 3/16384`: the square's three-unit error is scaled
by `2^(2*et)`, then `shr_round` contributes one unit. The error in `v`
relative to the square of the **represented** `u` is less than `3*d`.

Let `c[k] = SIN_FAST[k]*d` and
`P(u) = c[0]-u*c[1]+u^2*c[2]-u^3*c[3]+u^4*c[4]`.
The generator proves, using rational arithmetic and a Taylor remainder,

```
|P(z) - (1-sin(sqrt(z))/sqrt(z))/z| < 2^-114,  0 <= z <= U.
```

The value at zero is interpreted by continuity. Coefficient bounds give
`|P| < R`, `|(u*P)'| < R`, `c[2]+2*U^2*c[4] < 1/1000`, and
`c[3] < 1/100000`. For the even and odd chains in `sin_frac`, define

```
ea = 3 + 3*U^2 + 3/1000
eb = 3 + 3/100000
ep = ea + 3 + U*eb
```

These bound the errors of `A`, `B`, and `A-mhi(u,B)`, respectively, in Q128
units: each approximate multiply contributes three units, inner errors are
multiplied by `u` or `v`, and the derivative bounds account for `v`'s error.
The final two multiplies, the error in `u`, and the minimax remainder give

```
|sin_frac*d - sin(T)/2^et|
  < d * (3 + 3*R + U*ep + eu*R + U*2^14)
  < 4.3*d.
```

The generator evaluates the parenthesized bound with exact fractions
(approximately 4.283785). For direct sine there is no reduction error.
For reduced sine without table recombination the error is therefore less
than `3.6+4.3` units before final normalization, and less than **15.8** units
after the possible one-bit left shift.

## Cosine correction: less than 20 squared-residual units

The target of `cos_corr` is `(1-cos(T))/2^(2*et)`. The even and odd cosine
chains each have error less than `4.01*d`, including coefficient truncation
and the error in `v`. For example, bounds for the even and odd chains are
`4+4*U^2+U^4+3/1000` and `4+U^2+3/100000`.
The derivative of the Taylor polynomial for `(1-cos(sqrt(u)))/u` has
magnitude below `1/20` on `[0,U]`.

Put `ew = 4.01+3+U*4.01`. The two outer multiplies contribute three units
each; the error in `u1` contributes at most `3/2` because the correction
factor is at most `1/2`. The error in `u` contributes `eu/20`. Finally the
first omitted cosine term, `T^14/14!`, contributes at most
`U^6*2^128/14!` in these units. Hence

```
|cos_corr*d - (1-cos(T))/2^(2*et)|
  < d * (9/2 + 3 + U*ew + eu/20 + U^6*2^128/14!)
  < 20*d.
```

The parenthesized bound is approximately 18.668128. In direct cosine,
subtracting this correction from the saturated 256-bit representation of
one and taking the top word costs another `1+2^-128` Q128 units. Since
`et <= -8`, the total is less than `1+20/65536+2^-128 < 2` output units.

## Table recombination: less than 14 output units

Apart from the floating sine case already handled, the exact output
magnitude is at least `sin(pi/512)` minus the negligible reduction error,
or at least `cos(pi/512)` for `j=0`. The unnormalized recombination error
from the bounds above is below `2^-120`, so the computed frame remains
positive and greater than `2^-8`. Thus its normalization shift `L` is at
most 7. Since `et <= -7`, `2^(et+L) <= 1` and `2^(2*et+L) <= 1/128`.
The frame also stays below one: for nonzero breakpoints the angle is away
from either endpoint, and for `j=0` the correction is positive.

In final output units `2^(-L-128)`:

| Source | Bound |
| --- | ---: |
| Residual error (the derivative of sine/cosine has magnitude at most one) | `3.6` |
| Sine evaluation | `4.3` |
| Sine table top-word error and approximate multiply | `1+3` |
| Cosine evaluation, table top-word error and approximate multiply | `(20+1+3)/128` |
| Truncating the normalized result | `1` |
| Initial 256-bit table error | `2^-121` |

The sum is less than 14. This also justifies the normalization bound used
above without assuming the approximate frame equals the true output.

## Why the gate works in all four rounding modes

For every accepted fast reduction, the returned magnitude `Y` differs from
the exact magnitude by less than `15.8*h`, where `h = 2^(e2-128)`.
The result has 15 guard bits below its 113-bit significand. In that frame,
representable grid points repeat every 32768 units and nearest-rounding
midpoints lie at offset 16384.

`round_fast` rejects a **closed** window of radius 16 around the relevant
boundary, using modular arithmetic: midpoint for nearest, grid point for
all directed modes. Outside that window the error interval cannot cross a
rounding boundary. The retained guard and sign therefore determine the
same rounded result as the exact magnitude. At a binade endpoint the
nearest-rounding cell still contains this tiny interval; directed modes
reject the endpoint window. All results reaching this path are normal;
the separate tiny-input path handles possible subnormal results.

This analytic bound is deliberately more conservative than sampled MPFR
errors. It establishes the fast gate's error budget, not a fourfold margin
and not an exhaustive correct-rounding proof for the accurate fallback.
