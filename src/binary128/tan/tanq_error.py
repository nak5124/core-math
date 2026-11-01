#!/usr/bin/env python3
"""Exact rational bookkeeping for tanq_error.md; no sampled error estimates."""
from fractions import Fraction as F
from pathlib import Path
import re

# These constants and bounds come from the unchanged master sinq kernel.
root = Path(__file__).resolve().parent
sine = (root.parent/'sin/sinq.c').read_text()
tan = (root/'tanq.c').read_text()
for name in ('FRAC_2_PI', 'PIO2_128', 'PIO2_384', 'SINCOS', 'SIN_FAST', 'SIN_COEF', 'COS_COEF'):
    pattern = rf'static const [^;]*\b{name}\b[^;]*;'
    assert re.search(pattern, sine).group() == re.search(pattern, tan).group(), name

d = F(1, 2**128)
# a,b in [1/2,1) with q = a/b; |a/b - A/B| <= (ea + q*eb)*d/B, B >= b - eb*d.
# Table-recombined components err by < 14*d, including 3*2^-57 for the
# 192-bit frame; the floating sine, at most one of the two, by < 15.8*d.
table_fast = F(36, 10) + F(43, 10) + 4 + F(24, 128) + 1 + F(3, 2**57)
assert table_fast < 14
gate = 72
# Reciprocal division. With bh the divisor's high limb, bl < 2^64 its low one,
# (2^64+x0)*bh = 2^128-1-rho and 0 <= rho < bh, the scaled residual
# 2^128*eps0 = 1 + rho - (2^64+x0)*bl/2^64 lies in (-2^65, 2^64].
assert 1 + (2**64 - 2) <= 2**64            # rho <= bh-1 <= 2^64-2
assert F((2**65 - 1)*(2**64 - 1), 2**64) < 2**65  # 2^64+x0 <= 2^65-1
eps0_sq = F(2**65, 2**128)**2
newton = 2**128*eps0_sq  # q*eps0^2 with q < 2^128: downward
assert newton == 4
# Ceiling of the residual, times 1+x0/2^64 < 2, minus one product floor:
# the reciprocal errs in (-1, 2), scaled by r/2^128 < 1. The final
# approximate product is in (-3, 0].
division_low, division_high = newton + 1 + 3, 2
assert max(division_low, division_high) == 8
for ea, eb in ((F(14), F(158, 10)), (F(158, 10), F(14))):
    # a < b: the fraction is q < 1 itself, and a >= 1/2 gives b >= 1/(2q).
    # (ea + q*eb)*2q/(1 - 2q*eb*d) increases with q; evaluate at q = 1.
    below = 2*(ea + eb)/(1 - 2*eb*d) + 8
    # a >= b: q < 2 and b >= 1/2; the fraction q/2 halves both errors, and
    # the halving shift truncates less than one more unit.
    above = (ea + 2*eb)/(2*(F(1, 2) - eb*d)) + 4 + 1
    assert max(below, above) < 68 < gate
assert f'#define ZIV_GATE {gate}u' in tan
# The wide path keeps five 64-bit limbs: W = 320 fraction bits.
W = 320
# Ten-limb window tail and 448-bit extraction, with lz <= 160.
reduction = (1 + F(1, 2**14))*F(1, 2**448)*2**160*F(11, 7)
reduction += F(16, 2**W)  # pi/2 and residual truncation, normalization, product
assert reduction < F(1, 2**286)
# Tiered polynomial error as in ../sin/sinq_tables.py, with the top tier at
# W bits: truncated coefficient, product at most one unit short, under six
# units per step; then the caller's multiply by u (u < 2^-14.68).
u = F(3765, 100000000)
poly = F(2, 2**128)
for k in range(16, -1, -1):
    poly = u*poly + F(6, 2**(128 if k >= 11 else 256 if k >= 2 else W))
assert poly < F(1, 2**282)
poly = u*poly + F(5, 2**W)
assert poly < F(1, 2**296)
# Floating sine can normalize left by one bit.
floating = 2*(reduction + poly + F(32, 2**W))
# Nonzero table breakpoints have normalized magnitude >= 2^-8. Sum both
# polynomial contributions and all table/product/placement truncations.
table = 256*(2*poly + F(1, 2**440) + F(64, 2**W))
assert floating < F(1, 2**284)
assert table < F(1, 2**284)
wide = 3*F(1, 2**284)/(F(1, 2)-F(1, 2**284)) + F(1, 2**W)
assert wide < F(1, 2**270)
# The code's radius: 2^-270 in units of 2^-W.
assert 2**(W - 270) == 2**50 and 'u320 radius = {{(u64)1 << 50, 0, 0, 0, 0}};' in tan

print('Shared constants match; narrow and wide (320-bit) rational bounds passed.')
