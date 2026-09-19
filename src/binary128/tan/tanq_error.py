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
# Ten-limb window tail and 448-bit extraction, with lz <= 160.
reduction = (1 + F(1, 2**14))*F(1, 2**448)*2**160*F(11, 7)
reduction += F(16, 2**384)  # pi/2, normalization and product truncations
assert reduction < F(1, 2**286)
# Existing tiered polynomial error, after multiplying by u.
poly = F(1, 2**296)
# Floating sine can normalize left by one bit.
floating = 2*(reduction + poly + F(32, 2**384))
# Nonzero table breakpoints have normalized magnitude >= 2^-8. Sum both
# polynomial contributions and all table/product/placement truncations.
table = 256*(2*poly + F(1, 2**440) + F(64, 2**384))
assert floating < F(1, 2**284)
assert table < F(1, 2**284)
wide = 3*F(1, 2**284)/(F(1, 2)-F(1, 2**284)) + F(1, 2**384)
assert wide < F(1, 2**270)

# Adaptive term recurrence: u < (pi/4)^2 < 5/8 and denom >= 2.
# Product truncation followed by division costs < 1+1/denom.
assert (F(3)*F(5, 8)+1)/2 + 1 < 3
# When a represented term first vanishes, the exact term is <3 units;
# subsequent terms shrink by at least 5/16.
assert 3*F(5, 16)/(1-F(5, 16)) < 2
print('Shared constants match; narrow, wide and adaptive rational bounds passed.')
