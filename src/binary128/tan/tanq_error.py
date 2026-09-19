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
# a,b in [1/2,1), absolute component errors < 16*d.
# |a/b - A/B| <= (|a-A| + (a/b)*|b-B|)/B.
fast = 48/(F(1, 2)-16*d) + 1  # exact integer quotient truncation < one unit
assert fast < 128

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
