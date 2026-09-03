#!/usr/bin/env python3
"""Fixed-point tables for the CORE-MATH binary128 sinq/cosq.

Usage:  python3 sinq_tables.py > tables.inc

Prints C source: 2/pi as 64-bit limbs for the Payne-Hanek reduction, pi/2 for
the reconstruction, sin/cos of the breakpoints j*pi/256, and the Taylor
coefficients 1/(2k+3)! and 1/(2k+2)! of sin(t)/t and 1 - cos(t).  The mathematical constants are evaluated with mpmath (transcendentals) or exact
integer arithmetic (the factorial reciprocals). SIN_FAST is an independent
rminimax fit; its command and validation are below. Nothing is transcribed
from another library.
The including C file must typedef u64 (uint64_t) and u128 (unsigned __int128).
"""
from math import factorial
from mpmath import mp, mpf, pi, sin, cos, floor, nint

LIMBS = 264    # 64-bit limbs of 2/pi: a 10-limb window starting at (e-114)/64 reaches e = 16383
SECTORS = 128  # breakpoints j*pi/256 for j = 0..127
TERMS = 18     # Taylor terms: u^18/37! < 2^-407 for u < 2^-14.69, past the 384-bit frame
FRAME = 384    # fixed-point frame: value * 2^FRAME in three little-endian 128-bit limbs

mp.prec = 64 * LIMBS + 304  # 17200 bits: every bit of the 2/pi limbs plus a margin


def scaled(x, bits):
    """x * 2^bits rounded to nearest, as a Python int."""
    return int(nint(x * mpf(2) ** bits))


def ratio(num, den, bits):
    """num/den * 2^bits rounded to nearest, exactly (ties away from zero)."""
    return ((num << bits + 1) + den) // (2 * den)


def u128(v):
    """One 128-bit limb as a U128(hi, lo) macro call."""
    assert 0 <= v < 1 << 128, v
    return f"U128(0x{v >> 64:016x}, 0x{v & (1 << 64) - 1:016x})"


def limbs(v, n=3):
    """v as n little-endian 128-bit limbs in braces."""
    assert 0 <= v < 1 << 128 * n, v
    return "{" + ", ".join(u128(v >> 128 * i & (1 << 128) - 1) for i in range(n)) + "}"


print("#define U128(hi, lo) (((u128)(hi) << 64) | (u64)(lo))")

# floor(2/pi * 2^(64*LIMBS)): bit 64i+1 ..= 64i+64 after the binary point is limb i.
two_over_pi = int(floor(mpf(2) ** (64 * LIMBS + 1) / pi))
words = [two_over_pi >> 64 * i & (1 << 64) - 1 for i in reversed(range(LIMBS))]
print("// 2/pi as 64-bit limbs: FRAC_2_PI[i] holds bits 64i+1 ..= 64i+64 after the binary point,")
print("// most significant limb first (264 limbs: a window of 10 limbs must reach exponent 16383).")
print(f"static const u64 FRAC_2_PI[{LIMBS}] = {{")
for i in range(0, LIMBS, 4):
    print("  " + ", ".join(f"0x{w:016x}" for w in words[i:i + 4]) + ",")
print("};")

print("// pi/2 * 2^127, rounded to nearest.")
print(f"static const u128 PIO2_128 = {u128(scaled(pi / 2, 127))};")
print("// pi/2 * 2^383, rounded to nearest, little-endian 128-bit limbs.")
print(f"static const u128 PIO2_384[3] = {limbs(scaled(pi / 2, 383))};")

print(f"// [sin(j*pi/{2 * SECTORS}), cos(j*pi/{2 * SECTORS})] * 2^{FRAME}, rounded to nearest,"
      " little-endian 128-bit limbs;")
print(f"// cos(0) saturates to 2^{FRAME} - 1.")
print(f"static const u128 SINCOS[{SECTORS}][2][3] = {{")
for j in range(SECTORS):
    theta = j * pi / (2 * SECTORS)
    s, c = scaled(sin(theta), FRAME), scaled(cos(theta), FRAME)
    assert c < 1 << FRAME or j == 0, j
    print(f"  {{{limbs(s)}, {limbs(min(c, (1 << FRAME) - 1))}}},")
print("};")

# Independent rminimax fit, with 140-bit coefficients before Q128 rounding:
# ratapprox --function="(1-sin(sqrt(x))/sqrt(x))/x" --dom="[8.7e-19,3.765e-5]" --type=[4,0] --numF="[140]" --prec=400 --dispCoeff=hex
SIN_FAST_HEX = [
    "0x2.aaaaaaaaaaaaaaaaaaaaaaaaaaa8b20314cp-4",
    "-0x2.222222222222222222221fa2fd12429be5p-8",
    "0xd.00d00d00d00d00c7f4be42e1a6409d95ep-16",
    "-0x2.e3bc74aad85378e3e9d9c3058a37d5da418p-20",
    "0x6.b99115ea3fd77eea5743ff7be5de64bc63p-28",
]


def fast_sine_coefficients():
    """Round the fit to Q128 and check its error and positive partial sums."""
    with mp.workprec(400):
        coefficients = []
        for text in SIN_FAST_HEX:
            mantissa, exponent = text.lstrip("-")[2:].split("p")
            whole, fraction = mantissa.split(".")
            value = mpf(int(whole + fraction, 16)) * mpf(2)**(int(exponent) - 4*len(fraction))
            coefficients.append(int(nint(value * mpf(2)**128)))
        c = [mpf(v) / mpf(2)**128 for v in coefficients]
        end = (pi / 512)**2
        worst = mpf(0)
        for i in range(4097):
            u = end * i / 4096
            p = c[-1]
            for a in reversed(c[:-1]):
                p = a - u*p
                assert p > 0
            exact = (1 - sin(mp.sqrt(u))/mp.sqrt(u))/u if u else mpf(1)/6
            worst = max(worst, abs(exact - p))
        assert worst < mpf(2)**-114, worst
        # Certify between grid points using exact rational arithmetic. Compare
        # with Taylor through u^5; its alternating remainder is <= u^6/15!.
        from fractions import Fraction as F
        from math import factorial as fact
        end_bound, steps = F("3.765e-5"), 8192
        delta = [(-1)**k * (F(1, fact(2*k+3)) -
                 (F(coefficients[k], 1 << 128) if k < 5 else 0)) for k in range(6)]
        derivative = sum(k*abs(delta[k])*end_bound**(k-1) for k in range(1, 6))
        grid = F(0)
        for i in range(steps + 1):
            u = end_bound*i/steps
            p = delta[-1]
            for a in reversed(delta[:-1]):
                p = a + u*p
            grid = max(grid, abs(p))
        bound = grid + derivative*end_bound/(2*steps) + end_bound**6/fact(15)
        assert bound < F(1, 1 << 114), bound
        import sys
        print(f"SIN_FAST rational error bound: 2^{float(mp.log(mpf(bound.numerator)/bound.denominator, 2)):.4f}", file=sys.stderr)
        print(f"SIN_FAST sampled Q128 error: 2^{float(mp.log(worst, 2)):.4f}", file=sys.stderr)
        return coefficients


print("// Magnitudes of the degree-4 minimax for (1-sin(sqrt(u))/sqrt(u))/u, Q128, alternating signs.")
print('// ratapprox --function="(1-sin(sqrt(x))/sqrt(x))/x" --dom="[8.7e-19,3.765e-5]" --type=[4,0] --numF="[140]" --prec=400 --dispCoeff=hex')
print("// Certified Q128 error < 2^-114 on [0, 3.765e-5], containing (pi/512)^2.")
print("static const u128 SIN_FAST[5] = {")
for c in fast_sine_coefficients():
    print(f"  {u128(c)},")
print("};")

print(f"// 1/(2k+3)! * 2^{FRAME} for k = 0..{TERMS - 1}:"
      " sin(t)/t = 1 - sum_k (-1)^k SIN_COEF[k] (t^2)^(k+1).")
print(f"static const u128 SIN_COEF[{TERMS}][3] = {{")
for k in range(TERMS):
    print(f"  {limbs(ratio(1, factorial(2 * k + 3), FRAME))},")
print("};")

print(f"// 1/(2k+2)! * 2^{FRAME} for k = 0..{TERMS - 1}:"
      " 1 - cos(t) = sum_k (-1)^k COS_COEF[k] (t^2)^(k+1).")
print(f"static const u128 COS_COEF[{TERMS}][3] = {{")
for k in range(TERMS):
    print(f"  {limbs(ratio(1, factorial(2 * k + 2), FRAME))},")
print("};")
