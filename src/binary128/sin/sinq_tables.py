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
from fractions import Fraction
import sys
from mpmath import mp, mpf, pi, sin, cos, floor, nint

LIMBS = 264    # 64-bit limbs of 2/pi: a 10-limb window starting at (e-114)/64 reaches e = 16383
SECTORS = 128  # breakpoints j*pi/256 for j = 0..127
TERMS = 18     # Taylor terms: u^18/37! < 2^-407 for u < 2^-14.69, past the 384-bit frame
FRAME = 384    # fixed-point frame: value * 2^FRAME in three little-endian 128-bit limbs

def check_fast_error(coefficients):
    """Exact rational bookkeeping for sinq_error.md, in Q128 units.

    This checks the stated analytic inequalities, not an exhaustive execution
    of the C code. The reduction and recombination arguments are in the note.
    """
    F = Fraction
    u = F("3.765e-5")
    c = [F(x, 1 << 128) for x in coefficients]
    r = F(1, 6) + F(1, 1000000)
    assert sum(c[k] * u**k for k in range(5)) < r
    assert sum((k+1) * c[k] * u**k for k in range(5)) < r
    assert c[2] + 2*u*u*c[4] < F(1, 1000)
    assert c[3] < F(1, 100000)
    eu = 1 + F(3, 1 << 14)
    ea = 3 + 3*u*u + F(3, 1000)
    eb = 3 + F(3, 100000)
    ep = ea + 3 + u*eb
    sine = 3 + 3*r + u*ep + eu*r + u*(1 << 14)
    assert sine < F("4.3")

    # The cosine's E/O chains include coefficient truncation. Their errors
    # are each < 4.01 units; the omitted term is theta^14/14!.
    cc = [F(1, factorial(2*k+2)) for k in range(6)]
    assert cc[3] + 2*u*u*cc[5] < F(1, 1000)
    assert cc[4] < F(1, 100000)
    assert 1 + 3 + 4*u*u + u**4 + F(3, 1000) < F("4.01")
    assert 1 + 3 + u*u + F(3, 100000) < F("4.01")
    assert sum(k*cc[k]*u**(k-1) for k in range(1, 6)) < F(1, 20)
    ew = F("4.01") + 3 + u*F("4.01")
    cosine = F(9, 2) + 3 + u*ew + eu/20 + u**6*(1 << 128)/factorial(14)
    assert cosine < 20

    # pi/2 < 11/7; lz <= 56 and lzt <= 1 in reduce().
    reduction = (1 + F(1, 1 << 14))/256 + F(11, 7) + 1 + 1
    assert reduction < F("3.6")
    floating_sine = 2*(F("3.6") + F("4.3"))
    table = F("3.6") + F("4.3") + 4 + F(24, 128) + 1 + F(1, 1 << 121)
    direct_cosine = 1 + F(20, 1 << 16) + F(1, 1 << 128)
    assert table < 14
    assert direct_cosine < 2
    assert floating_sine < 16
    print("Fast-path error < 15.8 guard units (analytic bookkeeping; "
          "see sinq_error.md)", file=sys.stderr)


def check_tiered_error():
    """Bound alternating() against the exact Taylor polynomial at its input u.

    Check that all computed partial sums stay in [0,1). At width w,
    coefficient and argument truncation each cost less than one unit,
    and the approximate
    product costs less than three (128/256 bits) or five (384 bits). At
    384 bits there is no argument truncation and coefficient rounding costs
    at most half a unit. Six units therefore bound every local step.
    This bounds evaluation only, not argument reduction or final rounding.
    """
    u = Fraction(3765, 100000000)
    for offset in (2, 3):  # cosine and sine factorials
        coefficients = [((1 << 385) + factorial(2*k + offset))
                        // (2 * factorial(2*k + offset)) for k in range(18)]
        upper = Fraction(coefficients[17] >> 256, 1 << 128)
        for k in range(16, -1, -1):
            width = 128 if k >= 11 else 256 if k >= 2 else 384
            c = Fraction(coefficients[k] >> (384 - width), 1 << width)
            assert 0 <= c - u * upper <= c < 1
            upper = c  # every approximate product is nonnegative
    error = Fraction(2, 1 << 128)  # initial coefficient, including Q384 rounding
    for k in range(16, -1, -1):
        width = 128 if k >= 11 else 256 if k >= 2 else 384
        error = u * error + Fraction(6, 1 << width)
    assert error < Fraction(1, 1 << 282)
    assert u * error + Fraction(5, 1 << 384) < Fraction(1, 1 << 296)
    print("Accurate polynomial evaluation error after u multiplication < 2^-296"
          " (exact rational bound)", file=sys.stderr)


check_tiered_error()

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
        check_fast_error(coefficients)
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
