/* Correctly-rounded sine function for binary64 value.

Copyright (c) 2022-2023 Paul Zimmermann and Tom Hubrecht

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include <fenv.h>

// #define TRACE 0x1.1f1563d205f5fp+1023
// #define TRACE 0x1.005023d32fee5p+1
#define TRACE 0x1.0192deafac33fp+0

/******************** code copied from dint.h and pow.[ch] *******************/

typedef unsigned __int128 u128;

typedef union {
  struct {
    u128 r;
    int64_t _ex;
    uint64_t _sgn;
  };
  struct {
    uint64_t lo;
    uint64_t hi;
    int64_t ex;
    uint64_t sgn;
  };
} dint64_t;

typedef union {
  u128 r;
  struct {
    uint64_t l;
    uint64_t h;
  };
} uint128_t;

typedef union {
  double f;
  uint64_t u;
} f64_u;

// Extract both the mantissa and exponent of a double
static inline void fast_extract (int64_t *e, uint64_t *m, double x) {
  f64_u _x = {.f = x};

  *e = (_x.u >> 52) & 0x7ff;
  *m = (_x.u & (~0ul >> 12)) + (*e ? (1ul << 52) : 0);
  *e = *e - 0x3fe;
}

// Convert a non-zero double to the corresponding dint64_t value
static inline void dint_fromd (dint64_t *a, double b) {
  fast_extract (&a->ex, &a->hi, b);

  /* |b| = 2^(ex-52)*hi */

  uint32_t t = __builtin_clzl (a->hi);

  a->sgn = b < 0.0;
  a->hi = a->hi << t;
  a->ex = a->ex - (t > 11 ? t - 12 : 0);
  /* b = 2^ex*hi/2^64 where 1/2 <= hi/2^64 < 1 */
  a->lo = 0;
}

static inline void subnormalize_dint(dint64_t *a) {
  if (a->ex > -1023)
    return;

  uint64_t ex = -(1011 + a->ex);

  uint64_t hi = a->hi >> ex;
  uint64_t md = (a->hi >> (ex - 1)) & 0x1;
  uint64_t lo = (a->hi & (~0ul >> ex)) || a->lo;

  switch (fegetround()) {
  case FE_TONEAREST:
    hi += lo ? md : hi & md;
    break;
  case FE_DOWNWARD:
    hi += a->sgn & (md | lo);
    break;
  case FE_UPWARD:
    hi += (!a->sgn) & (md | lo);
    break;
  }

  a->hi = hi << ex;
  a->lo = 0;

  if (!a->hi) {
    a->ex++;
    a->hi = (1l << 63);
  }
}

// Convert a dint64_t value to a double
static inline double dint_tod(dint64_t *a) {
  subnormalize_dint (a);

  f64_u r = {.u = (a->hi >> 11) | (0x3ffl << 52)};

  double rd = 0.0;
  if ((a->hi >> 10) & 0x1)
    rd += 0x1p-53;

  if (a->hi & 0x3ff || a->lo)
    rd += 0x1p-54;

  if (a->sgn)
    rd = -rd;

  r.u = r.u | a->sgn << 63;
  r.f += rd;

  f64_u e;

  if (a->ex > -1023) { // The result is a normal double
    if (a->ex > 1023)
      if (a->ex == 1024) {
        r.f = r.f * 0x1p+1;
        e.f = 0x1p+1023;
      } else {
        r.f = 0x1.fffffffffffffp+1023;
        e.f = 0x1.fffffffffffffp+1023;
      }
    else
      e.u = ((a->ex + 1023) & 0x7ff) << 52;
  } else {
    if (a->ex < -1074) {
      if (a->ex == -1075) {
        r.f = r.f * 0x1p-1;
        e.f = 0x1p-1074;
      } else {
        r.f = 0x0.0000000000001p-1022;
        e.f = 0x0.0000000000001p-1022;
      }
    } else {
      e.u = 1l << (a->ex + 1074);
    }
  }

  return r.f * e.f;
}

/**************** end of code copied from dint.h and pow.[ch] ****************/

typedef union {double f; uint64_t u;} b64u64_u;

/* This table approximates 1/(2pi) downwards with precision 1216:
   1/(2*pi) ~ T[0]/2^64 + T[1]/2^128 + ... + T[i]/2^((i+1)*64) + ...
   Computed with computeT() from sin.sage. */
static uint64_t T[19] = {
   0x28be60db9391054a,
   0x7f09d5f47d4d3770,
   0x36d8a5664f10e410,
   0x7f9458eaf7aef158,
   0x6dc91b8e909374b8,
   0x1924bba82746487,
   0x3f877ac72c4a69cf,
   0xba208d7d4baed121,
   0x3a671c09ad17df90,
   0x4e64758e60d4ce7d,
   0x272117e2ef7e4a0e,
   0xc7fe25fff7816603,
   0xfbcbc462d6829b47,
   0xdb4d9fb3c9f2c26d,
   0xd3d18fd9a797fa8b,
   0x5d49eeb1faf97c5e, /* i=15 */
   0xcf41ce7de294a4ba,
   0x9afed7ec47e35742,
   0x1580cc11bf1edaea,
};

/* Approximate X/(2pi) mod 1. If Xin is the input value, and Xout the
   output value, we have:
   |Xout - (Xin/(2pi) mod 1)| < 2^-127
   Return non-zero if an actual reduction was done.
*/
static void
reduce (dint64_t *X, double x) // FIXME: remove argument x
{
  if (x == TRACE)
    printf ("hi=%lu ex=%ld\n", X->hi, X->ex);
  
  int e = X->ex;
  uint64_t c[5];
  u128 u;

  if (e <= 1) // simply multiply by T[0]/2^64 + T[1]/2^128
  {
    u = (u128) X->hi * (u128) T[1];
    c[1] = u >> 64;
    u = (u128) X->hi * (u128) T[0];
    X->lo = c[1] + u;
    X->hi = (u >> 64) + (X->lo < (uint64_t) u);
    /* the ignored part is at most 2^(ex-128) <= 2^-127 */
    if (x == TRACE)
      printf ("hi=%lu lo=%lu ex=%ld\n", X->hi, X->lo, X->ex);
    return;
  }

  // now 2 <= e <= 1024
  assert (2 <= e && e <= 1024);

  /* The upper 64-bit word X->hi corresponds to hi/2^64*2^e, if multiplied by
     T[i]/2^((i+1)*64) it yields hi*T[i]/2^128 * 2^(e-i*64).
     If e-64i <= -128, it contributes to less than 2^-128;
     if e-64i >= 128, it yields an integer, which is 0 modulo 1.
     We thus only consider the values of i such that -127 <= e-64i <= 127,
     i.e., (-127+e)/64 <= i <= (127+e)/64.
     Up to 4 consecutive values of T[i] can contribute (only 3 when e is a
     multiple of 64). */
  int i = (e < 127) ? 0 : (e - 127 + 64 - 1) / 64; // ceil((e-127)/64)
  if (x == TRACE) printf ("i=%d\n", i);
  // 0 <= i <= 15
  u = (u128) X->hi * (u128) T[i+3]; // i+3 <= 18
  /* we do not compute c[0] = u % 2^64, which does not contribute below */
  c[1] = u >> 64;
  u = (u128) X->hi * (u128) T[i+2];
  c[1] += u;
  c[2] = (u >> 64) + (c[1] < (uint64_t) u);
  u = (u128) X->hi * (u128) T[i+1];
  c[2] += u;
  c[3] = (u >> 64) + (c[2] < (uint64_t) u);
  if (x == TRACE) printf ("T[i]=%lx\n", T[i]);
  u = (u128) X->hi * (u128) T[i];
  if (x == TRACE) printf ("u=%lu:%lu\n", (uint64_t) (u>>64), (uint64_t) u);
  c[3] += u;
  c[4] = (u >> 64) + (c[3] < (uint64_t) u);

  int f = e - 64 * i; // hi*T[i]/2^128 is multiplied by 2^f
  /* {c, 5} = hi*(T[i]+T[i+1]/2^64+T[i+2]/2^128+T[i+3]/2^192) */
  assert (2 <= f && f <= 127);
  if (x == TRACE) printf ("f=%d\n", f);
  /* now shift c[0..4] by f bits to the left */
  if (f < 64)
  {
    X->hi = (c[4] << f) | (c[3] >> (64 - f));
    X->lo = (c[3] << f) | (c[2] >> (64 - f));
  }
  else if (f == 64)
  {
    X->hi = c[3];
    X->lo = c[2];
  }
  else /* 65 <= f <= 127 */
  {
    f -= 64; /* 1 <= f <= 63 */
    X->hi = (c[3] << f) | (c[2] >> (64 - f));
    X->lo = (c[2] << f) | (c[1] >> (64 - f));
  }
  /* the approximation error is at most 2 ulps = 2^-127:
     (a) the truncated part in the above shifts, which is less than 1 ulp,
         i.e., less than 2^-128
     (b) the ignored terms hi*T[i+4] + ..., which accumulate to less than
         1 ulp too.
  */
  X->ex = 0;
  if (x == TRACE)
    printf ("hi=%lu lo=%lu ex=%ld\n", X->hi, X->lo, X->ex);
}

/* return i and compute R such that X = i/2^8 + R */
static int
reduce2 (dint64_t *R, dint64_t *X)
{
  assert (X->ex <= 1);
  int i;
  if (X->ex <= -8)
  {
    i = 0;
    R->hi = X->hi;
  }
  else if (X->ex <= 0)
  {
    int sh = 64 - 8 - X->ex;
    i = X->hi >> sh;
    R->hi = X->hi & ((1ul << sh) - 1);
  }
  else /* X->ex = 1 */
  {
    i = X->hi >> 55;
    R->hi = X->hi & 0x7fffffffffffff;
  }
  R->lo = X->lo;
  R->ex = X->ex;
  R->sgn = X->sgn;
  return i;
}

static double
sin_accurate (double x)
{
  b64u64_u t = {.f = x};
  int e = (t.u >> 52) & 0x7ff;

  if (e == 0x7ff) /* NaN, +Inf and -Inf. */
    return 0.0 / 0.0;

  /* now x is a regular number */

  /* For |x| <= 0x1.7137449123ef6p-26, sin(x) rounds to x (to nearest):
     we can assume x >= 0 without loss of generality since sin(-x) = -sin(x),
     we have x - x^3/6 < sin(x) < x for say 0 < x <= 1 thus
     |sin(x) - x| < x^3/6.
     Write x = c*2^e with 1/2 <= c < 1.
     Then ulp(x)/2 = 2^(e-54), and x^3/6 = c^3/3*2^(3e), thus
     x^3/6 < ulp(x)/2 rewrites as c^3/6*2^(3e) < 2^(e-54),
     or c^3*2^(2e+53) < 3 (1).
     For e <= -26, since c^3 < 1, we have c^3*2^(2e+53) < 2 < 3.
     For e=-25, (1) rewrites 8*c^3 < 3 which yields c <= 0x1.7137449123ef6p-1.
  */
  uint64_t ux = t.u & 0x7fffffffffffffff;
  if (ux <= 0x3e57137449123ef6) // 0x3e57137449123ef6 = 0x1.7137449123ef6p-26
    return __builtin_fma (x, 0x1p-54, x);

  /* now |x| > 0x1.7137449123ef6p-26 */
  dint64_t X[1];
  dint_fromd (X, x);

  /* reduce argument */
  reduce (X, x);

  /* Now X = frac(x/(2pi)) + eps with |eps| < 2^-127, with |X| < 1.
     Write |X| = i/2^8 + r with r < 2^8. */
  dint64_t R[1];
  int i = reduce2 (R, X);

  return 0;
}

double
cr_sin (double x)
{
  return sin_accurate (x);
}
