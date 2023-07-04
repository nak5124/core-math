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
#include <stdlib.h>
#include <fenv.h>

#define TRACE -0x0p+0

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

// Return non-zero if a = 0
static inline int
dint_zero_p (const dint64_t *a)
{
  return a->hi == 0;
}

static inline char cmp(int64_t a, int64_t b) { return (a > b) - (a < b); }

static inline char cmpu128 (u128 a, u128 b) { return (a > b) - (a < b); }

/* ZERO is a dint64_t representation of 0, which ensures that
   dint_tod(ZERO) = 0 */
static const dint64_t ZERO = {.hi = 0x0, .lo = 0x0, .ex = -1076, .sgn = 0x0};
// MAGIC is a dint64_t representation of 1/2^11
static const dint64_t MAGIC = {.hi = 0x8000000000000000, .lo = 0x0, .ex = -10, .sgn = 0x0};

// Compare the absolute values of a and b
// Return -1 if |a| < |b|
// Return  0 if |a| = |b|
// Return +1 if |a| > |b|
static inline signed char
cmp_dint_abs (const dint64_t *a, const dint64_t *b) {
  if (dint_zero_p (a))
    return dint_zero_p (b) ? 0 : -1;
  if (dint_zero_p (b))
    return +1;
  char c1 = cmp (a->ex, b->ex);
  return c1 ? c1 : cmpu128 (a->r, b->r);
}

// Copy a dint64_t value
static inline void cp_dint(dint64_t *r, const dint64_t *a) {
  r->ex = a->ex;
  r->r = a->r;
  r->sgn = a->sgn;
}

// Add two dint64_t values, with error bounded by 2 ulps (ulp_128)
// (more precisely 1 ulp when a and b have same sign, 2 ulps otherwise)
// Moreover, when Sterbenz theorem applies, i.e., |b| <= |a| <= 2|b|
// and a,b are of different signs, there is no error, i.e., r = a-b.
static inline void
add_dint (dint64_t *r, const dint64_t *a, const dint64_t *b) {
  if (!(a->hi | a->lo)) {
    cp_dint (r, b);
    return;
  }

  switch (cmp_dint_abs (a, b)) {
  case 0:
    if (a->sgn ^ b->sgn) {
      cp_dint (r, &ZERO);
      return;
    }

    cp_dint (r, a);
    r->ex++;
    return;

  case -1: // |A| < |B|
    {
      // swap operands
      const dint64_t *tmp = a; a = b; b = tmp;
      break; // fall through the case |A| > |B|
    }
  }

  // From now on, |A| > |B| thus a->ex >= b->ex

  u128 A = a->r, B = b->r;
  uint64_t k = a->ex - b->ex;

  if (k > 0) {
    /* Warning: the right shift x >> k is only defined for 0 <= k < n
       where n is the bit-width of x. See for example
       https://developer.arm.com/documentation/den0024/a/The-A64-instruction-set/Data-processing-instructions/Shift-operations
       where is is said that k is interpreted modulo n. */
    B = (k < 128) ? B >> k : 0;
  }

  u128 C;
  unsigned char sgn = a->sgn;

  r->ex = a->ex; /* tentative exponent for the result */

  if (a->sgn ^ b->sgn) {
    /* a and b have different signs C = A + (-B)
       Sterbenz case |a|/2 <= |b| <= |a| can occur only when:
       * k=0: then B is not truncated, and C is exact below
       * k=1 and ex>0 below: then we ensure C is exact
     */
    C = A - B;
    uint64_t ch = C >> 64;
    /* We can't have C=0 here since we excluded the case |A| = |B|,
       thus __builtin_clzl(C) is well-defined below. */
    uint64_t ex = ch ? __builtin_clzl(ch) : 64 + __builtin_clzl(C);
    /* The error from the truncated part of B (1 ulp) is multiplied by 2^ex,
       thus by 2 ulps when ex <= 1. */
    if (ex > 0)
    {
      if (k == 1) /* Sterbenz case */
        C = (A << ex) - (b->r << (ex - 1));
      else
        C = (A << ex) - (B << ex);
      /* If C0 is the previous value of C, we have:
         (C0-1)*2^ex < A*2^ex-B*2^ex <= C0*2^ex
         since some neglected bits from B might appear which contribute
         a value less than ulp(C0)=1.
         As a consequence since 2^(127-ex) <= C0 < 2^(128-ex), because C0 had
         ex leading zero bits, we have 2^127-2^ex <= A*2^ex-B*2^ex < 2^128.
         Thus the value of C, which is truncated to 128 bits, is the right
         one (as if no truncation); moreover in some rare cases we need to
         shift by 1 bit to the left. */
      r->ex -= ex;
      ex = __builtin_clzl (C >> 64);
      /* Fall through with the code for ex = 0. */
    }
    C = C << ex;
    r->ex -= ex;
    /* The neglected part of B is bounded by 2 ulp(C) when ex=0, 1 ulp
       when ex > 0 but ex=0 at the end, and by 2*ulp(C) when ex > 0 and there
       is an extra shift at the end (in that case necessarily ex=1). */
  } else {
    C = A + B;
    if (C < A)
    {
      C = ((u128) 1 << 127) | (C >> 1);
      r->ex ++;
    }
  }

  /* In the addition case, we loose the truncated part of B, which
     contributes to at most 1 ulp. If there is an exponent shift, we
     might also loose the least significant bit of C, which counts as
     1/2 ulp, but the truncated part of B is now less than 1/2 ulp too,
     thus in all cases the error is less than 1 ulp(r). */

  r->sgn = sgn;
  r->r = C;
}

// Multiply two dint64_t numbers, with error bounded by 6 ulps
// on the 128-bit floating-point numbers.
// Overlap between r and a is allowed
static inline void
mul_dint (dint64_t *r, const dint64_t *a, const dint64_t *b) {
  u128 bh = b->hi, bl = b->lo;

  /* compute the two middle terms */
  u128 m1 = (u128)(a->hi) * bl;
  u128 m2 = (u128)(a->lo) * bh;

  /* put the 128-bit product of the high terms in r */
  r->r = (u128)(a->hi) * bh;

  /* there can be no overflow in the following addition since r <= (B-1)^2
     with B=2^64, (m1>>64) <= B-1 and (m2>>64) <= B-1, thus the sum is
     bounded by (B-1)^2+2*(B-1) = B^2-1 */
  r->r += (m1 >> 64) + (m2 >> 64);

  // Ensure that r->hi starts with a 1
  uint64_t ex = r->hi >> 63;
  r->r = r->r << (1 - ex);

  // Exponent and sign
  // if ex=1, then ex(r) = ex(a) + ex(b)
  // if ex=0, then ex(r) = ex(a) + ex(b) - 1
  r->ex = a->ex + b->ex + ex - 1;
  r->sgn = a->sgn ^ b->sgn;

  /* The ignored part can be as large as 3 ulps before the shift (one
     for the low part of a->hi * bl, one for the low part of a->lo * bh,
     and one for the neglected a->lo * bl term). After the shift this can
     be as large as 6 ulps. */
}

// Multiply two dint64_t numbers, assuming the low part of b is zero
// with error bounded by 2 ulps
static inline void
mul_dint_21 (dint64_t *r, const dint64_t *a, const dint64_t *b) {
  u128 bh = b->hi;
  u128 hi = (u128) (a->hi) * bh;
  u128 lo = (u128) (a->lo) * bh;

  /* put the 128-bit product of the high terms in r */
  r->r = hi;

  /* add the middle term */
  r->r += lo >> 64;

  // Ensure that r->hi starts with a 1
  uint64_t ex = r->hi >> 63;
  r->r = r->r << (1 - ex);

  // Exponent and sign
  r->ex = a->ex + b->ex + ex - 1;
  r->sgn = a->sgn ^ b->sgn;

  /* The ignored part can be as large as 1 ulp before the shift (truncated
     part of lo). After the shift this can be as large as 2 ulps. */
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

  if (a->ex > -1022) { // The result is a normal double
    if (a->ex > 1024)
      if (a->ex == 1025) {
        r.f = r.f * 0x1p+1;
        e.f = 0x1p+1023;
      } else {
        r.f = 0x1.fffffffffffffp+1023;
        e.f = 0x1.fffffffffffffp+1023;
      }
    else
      e.u = ((a->ex + 1022) & 0x7ff) << 52;
  } else {
    if (a->ex < -1073) {
      if (a->ex == -1074) {
        r.f = r.f * 0x1p-1;
        e.f = 0x1p-1074;
      } else {
        r.f = 0x0.0000000000001p-1022;
        e.f = 0x0.0000000000001p-1022;
      }
    } else {
      e.u = 1l << (a->ex + 1073);
    }
  }

  return r.f * e.f;
}

/**************** end of code copied from dint.h and pow.[ch] ****************/

typedef union {double f; uint64_t u;} b64u64_u;

#if 0
static void
print_dint (dint64_t *X)
{
  if (X->sgn == 0)
    printf ("2^%ld*(%lu/2^64+%lu/2^128)\n", X->ex, X->hi, X->lo);
  else
    printf ("-2^%ld*(%lu/2^64+%lu/2^128)\n", X->ex, X->hi, X->lo);
}
#endif

/* This table approximates 1/(2pi) downwards with precision 1280:
   1/(2*pi) ~ T[0]/2^64 + T[1]/2^128 + ... + T[i]/2^((i+1)*64) + ...
   Computed with computeT() from sin.sage. */
static const uint64_t T[20] = {
  0x28be60db9391054a, // i=0
   0x7f09d5f47d4d3770,
   0x36d8a5664f10e410,
   0x7f9458eaf7aef158,
   0x6dc91b8e909374b8,
   0x1924bba82746487, // i=5
   0x3f877ac72c4a69cf,
   0xba208d7d4baed121,
   0x3a671c09ad17df90,
   0x4e64758e60d4ce7d,
   0x272117e2ef7e4a0e, // i=10
   0xc7fe25fff7816603,
   0xfbcbc462d6829b47,
   0xdb4d9fb3c9f2c26d,
   0xd3d18fd9a797fa8b,
   0x5d49eeb1faf97c5e, // i=15
   0xcf41ce7de294a4ba,
   0x9afed7ec47e35742,
   0x1580cc11bf1edaea,
   0xfc33ef0826bd0d87, // i=19
};

/* Table containing 128-bit approximations of sin2pi(i/2^11) for 0 <= i < 256
   (to nearest).
   Each entry is to be interpreted as (hi/2^64+lo/2^128)*2^ex*(-1)*sgn.
   Generated with computeS() from sin.sage. */
static const dint64_t S[256] = {
  {.hi = 0x0, .lo = 0x0, .ex = 128, .sgn=0},
  {.hi = 0xc90fc5f66525d257, .lo = 0x480f7956b6470765, .ex = -8, .sgn=0},
  {.hi = 0xc90f87f3380388d5, .lo = 0xcb3ff35bd4d81baa, .ex = -7, .sgn=0},
  {.hi = 0x96cb587284b81770, .lo = 0xb767005691b9d9d1, .ex = -6, .sgn=0},
  {.hi = 0xc90e8fe6f63c2330, .lo = 0xf1d7d06db39ea9fc, .ex = -6, .sgn=0},
  {.hi = 0xfb514b55ccbe541a, .lo = 0xd784e031f9af76d6, .ex = -6, .sgn=0},
  {.hi = 0x96c9b5df1877e9b5, .lo = 0xf91ee371d6467dca, .ex = -5, .sgn=0},
  {.hi = 0xafea690fd5912ef3, .lo = 0xf56e3c87ae3c56df, .ex = -5, .sgn=0},
  {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=0},
  {.hi = 0xe22a7a6729d8e453, .lo = 0x850021e392744a4f, .ex = -5, .sgn=0},
  {.hi = 0xfb49b98e8e7807f6, .lo = 0xb21ccebc9caac3, .ex = -5, .sgn=0},
  {.hi = 0x8a342eda160bf5ae, .lo = 0xde5b1068d174be9c, .ex = -4, .sgn=0},
  {.hi = 0x96c32baca2ae68b4, .lo = 0x37b2dd49d5fca3c0, .ex = -4, .sgn=0},
  {.hi = 0xa351cb7fc30bc889, .lo = 0xb56007d16d4ad5a3, .ex = -4, .sgn=0},
  {.hi = 0xafe00694866a1b44, .lo = 0xcd34d2751c2e1da7, .ex = -4, .sgn=0},
  {.hi = 0xbc6dd52c3a342eb5, .lo = 0xf10bfca3d6464012, .ex = -4, .sgn=0},
  {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=0},
  {.hi = 0xd5880deafc18b534, .lo = 0x73d1472472f4a390, .ex = -4, .sgn=0},
  {.hi = 0xe214689606bf1676, .lo = 0x438b4a73aecd2541, .ex = -4, .sgn=0},
  {.hi = 0xeea037cc04764844, .lo = 0xc4e92d01a2f42935, .ex = -4, .sgn=0},
  {.hi = 0xfb2b73cfc106ff68, .lo = 0xf0a0e36a000c7350, .ex = -4, .sgn=0},
  {.hi = 0x83db0a7231831d8f, .lo = 0x60e782313f6161af, .ex = -3, .sgn=0},
  {.hi = 0x8a2009a6b84d9402, .lo = 0x77724a2b2a669bc4, .ex = -3, .sgn=0},
  {.hi = 0x9064b3a76a22640c, .lo = 0x56e0a8b0d177b55d, .ex = -3, .sgn=0},
  {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=0},
  {.hi = 0x9cecf8962d14c822, .lo = 0x50ffe4f5caa7f1fa, .ex = -3, .sgn=0},
  {.hi = 0xa3308bc93904ad69, .lo = 0xdec1b7f2768bdafa, .ex = -3, .sgn=0},
  {.hi = 0xa973ba526a6850d9, .lo = 0x76f8c63986598c79, .ex = -3, .sgn=0},
  {.hi = 0xafb68054d520c60b, .lo = 0xfdd2fc0936594c2d, .ex = -3, .sgn=0},
  {.hi = 0xb5f8d9f3cd8945d6, .lo = 0x924bef13600f9852, .ex = -3, .sgn=0},
  {.hi = 0xbc3ac352ead90abe, .lo = 0xeb13e106732687f1, .ex = -3, .sgn=0},
  {.hi = 0xc27c389609850433, .lo = 0xb228a03916371f6f, .ex = -3, .sgn=0},
  {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=0},
  {.hi = 0xcefdb7592542e1e9, .lo = 0x6b47b8c44e5b037e, .ex = -3, .sgn=0},
  {.hi = 0xd53db9224ae01bca, .lo = 0x7337412cf70716cb, .ex = -3, .sgn=0},
  {.hi = 0xdb7d3761c7b263b6, .lo = 0xbb286d23e11c8337, .ex = -3, .sgn=0},
  {.hi = 0xe1bc2e3cf616a7ac, .lo = 0x31883b30137c6e62, .ex = -3, .sgn=0},
  {.hi = 0xe7fa99d983ee098f, .lo = 0xeeb8f9c33340a2f2, .ex = -3, .sgn=0},
  {.hi = 0xee38765d74fe4897, .lo = 0xed16b994af6c18ae, .ex = -3, .sgn=0},
  {.hi = 0xf475bfef2551f5b9, .lo = 0x14e1a5488eaeab96, .ex = -3, .sgn=0},
  {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=0},
  {.hi = 0x8077456b7dc2d967, .lo = 0x3eac8308f1113e5e, .ex = -2, .sgn=0},
  {.hi = 0x8395023dd418e919, .lo = 0xdb1f70118c9c2198, .ex = -2, .sgn=0},
  {.hi = 0x86b26de5933c2e8e, .lo = 0xc5a9decdfaad4db5, .ex = -2, .sgn=0},
  {.hi = 0x89cf8676d7abb55b, .lo = 0x97965c9860c34e44, .ex = -2, .sgn=0},
  {.hi = 0x8cec4a05f12739e8, .lo = 0xdcdca90cc73b116a, .ex = -2, .sgn=0},
  {.hi = 0x9008b6a763de75b7, .lo = 0xa6e3df5975cca9da, .ex = -2, .sgn=0},
  {.hi = 0x9324ca6fe9a04b4e, .lo = 0x899c4de737feec22, .ex = -2, .sgn=0},
  {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=0},
  {.hi = 0x995bdfca28b53a54, .lo = 0x49c4863de522b217, .ex = -2, .sgn=0},
  {.hi = 0x9c76dd866c689dcc, .lo = 0xe7bc08111d0bfca4, .ex = -2, .sgn=0},
  {.hi = 0x9f917abeda4498df, .lo = 0xf3ff913a4aadb85e, .ex = -2, .sgn=0},
  {.hi = 0xa2abb58949f2ced7, .lo = 0xa5dbee6084ee1260, .ex = -2, .sgn=0},
  {.hi = 0xa5c58bfbcfd4436a, .lo = 0x69fcb11e19f58619, .ex = -2, .sgn=0},
  {.hi = 0xa8defc2cbe2f8fcc, .lo = 0xcd12a1f6ab6b095, .ex = -2, .sgn=0},
  {.hi = 0xabf80432a65ef190, .lo = 0x8c95c4c91179176b, .ex = -2, .sgn=0},
  {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=0},
  {.hi = 0xb228d418ec1869ad, .lo = 0x16031a34d4fc855d, .ex = -2, .sgn=0},
  {.hi = 0xb5409827b25591f0, .lo = 0xcd73fb5d8d45d302, .ex = -2, .sgn=0},
  {.hi = 0xb857ec684627fa4c, .lo = 0x187e26d290714d70, .ex = -2, .sgn=0},
  {.hi = 0xbb6ecef285f98a3a, .lo = 0xbddd8a0365d6b1d3, .ex = -2, .sgn=0},
  {.hi = 0xbe853dde9658dc60, .lo = 0xdfe1b074e22fc666, .ex = -2, .sgn=0},
  {.hi = 0xc19b3744e3262dcd, .lo = 0xad5a41de48f6b26f, .ex = -2, .sgn=0},
  {.hi = 0xc4b0b93e20c0213f, .lo = 0xdab4e426409b23a0, .ex = -2, .sgn=0},
  {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=0},
  {.hi = 0xcada4f4db157cf77, .lo = 0xfa6171200ab2efc3, .ex = -2, .sgn=0},
  {.hi = 0xcdee5f96e21b332c, .lo = 0x65a3132adfb7dfd5, .ex = -2, .sgn=0},
  {.hi = 0xd101f0d8c18ed1c1, .lo = 0xaadb580a1eba209f, .ex = -2, .sgn=0},
  {.hi = 0xd415012d802284f0, .lo = 0xdf4005ef6a64aa02, .ex = -2, .sgn=0},
  {.hi = 0xd7278eaf9dcd5b55, .lo = 0x1779df36d1cc8912, .ex = -2, .sgn=0},
  {.hi = 0xda399779eb391377, .lo = 0xcbabaeb97af8e8aa, .ex = -2, .sgn=0},
  {.hi = 0xdd4b19a78aed6515, .lo = 0xece7f445cecf1e28, .ex = -2, .sgn=0},
  {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=0},
  {.hi = 0xe36c829aeba6e720, .lo = 0x26a0eecdb4f16266, .ex = -2, .sgn=0},
  {.hi = 0xe67c659895943123, .lo = 0x82b0aecadf808123, .ex = -2, .sgn=0},
  {.hi = 0xe98bba6965ef725f, .lo = 0xb91caf23416e7e80, .ex = -2, .sgn=0},
  {.hi = 0xec9a7f2a2a188aeb, .lo = 0x7244ee20f591983b, .ex = -2, .sgn=0},
  {.hi = 0xefa8b1f8084ccdfc, .lo = 0x1050cdf22f34182f, .ex = -2, .sgn=0},
  {.hi = 0xf2b650f080d0da8d, .lo = 0x587f3fa044e2d27d, .ex = -2, .sgn=0},
  {.hi = 0xf5c35a316f1a3c80, .lo = 0x643720de93ba81bd, .ex = -2, .sgn=0},
  {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=0},
  {.hi = 0xfbdba405e9c00cca, .lo = 0xd24d3023da491920, .ex = -2, .sgn=0},
  {.hi = 0xfee6e0d6ff6fc5a4, .lo = 0x8b74fe2508ab8fc2, .ex = -2, .sgn=0},
  {.hi = 0x80f8c035cfee8d76, .lo = 0xfd958d68e8b49e6b, .ex = -1, .sgn=0},
  {.hi = 0x827dc071bfed6ffa, .lo = 0xfb4c92369f0cf008, .ex = -1, .sgn=0},
  {.hi = 0x8402702f5b30f2a9, .lo = 0xcb07b25a7b0372a7, .ex = -1, .sgn=0},
  {.hi = 0x8586ce7ededc809d, .lo = 0x9d3dc689006896f4, .ex = -1, .sgn=0},
  {.hi = 0x870ada70ba4e6d49, .lo = 0x9d52755ece3f70, .ex = -1, .sgn=0},
  {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=0},
  {.hi = 0x8a11f77e349bc245, .lo = 0xa66d1d936c38c329, .ex = -1, .sgn=0},
  {.hi = 0x8b9506bbb28bb922, .lo = 0x575f33366be0afef, .ex = -1, .sgn=0},
  {.hi = 0x8d17bfdf47921ac8, .lo = 0xcb590d74f64e77c9, .ex = -1, .sgn=0},
  {.hi = 0x8e9a21fa66d9ee8d, .lo = 0xf2be3ecae62789d4, .ex = -1, .sgn=0},
  {.hi = 0x901c2c1eb93dee39, .lo = 0x632b9cff5cfee724, .ex = -1, .sgn=0},
  {.hi = 0x919ddd5e1ddb8b33, .lo = 0x609c464b3dd676ec, .ex = -1, .sgn=0},
  {.hi = 0x931f34caaaa5d23a, .lo = 0x6a1ff8bfe6396e28, .ex = -1, .sgn=0},
  {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=0},
  {.hi = 0x9620d274aa290339, .lo = 0xe06a955a5b8e301d, .ex = -1, .sgn=0},
  {.hi = 0x97a116d7601c3515, .lo = 0xfc8b7184b21f2d50, .ex = -1, .sgn=0},
  {.hi = 0x9920fdb1c5d5783d, .lo = 0x9dd1eedf18a2e4df, .ex = -1, .sgn=0},
  {.hi = 0x9aa086170c0a8d86, .lo = 0x9ffa0d23f3c26c62, .ex = -1, .sgn=0},
  {.hi = 0x9c1faf1a9db554af, .lo = 0xdab6b478577e7be5, .ex = -1, .sgn=0},
  {.hi = 0x9d9e77d020a5bbe6, .lo = 0xdb895384528d0d60, .ex = -1, .sgn=0},
  {.hi = 0x9f1cdf4b76138b02, .lo = 0x98dbd3555ebcdefe, .ex = -1, .sgn=0},
  {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=0},
  {.hi = 0xa21886e449b78316, .lo = 0xd29d23a624acd00c, .ex = -1, .sgn=0},
  {.hi = 0xa395c52ab8829dfc, .lo = 0x2be036401ba87cc2, .ex = -1, .sgn=0},
  {.hi = 0xa5129e88dc17976a, .lo = 0x82d9495ead5be348, .ex = -1, .sgn=0},
  {.hi = 0xa68f1213c73b5124, .lo = 0x17218792857f4c5a, .ex = -1, .sgn=0},
  {.hi = 0xa80b1ee0cb823c27, .lo = 0x3269f4702b88324a, .ex = -1, .sgn=0},
  {.hi = 0xa986c40579e11c0a, .lo = 0x8e3bdf8085321556, .ex = -1, .sgn=0},
  {.hi = 0xab020097a33da341, .lo = 0xc1654b64a0081b46, .ex = -1, .sgn=0},
  {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=0},
  {.hi = 0xadf73c5ced9db0f3, .lo = 0x9a5318ac6fe94e4d, .ex = -1, .sgn=0},
  {.hi = 0xaf7139bcf5349ac6, .lo = 0x9fe5f4ea48965e2c, .ex = -1, .sgn=0},
  {.hi = 0xb0eacae4461013ed, .lo = 0x63c66682bae74898, .ex = -1, .sgn=0},
  {.hi = 0xb263eee9f93e3088, .lo = 0x695a5332090bb09b, .ex = -1, .sgn=0},
  {.hi = 0xb3dca4e56b1e54bb, .lo = 0x992d96e5021e3c37, .ex = -1, .sgn=0},
  {.hi = 0xb554ebee3bf0b58e, .lo = 0x971f4da709ad4378, .ex = -1, .sgn=0},
  {.hi = 0xb6ccc31c5065afee, .lo = 0x35ebacd79f209137, .ex = -1, .sgn=0},
  {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=0},
  {.hi = 0xb9bb1e4930848ead, .lo = 0xcdb0531c4e58484b, .ex = -1, .sgn=0},
  {.hi = 0xbb31a07920c7b256, .lo = 0x55b92083658bb897, .ex = -1, .sgn=0},
  {.hi = 0xbca7af309efd7182, .lo = 0xa4b0d21fc5036a5, .ex = -1, .sgn=0},
  {.hi = 0xbe1d4988ee67380c, .lo = 0xd1f90f79f46c7e01, .ex = -1, .sgn=0},
  {.hi = 0xbf926e9b9a0f2127, .lo = 0x91a1b5eb79658c67, .ex = -1, .sgn=0},
  {.hi = 0xc1071d8275561f9b, .lo = 0x721853f8e528a934, .ex = -1, .sgn=0},
  {.hi = 0xc27b55579c81f96d, .lo = 0xcdc2bd470675104d, .ex = -1, .sgn=0},
  {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=0},
  {.hi = 0xc5625c36af6a222f, .lo = 0xf4ff2895ab6ebe89, .ex = -1, .sgn=0},
  {.hi = 0xc6d5297645257e8d, .lo = 0x14d24739de27e2e9, .ex = -1, .sgn=0},
  {.hi = 0xc8477c0f7bde8a98, .lo = 0x4ce0246ad4fa74, .ex = -1, .sgn=0},
  {.hi = 0xc9b9531de49eb968, .lo = 0x4319e5ad5b0dcb84, .ex = -1, .sgn=0},
  {.hi = 0xcb2aadbd5ca47af5, .lo = 0xfaa3dfe675a65ee2, .ex = -1, .sgn=0},
  {.hi = 0xcc9b8b0a0deff5d4, .lo = 0x2e663b3c7555a6c3, .ex = -1, .sgn=0},
  {.hi = 0xce0bea206fcf9192, .lo = 0x3c540a9eec47af38, .ex = -1, .sgn=0},
  {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=0},
  {.hi = 0xd0eb2a1da855fefd, .lo = 0xb9302788604e88f1, .ex = -1, .sgn=0},
  {.hi = 0xd25a093ef50f2482, .lo = 0x721fc87ba1d42456, .ex = -1, .sgn=0},
  {.hi = 0xd3c8669edf98d680, .lo = 0x87967926fdcecec4, .ex = -1, .sgn=0},
  {.hi = 0xd536415b69fe4c54, .lo = 0x1df22346611c6b4b, .ex = -1, .sgn=0},
  {.hi = 0xd6a39892e6e04764, .lo = 0x3090d44db12c418c, .ex = -1, .sgn=0},
  {.hi = 0xd8106b63fa0048a0, .lo = 0xa573f2aa90434ba5, .ex = -1, .sgn=0},
  {.hi = 0xd97cb8ed98cb93f5, .lo = 0x2e349483e3fb2a6a, .ex = -1, .sgn=0},
  {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=0},
  {.hi = 0xdc53c0a7eab49b35, .lo = 0x3ccca3982328ed8b, .ex = -1, .sgn=0},
  {.hi = 0xddbe791825e8099e, .lo = 0x1a5bd9269d408d7e, .ex = -1, .sgn=0},
  {.hi = 0xdf28a8bffe06ca56, .lo = 0xcce2634be2bf54df, .ex = -1, .sgn=0},
  {.hi = 0xe0924ec008f734fd, .lo = 0x8aa895d5bf3e84ea, .ex = -1, .sgn=0},
  {.hi = 0xe1fb6a3931894b38, .lo = 0xf7a1f9bd9ba13b6b, .ex = -1, .sgn=0},
  {.hi = 0xe363fa4cb8005482, .lo = 0x7b32c72e31824e51, .ex = -1, .sgn=0},
  {.hi = 0xe4cbfe1c329c453a, .lo = 0xd40e9e6b989f89e5, .ex = -1, .sgn=0},
  {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=0},
  {.hi = 0xe79a5d770e6905dc, .lo = 0xf1b65cc5fd780262, .ex = -1, .sgn=0},
  {.hi = 0xe900b7474edad637, .lo = 0x431626c10485bdda, .ex = -1, .sgn=0},
  {.hi = 0xea66815d4304e6c8, .lo = 0xcc39cfcc29960b1, .ex = -1, .sgn=0},
  {.hi = 0xebcbbadc371c4aaa, .lo = 0x1d90f780ae951140, .ex = -1, .sgn=0},
  {.hi = 0xed3062e7d086c6f0, .lo = 0xc71debc372b6f9d4, .ex = -1, .sgn=0},
  {.hi = 0xee9478a40e62bf86, .lo = 0x2a24164daec85ccb, .ex = -1, .sgn=0},
  {.hi = 0xeff7fb354a0eecb1, .lo = 0x527233b40d3432bb, .ex = -1, .sgn=0},
  {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=0},
  {.hi = 0xf2bd4369e6c126d3, .lo = 0x7f232aee178c6323, .ex = -1, .sgn=0},
  {.hi = 0xf41f0757c2889e84, .lo = 0x3c7f10db458c337c, .ex = -1, .sgn=0},
  {.hi = 0xf58034af92b102a7, .lo = 0x93fa6107c4327527, .ex = -1, .sgn=0},
  {.hi = 0xf6e0ca977bc6ac45, .lo = 0xe1079824233fef46, .ex = -1, .sgn=0},
  {.hi = 0xf840c835ffbfed66, .lo = 0xa9a56012067c570c, .ex = -1, .sgn=0},
  {.hi = 0xf9a02cb1fe833a0d, .lo = 0x8da894471de1a18, .ex = -1, .sgn=0},
  {.hi = 0xfafef732b66d1742, .lo = 0x343fbf4a7d42af3, .ex = -1, .sgn=0},
  {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=0},
  {.hi = 0xfdbabae12696eea4, .lo = 0x2377c3799c052fa, .ex = -1, .sgn=0},
  {.hi = 0xff17b25f38907dad, .lo = 0xa9c6ba50490539f, .ex = -1, .sgn=0},
  {.hi = 0x803a06415c170525, .lo = 0x6f53873e2f1477ff, .ex = 0, .sgn=0},
  {.hi = 0x80e7e43a61f5b6cb, .lo = 0x5ca183dc973abc22, .ex = 0, .sgn=0},
  {.hi = 0x819572af6decac84, .lo = 0x9fba97fdf0c4d24c, .ex = 0, .sgn=0},
  {.hi = 0x8242b1357110d372, .lo = 0x6fb2123fedfa6e22, .ex = 0, .sgn=0},
  {.hi = 0x82ef9f618dc5b70e, .lo = 0x91a965931f1a200a, .ex = 0, .sgn=0},
  {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=0},
  {.hi = 0x8448890195846099, .lo = 0x246efcff30cb064a, .ex = 0, .sgn=0},
  {.hi = 0x84f483a0be2f0403, .lo = 0x51917cac857fd5f5, .ex = 0, .sgn=0},
  {.hi = 0x85a02c3c7c2f5ca5, .lo = 0x327888fe4b62687b, .ex = 0, .sgn=0},
  {.hi = 0x864b826aec4c74e5, .lo = 0x85043222c9bdd18d, .ex = 0, .sgn=0},
  {.hi = 0x86f685c25e25acf5, .lo = 0x7e0b9b07548471a2, .ex = 0, .sgn=0},
  {.hi = 0x87a135d95473ec89, .lo = 0x4e091160e2430712, .ex = 0, .sgn=0},
  {.hi = 0x884b9246854ab50b, .lo = 0x4f14c8afe4560291, .ex = 0, .sgn=0},
  {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=0},
  {.hi = 0x899f4e7f712a765e, .lo = 0xc88302a31afce54a, .ex = 0, .sgn=0},
  {.hi = 0x8a48ad799b6759f3, .lo = 0x660558a02136130a, .ex = 0, .sgn=0},
  {.hi = 0x8af1b726df15e13c, .lo = 0x545f7d79ead8fa19, .ex = 0, .sgn=0},
  {.hi = 0x8b9a6b1ef6da4502, .lo = 0x21a6675f51580bc4, .ex = 0, .sgn=0},
  {.hi = 0x8c42c8f9d2372644, .lo = 0x101a5adbcb9ffb43, .ex = 0, .sgn=0},
  {.hi = 0x8cead04f95cdbf66, .lo = 0x4d49cbaf15aecd80, .ex = 0, .sgn=0},
  {.hi = 0x8d9280b89b9df49b, .lo = 0xde2d43c6b67a7cbe, .ex = 0, .sgn=0},
  {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=0},
  {.hi = 0x8ee0db26e24390f8, .lo = 0xaf0e2345f3bd24b4, .ex = 0, .sgn=0},
  {.hi = 0x8f87845de430d777, .lo = 0x9311a82459aa0f72, .ex = 0, .sgn=0},
  {.hi = 0x902dd50bab06b1b7, .lo = 0xb144016c7a30b39a, .ex = 0, .sgn=0},
  {.hi = 0x90d3ccc99f5ac58b, .lo = 0x9d1072e09b72292, .ex = 0, .sgn=0},
  {.hi = 0x91796b31609f0c54, .lo = 0x6714fe6925b78cc4, .ex = 0, .sgn=0},
  {.hi = 0x921eafdcc560f9c5, .lo = 0x33d0a284a8c954ad, .ex = 0, .sgn=0},
  {.hi = 0x92c39a65db88809d, .lo = 0x1f8481e704e4a767, .ex = 0, .sgn=0},
  {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=0},
  {.hi = 0x940c5f7a69e5ce1c, .lo = 0x1489a97671a42, .ex = 0, .sgn=0},
  {.hi = 0x94b0393b14e54156, .lo = 0xd6c7af02d5c16fd9, .ex = 0, .sgn=0},
  {.hi = 0x9553b743d75ac03f, .lo = 0xac0106650f4ef023, .ex = 0, .sgn=0},
  {.hi = 0x95f6d92fd79f4fba, .lo = 0xd9f8e1a446e973b9, .ex = 0, .sgn=0},
  {.hi = 0x96999e9a74ddbde3, .lo = 0xa7a7556c3b33abc1, .ex = 0, .sgn=0},
  {.hi = 0x973c071f4750b49c, .lo = 0xc0a03934f0cce19b, .ex = 0, .sgn=0},
  {.hi = 0x97de125a2080a8ed, .lo = 0xd243aa0843a2c144, .ex = 0, .sgn=0},
  {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=0},
  {.hi = 0x99210f624d30facb, .lo = 0xc4b992a37fb9b9bd, .ex = 0, .sgn=0},
  {.hi = 0x99c200686472b4a8, .lo = 0x1ab42d43235757b6, .ex = 0, .sgn=0},
  {.hi = 0x9a6292960a6f0ab0, .lo = 0x7e92c655656e6b85, .ex = 0, .sgn=0},
  {.hi = 0x9b02c58832cf95c0, .lo = 0x698b94f50326a043, .ex = 0, .sgn=0},
  {.hi = 0x9ba298dc0bfc6a88, .lo = 0x9a5614e8ffbeac6f, .ex = 0, .sgn=0},
  {.hi = 0x9c420c2eff590e5f, .lo = 0xc7fd954194e6d8aa, .ex = 0, .sgn=0},
  {.hi = 0x9ce11f1eb18147b1, .lo = 0x3e93627de8fd5779, .ex = 0, .sgn=0},
  {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=0},
  {.hi = 0x9e1e224c0e28bc94, .lo = 0x2cad377d5c9c35d8, .ex = 0, .sgn=0},
  {.hi = 0x9ebc11c62c1a1dfb, .lo = 0xcc141e10c6460c8b, .ex = 0, .sgn=0},
  {.hi = 0x9f599f55f0340061, .lo = 0xa88d5f46834bbf8d, .ex = 0, .sgn=0},
  {.hi = 0x9ff6ca9a2ab6a26d, .lo = 0x22cc118a0c118aa0, .ex = 0, .sgn=0},
  {.hi = 0xa0939331e8846237, .lo = 0x7cec6df5bea167cf, .ex = 0, .sgn=0},
  {.hi = 0xa12ff8bc735d8af6, .lo = 0x71acea2819360c35, .ex = 0, .sgn=0},
  {.hi = 0xa1cbfad9521bfd1b, .lo = 0x166c36e7bb3c402f, .ex = 0, .sgn=0},
  {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=0},
  {.hi = 0xa302d34959951243, .lo = 0x9443372e20d4377c, .ex = 0, .sgn=0},
  {.hi = 0xa39da8dcc39a38e5, .lo = 0xca9a8a720d4c69c, .ex = 0, .sgn=0},
  {.hi = 0xa4381983048ff747, .lo = 0xbf623cf5301a2dde, .ex = 0, .sgn=0},
  {.hi = 0xa4d224dcd849c5b0, .lo = 0x23d251cc8d7975cc, .ex = 0, .sgn=0},
  {.hi = 0xa56bca8b391785db, .lo = 0x189d39ffe11aaa2b, .ex = 0, .sgn=0},
  {.hi = 0xa6050a2f60002049, .lo = 0x8c33ebf3aa8501fb, .ex = 0, .sgn=0},
  {.hi = 0xa69de36ac4fbfadc, .lo = 0x9b3ad6e4022183d9, .ex = 0, .sgn=0},
  {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=0},
  {.hi = 0xa7ce612e65243291, .lo = 0x6b2a39f856a69781, .ex = 0, .sgn=0},
  {.hi = 0xa86604facd04d969, .lo = 0x3463a2c2e6e9cc55, .ex = 0, .sgn=0},
  {.hi = 0xa8fd40e6ccd52ffd, .lo = 0x6cc14c4f53e2e82d, .ex = 0, .sgn=0},
  {.hi = 0xa99414951aacae5e, .lo = 0xd147625fda929af8, .ex = 0, .sgn=0},
  {.hi = 0xaa2a7fa8acefdd63, .lo = 0xb714ee81b53b4b9d, .ex = 0, .sgn=0},
  {.hi = 0xaac081c4ba89ba8a, .lo = 0xe1b3dfc4dbda9bfd, .ex = 0, .sgn=0},
  {.hi = 0xab561a8cbb24f410, .lo = 0xf17cee69b0d2ecde, .ex = 0, .sgn=0},
  {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=0},
  {.hi = 0xac800eafb91ef9a9, .lo = 0xf86ba0dde982fb59, .ex = 0, .sgn=0},
  {.hi = 0xad146952eb9282af, .lo = 0x44bf16268608db96, .ex = 0, .sgn=0},
  {.hi = 0xada859327ba24151, .lo = 0x9d30d4cfeb04f1fb, .ex = 0, .sgn=0},
  {.hi = 0xae3bddf3280c620d, .lo = 0x3d53817865422565, .ex = 0, .sgn=0},
  {.hi = 0xaecef739f1a2df10, .lo = 0xf74d099042e8f326, .ex = 0, .sgn=0},
  {.hi = 0xaf61a4ac1b83a1de, .lo = 0xa89a9b8f726b95bf, .ex = 0, .sgn=0},
  {.hi = 0xaff3e5ef2b507c06, .lo = 0x8c679e67fc462d51, .ex = 0, .sgn=0},
  {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=0},
  {.hi = 0xb117227f6117f9f9, .lo = 0x8d8be132d576e614, .ex = 0, .sgn=0},
  {.hi = 0xb1a81d18e0df4889, .lo = 0x24784f32c3e3e5bd, .ex = 0, .sgn=0},
  {.hi = 0xb238aa1bfa9ad507, .lo = 0x8cc7d4bd05ffd5ae, .ex = 0, .sgn=0},
  {.hi = 0xb2c8c92f83c1eb87, .lo = 0xac9f7ebbc469ef59, .ex = 0, .sgn=0},
  {.hi = 0xb35879fa959c323c, .lo = 0x5d6635109164f740, .ex = 0, .sgn=0},
  {.hi = 0xb3e7bc248d78802e, .lo = 0xa156468ef6c18c60, .ex = 0, .sgn=0},
  {.hi = 0xb4768f550ce389fd, .lo = 0x4a85350f69018c55, .ex = 0, .sgn=0},
};

/* Table containing 128-bit approximations of cos2pi(i/2^11) for 0 <= i < 256
   (to nearest).
   Each entry is to be interpreted as (hi/2^64+lo/2^128)*2^ex*(-1)*sgn.
   Generated with computeC() from sin.sage. */
static const dint64_t C[256] = {
  {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=0},
  {.hi = 0xffffb10b10e80e95, .lo = 0x3031437d7eccb9df, .ex = 0, .sgn=0},
  {.hi = 0xfffec42c7454926b, .lo = 0x38e310779edfec68, .ex = 0, .sgn=0},
  {.hi = 0xfffd3964bc6275ba, .lo = 0x69fff9ae0dedb047, .ex = 0, .sgn=0},
  {.hi = 0xfffb10b4dc96dabb, .lo = 0xb47903f7a19f8ee2, .ex = 0, .sgn=0},
  {.hi = 0xfff84a1e29de8571, .lo = 0x8cc193c5d508e13f, .ex = 0, .sgn=0},
  {.hi = 0xfff4e5a25a8d095b, .lo = 0x43366df666fd54ff, .ex = 0, .sgn=0},
  {.hi = 0xfff0e343865bbb13, .lo = 0x5428ed0647c9e5d1, .ex = 0, .sgn=0},
  {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=0},
  {.hi = 0xffe704e71533c508, .lo = 0x53aa9423bb0adc21, .ex = 0, .sgn=0},
  {.hi = 0xffe128ef8e9fc17a, .lo = 0x7d209f32d42d864e, .ex = 0, .sgn=0},
  {.hi = 0xffdaaf212fed72db, .lo = 0x4fd8f038449ec436, .ex = 0, .sgn=0},
  {.hi = 0xffd3977ff7bae4e9, .lo = 0x664649b4d541b9c5, .ex = 0, .sgn=0},
  {.hi = 0xffcbe2104600a0a9, .lo = 0x5595ca3f421ae09c, .ex = 0, .sgn=0},
  {.hi = 0xffc38ed6dc0ef98b, .lo = 0x1c676208aa3be545, .ex = 0, .sgn=0},
  {.hi = 0xffba9dd8dc8b1e83, .lo = 0xccfed60a91097c48, .ex = 0, .sgn=0},
  {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=0},
  {.hi = 0xffa6e2a58df6947d, .lo = 0xd2c665c2da3e7844, .ex = 0, .sgn=0},
  {.hi = 0xff9c187c6abade6a, .lo = 0x1e1862cca089938b, .ex = 0, .sgn=0},
  {.hi = 0xff90b0a7098f6443, .lo = 0x2dabd3195a05710f, .ex = 0, .sgn=0},
  {.hi = 0xff84ab2c738d6a03, .lo = 0x519c314973ccae6b, .ex = 0, .sgn=0},
  {.hi = 0xff780814130c893c, .lo = 0x3ea4f30adda3016f, .ex = 0, .sgn=0},
  {.hi = 0xff6ac765b39e1e19, .lo = 0x1b9d5851979f28fb, .ex = 0, .sgn=0},
  {.hi = 0xff5ce92982087867, .lo = 0x50a7bb6a6ee3b0f1, .ex = 0, .sgn=0},
  {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=0},
  {.hi = 0xff3f542a416b0134, .lo = 0xb085c1828f69296a, .ex = 0, .sgn=0},
  {.hi = 0xff2f9d7971ca0364, .lo = 0x27e31939e2eec09c, .ex = 0, .sgn=0},
  {.hi = 0xff1f495f4ec430d7, .lo = 0xf5971326a3540ea9, .ex = 0, .sgn=0},
  {.hi = 0xff0e57e5ead848d1, .lo = 0x1f1901544271c3f8, .ex = 0, .sgn=0},
  {.hi = 0xfefcc917b99839a5, .lo = 0xe0abd3a9b64df725, .ex = 0, .sgn=0},
  {.hi = 0xfeea9cff8fa2ae54, .lo = 0xec34413e87ef2740, .ex = 0, .sgn=0},
  {.hi = 0xfed7d3a8a29c603b, .lo = 0x2f88b949a72ff96c, .ex = 0, .sgn=0},
  {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=0},
  {.hi = 0xfeb0696d3ae4f04d, .lo = 0xb7b6cc53c3abc817, .ex = 0, .sgn=0},
  {.hi = 0xfe9bc8a1105c22a5, .lo = 0xd3af6ee4f2101c20, .ex = 0, .sgn=0},
  {.hi = 0xfe868ac6c3043b2e, .lo = 0xb4f70c910505e10, .ex = 0, .sgn=0},
  {.hi = 0xfe70afeb6d33d6a2, .lo = 0x2907cf2b3f6feac2, .ex = 0, .sgn=0},
  {.hi = 0xfe5a381c8a1aa224, .lo = 0xd54faa364b7da8f6, .ex = 0, .sgn=0},
  {.hi = 0xfe432367f5b90a62, .lo = 0x87b8875373a818a4, .ex = 0, .sgn=0},
  {.hi = 0xfe2b71dbecd7aefc, .lo = 0x8598c2c429caf7, .ex = 0, .sgn=0},
  {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=0},
  {.hi = 0xfdfa3878546c3d28, .lo = 0x9bfe5c51e91cbdcd, .ex = 0, .sgn=0},
  {.hi = 0xfde0b0bf220c2fd4, .lo = 0xe276d247626a23fd, .ex = 0, .sgn=0},
  {.hi = 0xfdc68c6b356db62f, .lo = 0x499ddb331d19539d, .ex = 0, .sgn=0},
  {.hi = 0xfdabcb8caeba091b, .lo = 0xfac7397cc07a6470, .ex = 0, .sgn=0},
  {.hi = 0xfd906e340eaa6401, .lo = 0xd6e270740a186977, .ex = 0, .sgn=0},
  {.hi = 0xfd747472367dd6c5, .lo = 0x61beb8cd2696fc78, .ex = 0, .sgn=0},
  {.hi = 0xfd57de5867eedc39, .lo = 0x6c696582f346fd91, .ex = 0, .sgn=0},
  {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=0},
  {.hi = 0xfd1cdd63d0bc8735, .lo = 0x863b87258f11ad7e, .ex = 0, .sgn=0},
  {.hi = 0xfcfe72ad6d9641f2, .lo = 0xa06fab9f9d106709, .ex = 0, .sgn=0},
  {.hi = 0xfcdf6be7def1464c, .lo = 0xa4e064308f4999f4, .ex = 0, .sgn=0},
  {.hi = 0xfcbfc926484cd43a, .lo = 0xa3e22b4d38917e73, .ex = 0, .sgn=0},
  {.hi = 0xfc9f8a7c2d603c60, .lo = 0x5d582cac7cb4391c, .ex = 0, .sgn=0},
  {.hi = 0xfc7eaffd720ed673, .lo = 0x2880268f2e62955, .ex = 0, .sgn=0},
  {.hi = 0xfc5d39be5a5bbc4b, .lo = 0x1c0d254b6c8da4bd, .ex = 0, .sgn=0},
  {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=0},
  {.hi = 0xfc187a52063060c2, .lo = 0x9433b49289417ea2, .ex = 0, .sgn=0},
  {.hi = 0xfbf5314f31eb7375, .lo = 0x25aafd7fdba12c5f, .ex = 0, .sgn=0},
  {.hi = 0xfbd14ce0d191516e, .lo = 0x7190c94899dff1b8, .ex = 0, .sgn=0},
  {.hi = 0xfbaccd1d0903bb09, .lo = 0xe63ae8632b84473c, .ex = 0, .sgn=0},
  {.hi = 0xfb87b21a5bf5b917, .lo = 0x75df66f0ec3dd459, .ex = 0, .sgn=0},
  {.hi = 0xfb61fbefadddb985, .lo = 0x61ce9d5ef5a81487, .ex = 0, .sgn=0},
  {.hi = 0xfb3baab441e770f7, .lo = 0xb4b54683879c9c17, .ex = 0, .sgn=0},
  {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=0},
  {.hi = 0xfaed376a1b42e559, .lo = 0x2079880c450348ac, .ex = 0, .sgn=0},
  {.hi = 0xfac5158bc4f4211f, .lo = 0x4a188aa367f90ab1, .ex = 0, .sgn=0},
  {.hi = 0xfa9c58fd796837d4, .lo = 0x10655ecd5cc771d8, .ex = 0, .sgn=0},
  {.hi = 0xfa7301d859796671, .lo = 0x1fe196a53fb5b237, .ex = 0, .sgn=0},
  {.hi = 0xfa491035e55da3a3, .lo = 0xd24377c77a591e24, .ex = 0, .sgn=0},
  {.hi = 0xfa1e842ffc96e4e0, .lo = 0x431c393c7f62da65, .ex = 0, .sgn=0},
  {.hi = 0xf9f35de0dde328ab, .lo = 0xba5dbf4510eddc8f, .ex = 0, .sgn=0},
  {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=0},
  {.hi = 0xf99b42d1d57781eb, .lo = 0x78685d850f80ecdc, .ex = 0, .sgn=0},
  {.hi = 0xf96e4e4844d4e82a, .lo = 0x80e8c17bf80e8f02, .ex = 0, .sgn=0},
  {.hi = 0xf940bfe2304e6c45, .lo = 0xc0e2a1352ed7f292, .ex = 0, .sgn=0},
  {.hi = 0xf91297bbb1d6cdbe, .lo = 0x68fc6e4d6a920bd2, .ex = 0, .sgn=0},
  {.hi = 0xf8e3d5f1423842a0, .lo = 0x9701914c7f8fbcd7, .ex = 0, .sgn=0},
  {.hi = 0xf8b47a9fb902e76c, .lo = 0xac9f07f54ff5bc14, .ex = 0, .sgn=0},
  {.hi = 0xf88485e44c7af48a, .lo = 0xb36a9dfaadafc1e1, .ex = 0, .sgn=0},
  {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=0},
  {.hi = 0xf822d0a67b9c5cb5, .lo = 0xa776175bd284fe05, .ex = 0, .sgn=0},
  {.hi = 0xf7f110605caf6390, .lo = 0xa76f7efc19aed41c, .ex = 0, .sgn=0},
  {.hi = 0xf7beb728e51dfcb8, .lo = 0x730785813f78aa1e, .ex = 0, .sgn=0},
  {.hi = 0xf78bc51f239e12c6, .lo = 0x214cffcee9dd33ca, .ex = 0, .sgn=0},
  {.hi = 0xf7583a62852a23b2, .lo = 0x4becad887680c197, .ex = 0, .sgn=0},
  {.hi = 0xf7241712d4edde49, .lo = 0xf99107e50d631330, .ex = 0, .sgn=0},
  {.hi = 0xf6ef5b503c328589, .lo = 0x50ca117eb18beed7, .ex = 0, .sgn=0},
  {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=0},
  {.hi = 0xf6841af4cc8048a4, .lo = 0xce8c455197cdf8a7, .ex = 0, .sgn=0},
  {.hi = 0xf64d969e1dfc2119, .lo = 0x119d358de0493956, .ex = 0, .sgn=0},
  {.hi = 0xf6167a58d7b59026, .lo = 0x9dc7e5954c5a8f24, .ex = 0, .sgn=0},
  {.hi = 0xf5dec646f85ba1c6, .lo = 0xc8c615e72768d6b5, .ex = 0, .sgn=0},
  {.hi = 0xf5a67a8adc4088ca, .lo = 0xed0dd4bf62edd13f, .ex = 0, .sgn=0},
  {.hi = 0xf56d97473d446cda, .lo = 0x275a2bbb2bab6c8a, .ex = 0, .sgn=0},
  {.hi = 0xf5341c9f32bffeb9, .lo = 0x8da64484aaa0febc, .ex = 0, .sgn=0},
  {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=0},
  {.hi = 0xf4bf61b00b5982b7, .lo = 0x890ac4aafa6a37bf, .ex = 0, .sgn=0},
  {.hi = 0xf48421b0efbf939b, .lo = 0xf8f9d3b87d11fd52, .ex = 0, .sgn=0},
  {.hi = 0xf4484add6b01254b, .lo = 0x667e06866c07c369, .ex = 0, .sgn=0},
  {.hi = 0xf40bdd5a6688662f, .lo = 0x5019794a1f5896e5, .ex = 0, .sgn=0},
  {.hi = 0xf3ced94d28b2ce8a, .lo = 0x18ef535a7ffa7a3d, .ex = 0, .sgn=0},
  {.hi = 0xf3913edb54ba2242, .lo = 0x50f29b4b49f31c37, .ex = 0, .sgn=0},
  {.hi = 0xf3530e2aea9d3966, .lo = 0xd981acdcf6bc3e4, .ex = 0, .sgn=0},
  {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=0},
  {.hi = 0xf2d4eaa8233e997d, .lo = 0x431be53f92ece9e6, .ex = 0, .sgn=0},
  {.hi = 0xf294f82394ffe320, .lo = 0xebadcdbf915e8f6c, .ex = 0, .sgn=0},
  {.hi = 0xf2546ffc0e72f286, .lo = 0xaf0eed81e8c51e55, .ex = 0, .sgn=0},
  {.hi = 0xf21352595e0bf350, .lo = 0xe7112e89103cc0c7, .ex = 0, .sgn=0},
  {.hi = 0xf1d19f63ae7428a2, .lo = 0x844e6a35ddc2b713, .ex = 0, .sgn=0},
  {.hi = 0xf18f574386712643, .lo = 0x8f6bac72988088b0, .ex = 0, .sgn=0},
  {.hi = 0xf14c7a21c8cbd0f4, .lo = 0x2730081c758fb42b, .ex = 0, .sgn=0},
  {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=0},
  {.hi = 0xf0c5017ee336ca0f, .lo = 0xc4e557b119ef3185, .ex = 0, .sgn=0},
  {.hi = 0xf08066514c055f7e, .lo = 0x973ea9903ed5125f, .ex = 0, .sgn=0},
  {.hi = 0xf03b36c9407aa3e8, .lo = 0x992d39ec5c561d28, .ex = 0, .sgn=0},
  {.hi = 0xeff573116df1555d, .lo = 0x62aef7b55319d1d4, .ex = 0, .sgn=0},
  {.hi = 0xefaf1b54dd2cdf0f, .lo = 0xf03a18a5e16ab641, .ex = 0, .sgn=0},
  {.hi = 0xef682fbef23ecda6, .lo = 0x767c0e8ad33bc085, .ex = 0, .sgn=0},
  {.hi = 0xef20b07b6c6c0b37, .lo = 0xe2398bf0eeb28cde, .ex = 0, .sgn=0},
  {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=0},
  {.hi = 0xee8ff79c548acd0f, .lo = 0xa1d2c3d018a9279f, .ex = 0, .sgn=0},
  {.hi = 0xee46be5a0813016b, .lo = 0x7872773830d368be, .ex = 0, .sgn=0},
  {.hi = 0xedfcf21cabacd3b1, .lo = 0xfee6a1eebfa13b4a, .ex = 0, .sgn=0},
  {.hi = 0xedb29311c504d652, .lo = 0x11815196b9fbf5df, .ex = 0, .sgn=0},
  {.hi = 0xed67a1673455c601, .lo = 0x7289102076a125e5, .ex = 0, .sgn=0},
  {.hi = 0xed1c1d4b344c3d4f, .lo = 0xddffe98c4f8aa031, .ex = 0, .sgn=0},
  {.hi = 0xecd006ec59ea306f, .lo = 0xa8392eb238578ab0, .ex = 0, .sgn=0},
  {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=0},
  {.hi = 0xec3624222d227bd1, .lo = 0x278047ae3dd0889, .ex = 0, .sgn=0},
  {.hi = 0xebe85815c767cb00, .lo = 0x1e99ccb9adc62ca6, .ex = 0, .sgn=0},
  {.hi = 0xeb99fa84606ff5ff, .lo = 0xdae311e656e0661, .ex = 0, .sgn=0},
  {.hi = 0xeb4b0b9e4f345617, .lo = 0x39e39c6c2ab3655d, .ex = 0, .sgn=0},
  {.hi = 0xeafb8b944453f52f, .lo = 0x3383bbb5156bf1d7, .ex = 0, .sgn=0},
  {.hi = 0xeaab7a9749f584fe, .lo = 0x24db98ad3a0647a1, .ex = 0, .sgn=0},
  {.hi = 0xea5ad8d8c3a91f05, .lo = 0x4a0ca5ea449b1c83, .ex = 0, .sgn=0},
  {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=0},
  {.hi = 0xe9b7e3de5fdedc8b, .lo = 0xcd24d4bd1056c826, .ex = 0, .sgn=0},
  {.hi = 0xe9659107077cf60f, .lo = 0x89a92b199adfbafa, .ex = 0, .sgn=0},
  {.hi = 0xe912ae372d27045d, .lo = 0xacb1c26a06e5ae02, .ex = 0, .sgn=0},
  {.hi = 0xe8bf3ba1f1aedfbb, .lo = 0xf8972affb3d98e1f, .ex = 0, .sgn=0},
  {.hi = 0xe86b397ace95c46f, .lo = 0x9fec1e78c4376186, .ex = 0, .sgn=0},
  {.hi = 0xe816a7f595ec9232, .lo = 0xbfe8378abfb87b6f, .ex = 0, .sgn=0},
  {.hi = 0xe7c187467233d508, .lo = 0xdbfb0fe56c6f80fe, .ex = 0, .sgn=0},
  {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=0},
  {.hi = 0xe715993ccd02fe9c, .lo = 0xe2ba81b9ce96e02e, .ex = 0, .sgn=0},
  {.hi = 0xe6becc4c5997af06, .lo = 0x82fcedb4c6434d76, .ex = 0, .sgn=0},
  {.hi = 0xe667710616f4fc59, .lo = 0xdd2a3e32c3859960, .ex = 0, .sgn=0},
  {.hi = 0xe60f879fe7e2e1e5, .lo = 0x7613b68f6ab03130, .ex = 0, .sgn=0},
  {.hi = 0xe5b7105006d4c560, .lo = 0x9b695cd67c93bd79, .ex = 0, .sgn=0},
  {.hi = 0xe55e0b4d05c80388, .lo = 0x5a7c210a3a15e7ea, .ex = 0, .sgn=0},
  {.hi = 0xe50478cdce2246bc, .lo = 0xe1f5a58c80292554, .ex = 0, .sgn=0},
  {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=0},
  {.hi = 0xe44fac3814e09856, .lo = 0x20d63b5b9e3cd6ac, .ex = 0, .sgn=0},
  {.hi = 0xe3f4729119e798d9, .lo = 0x56992551ae074e99, .ex = 0, .sgn=0},
  {.hi = 0xe398ac4cf556b732, .lo = 0xd1197dc12c63176, .ex = 0, .sgn=0},
  {.hi = 0xe33c59a4439cd8ec, .lo = 0x36563e2ffad8351a, .ex = 0, .sgn=0},
  {.hi = 0xe2df7acff7c2cf83, .lo = 0xd6fe4dd22e60a4a2, .ex = 0, .sgn=0},
  {.hi = 0xe28210095b483751, .lo = 0xfd39138aa2d508ed, .ex = 0, .sgn=0},
  {.hi = 0xe224198a0e002123, .lo = 0xe0521df01a1be6f5, .ex = 0, .sgn=0},
  {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=0},
  {.hi = 0xe1668a498f1f892c, .lo = 0xe2f9d4600f4d0325, .ex = 0, .sgn=0},
  {.hi = 0xe106f1fd4b8d7c96, .lo = 0x6ba8a9d9ba877899, .ex = 0, .sgn=0},
  {.hi = 0xe0a6cee232f2bb9c, .lo = 0x6d6c98fe79817946, .ex = 0, .sgn=0},
  {.hi = 0xe046213392aa486c, .lo = 0x55ff6038a5197367, .ex = 0, .sgn=0},
  {.hi = 0xdfe4e92d0d8a37f5, .lo = 0x720588ff6547d884, .ex = 0, .sgn=0},
  {.hi = 0xdf83270a9bbee890, .lo = 0xab01350f013d78dd, .ex = 0, .sgn=0},
  {.hi = 0xdf20db088aa60404, .lo = 0x64a58b2f103485dd, .ex = 0, .sgn=0},
  {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=0},
  {.hi = 0xde5aa65869193805, .lo = 0x4248f15548f69ca, .ex = 0, .sgn=0},
  {.hi = 0xddf6be249c075037, .lo = 0xd597b10a01676659, .ex = 0, .sgn=0},
  {.hi = 0xdd924d05b620678a, .lo = 0x739c45b982193b5e, .ex = 0, .sgn=0},
  {.hi = 0xdd2d5339ac8692fd, .lo = 0x49c6e0ea76cbcaac, .ex = 0, .sgn=0},
  {.hi = 0xdcc7d0fec8aaf2aa, .lo = 0xb2069fd0b482b4e8, .ex = 0, .sgn=0},
  {.hi = 0xdc61c693a82745d5, .lo = 0xaca8017e375b64e5, .ex = 0, .sgn=0},
  {.hi = 0xdbfb34373c974b0e, .lo = 0xccb7fd40d543f4a1, .ex = 0, .sgn=0},
  {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=0},
  {.hi = 0xdb2c78a7ede238a9, .lo = 0x5a98479cbef2ecbc, .ex = 0, .sgn=0},
  {.hi = 0xdac44ff490a02710, .lo = 0x5b267c1bcff0ab62, .ex = 0, .sgn=0},
  {.hi = 0xda5ba04ef3c929f4, .lo = 0xe257bde73d83dc1a, .ex = 0, .sgn=0},
  {.hi = 0xd9f269f7aab88c29, .lo = 0x28e81dcb6dab91ac, .ex = 0, .sgn=0},
  {.hi = 0xd988ad2f9bdf9bbb, .lo = 0xc4e4dc69fc2fff6f, .ex = 0, .sgn=0},
  {.hi = 0xd91e6a38009da15a, .lo = 0x1bb35ad6d2e74b67, .ex = 0, .sgn=0},
  {.hi = 0xd8b3a1526517a48b, .lo = 0x1ed1a8ff78f1b632, .ex = 0, .sgn=0},
  {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=0},
  {.hi = 0xd7dc7ec4fabdb011, .lo = 0xced12d2899b803db, .ex = 0, .sgn=0},
  {.hi = 0xd77025a1e0a39d8b, .lo = 0xcb78e80e67ba1b8, .ex = 0, .sgn=0},
  {.hi = 0xd703479a2f6776cc, .lo = 0x6cb3bfd65b38562b, .ex = 0, .sgn=0},
  {.hi = 0xd695e4f10ea88570, .lo = 0x83f082b570611d7, .ex = 0, .sgn=0},
  {.hi = 0xd627fde9f7d63e7e, .lo = 0x7afbefc05e9f7d99, .ex = 0, .sgn=0},
  {.hi = 0xd5b992c8b606a351, .lo = 0x7190b755535d4f18, .ex = 0, .sgn=0},
  {.hi = 0xd54aa3d165cc7018, .lo = 0x7d00ae97abaa4096, .ex = 0, .sgn=0},
  {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=0},
  {.hi = 0xd46b3b72a2d68fc9, .lo = 0xdc4663a3168698d2, .ex = 0, .sgn=0},
  {.hi = 0xd3fac294ff34e4d0, .lo = 0xb77d4f6bd0ee8591, .ex = 0, .sgn=0},
  {.hi = 0xd389c6f4eb07a41c, .lo = 0xa8faac741a6394dc, .ex = 0, .sgn=0},
  {.hi = 0xd31848d817d70e16, .lo = 0xeeeaddb72f00e0dd, .ex = 0, .sgn=0},
  {.hi = 0xd2a6488487a91918, .lo = 0x4300fd1c1ce507e5, .ex = 0, .sgn=0},
  {.hi = 0xd233c6408cd64236, .lo = 0x981ba7e42537275f, .ex = 0, .sgn=0},
  {.hi = 0xd1c0c252c9de2c86, .lo = 0xda7485a5aeffeb4c, .ex = 0, .sgn=0},
  {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=0},
  {.hi = 0xd0d93696053af098, .lo = 0x77a18eb13d2ecde5, .ex = 0, .sgn=0},
  {.hi = 0xd064af55d7c9b43e, .lo = 0x6b8a685f6cb61c21, .ex = 0, .sgn=0},
  {.hi = 0xcfefa7898a4ef23c, .lo = 0xdaf200dd81212d10, .ex = 0, .sgn=0},
  {.hi = 0xcf7a1f794d7ca1b1, .lo = 0xdfcb60445c1bf973, .ex = 0, .sgn=0},
  {.hi = 0xcf04176da12390ac, .lo = 0x4d27090f10c454e, .ex = 0, .sgn=0},
  {.hi = 0xce8d8faf5406ab8b, .lo = 0xf5babff66def7892, .ex = 0, .sgn=0},
  {.hi = 0xce16888783ae13b3, .lo = 0x93e391861a034684, .ex = 0, .sgn=0},
  {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=0},
  {.hi = 0xcd26fd2158358e7d, .lo = 0x649474e36b8db9d3, .ex = 0, .sgn=0},
  {.hi = 0xccae7976c0691177, .lo = 0x83e907fbd7aaf0b0, .ex = 0, .sgn=0},
  {.hi = 0xcc35778a2bac9ca1, .lo = 0xf839ce18e08bfb50, .ex = 0, .sgn=0},
  {.hi = 0xcbbbf7a63eba0dd5, .lo = 0x70cbb7f3343451be, .ex = 0, .sgn=0},
  {.hi = 0xcb41fa15ebff0777, .lo = 0x2293661be51140ab, .ex = 0, .sgn=0},
  {.hi = 0xcac77f24736eb553, .lo = 0xd9944be1631846d8, .ex = 0, .sgn=0},
  {.hi = 0xca4c871d625361a9, .lo = 0x5328edeb3e6784de, .ex = 0, .sgn=0},
  {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=0},
  {.hi = 0xc95520fe2d40a74b, .lo = 0x83b0e96e1249c2b0, .ex = 0, .sgn=0},
  {.hi = 0xc8d8b37ea4ed0f62, .lo = 0xb562c00b34ee771, .ex = 0, .sgn=0},
  {.hi = 0xc85bca1abaf7f0a7, .lo = 0x65862939b83382e0, .ex = 0, .sgn=0},
  {.hi = 0xc7de651f7ca06749, .lo = 0x2b31bc86877fd2c, .ex = 0, .sgn=0},
  {.hi = 0xc76084da43624634, .lo = 0xd5c149509e9059f1, .ex = 0, .sgn=0},
  {.hi = 0xc6e22998b4c6608e, .lo = 0xcfe6c1b1a6b4e2a4, .ex = 0, .sgn=0},
  {.hi = 0xc66353a8c232a43c, .lo = 0xe993503baf5afb41, .ex = 0, .sgn=0},
  {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=0},
  {.hi = 0xc56438f6f0ec3cca, .lo = 0xab4906075507e74, .ex = 0, .sgn=0},
  {.hi = 0xc4e3f4d26ea553b6, .lo = 0xdd40950cf1ed92fa, .ex = 0, .sgn=0},
  {.hi = 0xc463373a40dd06a3, .lo = 0x9dd768f30ca8e85c, .ex = 0, .sgn=0},
  {.hi = 0xc3e2007dd175f5a4, .lo = 0xa87e78136665cdb2, .ex = 0, .sgn=0},
  {.hi = 0xc36050ecd50ca830, .lo = 0x8ac9e1386e4cbabb, .ex = 0, .sgn=0},
  {.hi = 0xc2de28d74ac6628b, .lo = 0x74c8f010d986a9e0, .ex = 0, .sgn=0},
  {.hi = 0xc25b888d7c1fcd38, .lo = 0xb7041e9bc8c18b0d, .ex = 0, .sgn=0},
  {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=0},
  {.hi = 0xc154e09faa2ff69a, .lo = 0x17858573216e0a22, .ex = 0, .sgn=0},
  {.hi = 0xc0d0d99dabd65d44, .lo = 0x2bda5328933c854a, .ex = 0, .sgn=0},
  {.hi = 0xc04c5bab7297d322, .lo = 0x6dd06968e0ed1957, .ex = 0, .sgn=0},
  {.hi = 0xbfc7671ab8bb84c6, .lo = 0xe4e62d86dd136e78, .ex = 0, .sgn=0},
  {.hi = 0xbf41fc3d81b430db, .lo = 0xd46655d6b012455, .ex = 0, .sgn=0},
  {.hi = 0xbebc1b6619ed9116, .lo = 0x2715ef03f8543355, .ex = 0, .sgn=0},
  {.hi = 0xbe35c4e716999630, .lo = 0x29d7f7b67d43b177, .ex = 0, .sgn=0},
  {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=0},
  {.hi = 0xbd27b83dfcbe9279, .lo = 0x2ea36923d5d8e213, .ex = 0, .sgn=0},
  {.hi = 0xbca002ba7aaf25ea, .lo = 0x4a48496734be336d, .ex = 0, .sgn=0},
  {.hi = 0xbc17d8dc859ad583, .lo = 0x727c405ffc73af56, .ex = 0, .sgn=0},
  {.hi = 0xbb8f3af81b93095c, .lo = 0xfce8d84068e825b6, .ex = 0, .sgn=0},
  {.hi = 0xbb062961823b1ddc, .lo = 0x5120e35e1c1a250c, .ex = 0, .sgn=0},
  {.hi = 0xba7ca46d46946802, .lo = 0x33201477347447d8, .ex = 0, .sgn=0},
  {.hi = 0xb9f2ac703cca0db3, .lo = 0x39db32d014440024, .ex = 0, .sgn=0},
  {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=0},
  {.hi = 0xb8dd64b0720df647, .lo = 0xa726f4f0828585c9, .ex = 0, .sgn=0},
  {.hi = 0xb8521598bb6bce26, .lo = 0x1c041d1ea5fb3fdb, .ex = 0, .sgn=0},
  {.hi = 0xb7c654ce4adba9f2, .lo = 0x2e7a35723f3ed035, .ex = 0, .sgn=0},
  {.hi = 0xb73a22a755457448, .lo = 0x7f86f63bb23f496a, .ex = 0, .sgn=0},
  {.hi = 0xb6ad7f7a557e64f2, .lo = 0xeb2d28ef943dc88c, .ex = 0, .sgn=0},
  {.hi = 0xb6206b9e0c13a892, .lo = 0xea7c015f12b987f7, .ex = 0, .sgn=0},
  {.hi = 0xb592e7697f14dd4a, .lo = 0x737dd2824b608d13, .ex = 0, .sgn=0},
};

/* The following is a degree-7 polynomial with odd coefficients
   approximating sin2pi(x) for 0 <= x < 2^-11 with relative error 2^-77.307.
   Generated with sin_fast.sollya. */
static const  double PSfast[] = {
  0x1.921fb54442d18p+2, 0x1.1a62645458ee1p-52, // degree 1 (h+l)
  -0x1.4abbce625be53p5,                        // degree 3
  0x1.466bc678d8e3fp6,                         // degree 5
  -0x1.33155a7aff959p6,                        // degree 7
};

/* The following is a degree-6 polynomial with even coefficients
   approximating cos2pi(x) for 0 <= x < 2^-11 with relative error 2^-75.189.
   Generated with cos_fast.sollya. */
static const  double PCfast[] = {
  0x1p+0, -0x1.9249c1ep-77,                    // degree 0
  -0x1.3bd3cc9be45dep4,                        // degree 2
  0x1.03c1f080ad7f9p6,                         // degree 4
  -0x1.55a5c19e443dcp6,                        // degree 6
};

/* The following is a degree-11 polynomial with odd coefficients
   approximating sin2pi(x) for 0 <= x < 2^-11 with relative error 2^-127.75.
   Generated with sin_accurate.sollya. */
static const dint64_t PS[] = {
  {.hi = 0xc90fdaa22168c234, .lo = 0xc4c6628b80dc1cd1, .ex = 3, .sgn=0}, // 1
  {.hi = 0xa55de7312df295f5, .lo = 0x5dc72f712aa57db4, .ex = 6, .sgn=1}, // 3
  {.hi = 0xa335e33bad570e92, .lo = 0x3f33be0021aa54d2, .ex = 7, .sgn=0}, // 5
  {.hi = 0x9969667315ec2d9d, .lo = 0xe59d6ab8509a2025, .ex = 7, .sgn=1}, // 7
  {.hi = 0xa83c1a43bf1c6485, .lo = 0x7d5f8f76fa7d74ed, .ex = 6, .sgn=0}, // 9
  {.hi = 0xf16ab2898eae62f9, .lo = 0xa7f0339113b8b3c5, .ex = 4, .sgn=1}, // 11
};

/* The following is a degree-10 polynomial with even coefficients
   approximating cos2pi(x) for 0 <= x < 2^-11 with relative error 2^-137.246.
   Generated with cos_accurate.sollya. */
static const dint64_t PC[] = {
  {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=0}, // degree 0
  {.hi = 0x9de9e64df22ef2d2, .lo = 0x56e26cd9808c1949, .ex = 5, .sgn=1}, // 2
  {.hi = 0x81e0f840dad61d9a, .lo = 0x9980f00630cb655e, .ex = 7, .sgn=0}, // 4
  {.hi = 0xaae9e3f1e5ffcfe2, .lo = 0xa508509534006249, .ex = 7, .sgn=1}, // 6
  {.hi = 0xf0fa83448dd1e094, .lo = 0xe0603ce7044eeba, .ex = 6, .sgn=0},  // 8
  {.hi = 0xd368f6f4207cfe49, .lo = 0xec63157807ebffa, .ex = 5, .sgn=1},  // 10
};

/* Table generated with ./buildSC 11 using accompanying buildSC.c.
   For each i, 0 <= i < 256, xi=i/2^11+SC[i][0], with
   SC[i][1] and SC[i][2] approximating sin2pi(xi) and cos2pi(xi)
   respectively, both with 53+11 bits of accuracy. */
static const double SC[256][3] = {
   {0x0p+0, 0x0p+0, 0x1p+0}, /* 0 */
   {-0x1.a834cp-45, 0x1.921f8bec23b6p-9, 0x1.ffff621621d0ap-1}, /* 1 */
   {-0x1.cef8p-44, 0x1.921f0fe5ba394p-8, 0x1.fffd8858e8ab6p-1}, /* 2 */
   {-0x1.a8d5p-46, 0x1.2d96b0e4f495dp-7, 0x1.fffa72c978c5bp-1}, /* 3 */
   {-0x1.e8823p-41, 0x1.921d1fcaed2e2p-7, 0x1.fff62169b9536p-1}, /* 4 */
   {0x1.e9d64p-42, 0x1.f6a296ad1a28ap-7, 0x1.fff0943c53a57p-1}, /* 5 */
   {-0x1.db31cp-43, 0x1.2d936bbdd3a6p-6, 0x1.ffe9cb44b527dp-1}, /* 6 */
   {-0x1.4030ep-42, 0x1.5fd4d21f2d6cbp-6, 0x1.ffe1c6870ccd1p-1}, /* 7 */
   {0x1.ec708p-43, 0x1.92155f7a97112p-6, 0x1.ffd886084cbddp-1}, /* 8 */
   {0x1.ccef8p-40, 0x1.c454f4d127745p-6, 0x1.ffce09ce29c7ap-1}, /* 9 */
   {0x1.31659cp-38, 0x1.f69373249ae67p-6, 0x1.ffc251df1b68ap-1}, /* 10 */
   {-0x1.0a03cp-38, 0x1.14685db0e8dc1p-5, 0x1.ffb55e42619e1p-1}, /* 11 */
   {-0x1.6d229p-39, 0x1.2d8657570832dp-5, 0x1.ffa72efff0c79p-1}, /* 12 */
   {0x1.8efe09p-36, 0x1.46a3971317fd8p-5, 0x1.ff97c4207f829p-1}, /* 13 */
   {0x1.7d1b7p-40, 0x1.5fc00d2a37dfbp-5, 0x1.ff871dadb7505p-1}, /* 14 */
   {-0x1.6fb75p-40, 0x1.78dbaa5753e8fp-5, 0x1.ff753bb1b9eadp-1}, /* 15 */
   {-0x1.860bdp-39, 0x1.91f65f0e798f5p-5, 0x1.ff621e3798b8ep-1}, /* 16 */
   {-0x1.1f85ap-39, 0x1.ab101bd4352b3p-5, 0x1.ff4dc54b1d65ep-1}, /* 17 */
   {-0x1.950eep-40, 0x1.c428d12acfd89p-5, 0x1.ff3830f8d68ebp-1}, /* 18 */
   {0x1.19ac5p-38, 0x1.dd406f9b7c52cp-5, 0x1.ff21614e0fe5fp-1}, /* 19 */
   {0x1.80b7p-40, 0x1.f656e7a0afa3fp-5, 0x1.ff095658e5f26p-1}, /* 20 */
   {0x1.bb805p-39, 0x1.07b614e5bea07p-4, 0x1.fef01028234b7p-1}, /* 21 */
   {-0x1.04b098p-38, 0x1.1440134bd80c4p-4, 0x1.fed58ecb6abp-1}, /* 22 */
   {0x1.d4324p-40, 0x1.20c9674f8babep-4, 0x1.feb9d25302722p-1}, /* 23 */
   {-0x1.272734p-36, 0x1.2d520925a8254p-4, 0x1.fe9cdad029914p-1}, /* 24 */
   {-0x1.cf2eap-40, 0x1.39d9f12ba4ce5p-4, 0x1.fe7ea8548494p-1}, /* 25 */
   {-0x1.42bfap-39, 0x1.46611791755b4p-4, 0x1.fe5f3af2e61a6p-1}, /* 26 */
   {0x1.c7aa8p-42, 0x1.52e774a501658p-4, 0x1.fe3e92be9d11fp-1}, /* 27 */
   {0x1.8f706p-39, 0x1.5f6d00aae2d16p-4, 0x1.fe1cafcbd2534p-1}, /* 28 */
   {0x1.0f9a6p-40, 0x1.6bf1b3e8054f1p-4, 0x1.fdf9922f72013p-1}, /* 29 */
   {0x1.6e2e4cp-36, 0x1.787586aec8bd5p-4, 0x1.fdd539ff04d69p-1}, /* 30 */
   {-0x1.7f29p-42, 0x1.84f8712bed97bp-4, 0x1.fdafa75145ab1p-1}, /* 31 */
   {-0x1.79566p-40, 0x1.917a6bc207cb5p-4, 0x1.fd88da3d14232p-1}, /* 32 */
   {-0x1.481a8p-38, 0x1.9dfb6eb049c72p-4, 0x1.fd60d2da7c4ccp-1}, /* 33 */
   {0x1.c3dp-43, 0x1.aa7b7244abcfp-4, 0x1.fd379142206e6p-1}, /* 34 */
   {0x1.a74f4p-38, 0x1.b6fa6ec6247f4p-4, 0x1.fd0d158d7d201p-1}, /* 35 */
   {-0x1.014fp-41, 0x1.c3785c79b9f66p-4, 0x1.fce15fd6db19ep-1}, /* 36 */
   {-0x1.601ap-42, 0x1.cff533b2e583cp-4, 0x1.fcb4703914b29p-1}, /* 37 */
   {0x1.1f77bp-38, 0x1.dc70ecbcaa797p-4, 0x1.fc8646cfe4e15p-1}, /* 38 */
   {-0x1.cbcdp-40, 0x1.e8eb7fdd975e3p-4, 0x1.fc56e3b7dc611p-1}, /* 39 */
   {-0x1.2dfacp-36, 0x1.f564e5633c113p-4, 0x1.fc26470e37059p-1}, /* 40 */
   {-0x1.0548ap-39, 0x1.00ee8ad695ba2p-3, 0x1.fbf470f0ac105p-1}, /* 41 */
   {-0x1.f4fedp-38, 0x1.072a047a21f9dp-3, 0x1.fbc1617e50bc5p-1}, /* 42 */
   {0x1.75bf8p-40, 0x1.0d64dbcb6f37cp-3, 0x1.fb8d18d66871ap-1}, /* 43 */
   {-0x1.c8548p-40, 0x1.139f0ced568e7p-3, 0x1.fb5797196077dp-1}, /* 44 */
   {-0x1.e70ccp-39, 0x1.19d8940b24dccp-3, 0x1.fb20dc6823e96p-1}, /* 45 */
   {-0x1.67e3bp-38, 0x1.20116d4dafe44p-3, 0x1.fae8e8e476ec3p-1}, /* 46 */
   {-0x1.b2992p-39, 0x1.264994df2a5b7p-3, 0x1.faafbcb0d5ff4p-1}, /* 47 */
   {0x1.1bd37p-38, 0x1.2c8106e9c2949p-3, 0x1.fa7557f082241p-1}, /* 48 */
   {-0x1.fde3bp-37, 0x1.32b7bf9139844p-3, 0x1.fa39bac7bf75fp-1}, /* 49 */
   {-0x1.08e9ep-37, 0x1.38edbb0b3d941p-3, 0x1.f9fce55aeb11dp-1}, /* 50 */
   {-0x1.230a2p-39, 0x1.3f22f57d43a42p-3, 0x1.f9bed7cfc2566p-1}, /* 51 */
   {-0x1.d28a8p-40, 0x1.45576b1239748p-3, 0x1.f97f924c943dp-1}, /* 52 */
   {-0x1.13ecap-38, 0x1.4b8b17f6c9ce7p-3, 0x1.f93f14f86385cp-1}, /* 53 */
   {0x1.7f018p-41, 0x1.51bdf859a1755p-3, 0x1.f8fd5ffae290dp-1}, /* 54 */
   {-0x1.3e262p-38, 0x1.57f00864566acp-3, 0x1.f8ba737cbf352p-1}, /* 55 */
   {-0x1.61388p-41, 0x1.5e21444891d1bp-3, 0x1.f8764fa71636p-1}, /* 56 */
   {0x1.1ecd5p-38, 0x1.6451a832b6025p-3, 0x1.f830f4a402949p-1}, /* 57 */
   {0x1.f5b61p-37, 0x1.6a8130526c4f9p-3, 0x1.f7ea629e40f74p-1}, /* 58 */
   {-0x1.e8064p-39, 0x1.70afd8cfcfcbfp-3, 0x1.f7a299c1abc2bp-1}, /* 59 */
   {0x1.ef77cp-40, 0x1.76dd9de56b974p-3, 0x1.f7599a3a0d93dp-1}, /* 60 */
   {0x1.ada62p-39, 0x1.7d0a7bbdd2789p-3, 0x1.f70f6434b0126p-1}, /* 61 */
   {-0x1.d726p-43, 0x1.83366e89baf17p-3, 0x1.f6c3f7df5c476p-1}, /* 62 */
   {0x1.1b401p-37, 0x1.8961727df625p-3, 0x1.f67755686e715p-1}, /* 63 */
   {-0x1.250cp-40, 0x1.8f8b83c661f19p-3, 0x1.f6297cff78997p-1}, /* 64 */
   {0x1.d4e02p-37, 0x1.95b49e9e34991p-3, 0x1.f5da6ed4120c2p-1}, /* 65 */
   {-0x1.0db04p-38, 0x1.9bdcbf2cf4ba2p-3, 0x1.f58a2b17948e7p-1}, /* 66 */
   {0x1.474318p-34, 0x1.a203e1c13d097p-3, 0x1.f538b1fa20f7p-1}, /* 67 */
   {-0x1.9634ap-38, 0x1.a82a0259c8277p-3, 0x1.f4e603b0c377ap-1}, /* 68 */
   {0x1.9ef77p-36, 0x1.ae4f1d643628cp-3, 0x1.f492206b8640dp-1}, /* 69 */
   {-0x1.ec44p-41, 0x1.b4732ef3a73a6p-3, 0x1.f43d085ffbc0fp-1}, /* 70 */
   {-0x1.c09888p-36, 0x1.ba963349b5dap-3, 0x1.f3e6bbc207ea3p-1}, /* 71 */
   {-0x1.2bbe1p-37, 0x1.c0b826a619915p-3, 0x1.f38f3ac66822dp-1}, /* 72 */
   {0x1.2274p-38, 0x1.c6d90536b5ba6p-3, 0x1.f33685a39e448p-1}, /* 73 */
   {0x1.3fc51p-36, 0x1.ccf8cb34fdf74p-3, 0x1.f2dc9c9051226p-1}, /* 74 */
   {-0x1.05094p-39, 0x1.d31774d268105p-3, 0x1.f2817fc466752p-1}, /* 75 */
   {0x1.64a38p-39, 0x1.d934fe54dc742p-3, 0x1.f2252f775b969p-1}, /* 76 */
   {0x1.51da3p-37, 0x1.df5163f2148edp-3, 0x1.f1c7abe265636p-1}, /* 77 */
   {0x1.8aabp-39, 0x1.e56ca1e198337p-3, 0x1.f168f53f68d6cp-1}, /* 78 */
   {-0x1.2c34ep-38, 0x1.eb86b461f9509p-3, 0x1.f1090bc8a71b9p-1}, /* 79 */
   {0x1.88c1ep-38, 0x1.f19f97b3412bap-3, 0x1.f0a7efb910508p-1}, /* 80 */
   {0x1.11cf1p-37, 0x1.f7b7480d7462ep-3, 0x1.f045a14cdcc68p-1}, /* 81 */
   {-0x1.1ea75p-36, 0x1.fdcdc1aa96acp-3, 0x1.efe220c0f169cp-1}, /* 82 */
   {-0x1.7bcc98p-36, 0x1.01f180695e83ap-2, 0x1.ef7d6e52155fbp-1}, /* 83 */
   {-0x1.590bp-41, 0x1.04fb80e36f7a2p-2, 0x1.ef178a3e4964bp-1}, /* 84 */
   {-0x1.2b2cp-41, 0x1.0804e05ea8318p-2, 0x1.eeb074c50c38fp-1}, /* 85 */
   {-0x1.fb9f4p-38, 0x1.0b0d9cfcfd474p-2, 0x1.ee482e25c3da2p-1}, /* 86 */
   {0x1.042ccp-37, 0x1.0e15b4e239b7cp-2, 0x1.eddeb6a05d726p-1}, /* 87 */
   {0x1.84f19p-37, 0x1.111d262c45d07p-2, 0x1.ed740e765bd98p-1}, /* 88 */
   {-0x1.22a4ep-38, 0x1.1423eefbfb4f4p-2, 0x1.ed0835e9a8644p-1}, /* 89 */
   {-0x1.d144p-38, 0x1.172a0d76b54d9p-2, 0x1.ec9b2d3c54ep-1}, /* 90 */
   {0x1.529b1p-37, 0x1.1a2f7fbf8ec8ep-2, 0x1.ec2cf4b18ac69p-1}, /* 91 */
   {0x1.7ed7ap-38, 0x1.1d3443f55e187p-2, 0x1.ebbd8c8ddbc78p-1}, /* 92 */
   {-0x1.30aeb8p-36, 0x1.2038583ba73d1p-2, 0x1.eb4cf515fbdbap-1}, /* 93 */
   {-0x1.f9486p-38, 0x1.233bbabb7d7c4p-2, 0x1.eadb2e8e96c05p-1}, /* 94 */
   {0x1.6b72p-39, 0x1.263e699599a5fp-2, 0x1.ea68393e5b3f4p-1}, /* 95 */
   {-0x1.20e4cp-37, 0x1.294062ec80dp-2, 0x1.e9f4156c83cc5p-1}, /* 96 */
   {-0x1.4ec28p-37, 0x1.2c41a4e858f51p-2, 0x1.e97ec3603d3efp-1}, /* 97 */
   {-0x1.64f65p-37, 0x1.2f422dadb4708p-2, 0x1.e90843620902cp-1}, /* 98 */
   {0x1.5a13ccp-35, 0x1.3241fb6799231p-2, 0x1.e89095ba336c9p-1}, /* 99 */
   {-0x1.c3792p-38, 0x1.35410c2d6f116p-2, 0x1.e817bab4e7d66p-1}, /* 100 */
   {-0x1.d2d84p-39, 0x1.383f5e34e41e1p-2, 0x1.e79db29a5f5f6p-1}, /* 101 */
   {-0x1.38398p-40, 0x1.3b3cefa024218p-2, 0x1.e7227db6ae2c2p-1}, /* 102 */
   {-0x1.492d1p-36, 0x1.3e39be9500afbp-2, 0x1.e6a61c5625928p-1}, /* 103 */
   {-0x1.e62bcp-39, 0x1.4135c9411bbbep-2, 0x1.e6288ec49d09fp-1}, /* 104 */
   {-0x1.b99e9p-37, 0x1.44310dc74a6d9p-2, 0x1.e5a9d5507d64dp-1}, /* 105 */
   {0x1.5c2b08p-36, 0x1.472b8a5777419p-2, 0x1.e529f046d2a19p-1}, /* 106 */
   {-0x1.73f78p-39, 0x1.4a253d11730c8p-2, 0x1.e4a8dff828abfp-1}, /* 107 */
   {-0x1.5020cp-39, 0x1.4d1e2427500ep-2, 0x1.e426a4b2c6d41p-1}, /* 108 */
   {-0x1.14bcp-41, 0x1.50163dc18a2f9p-2, 0x1.e3a33ec75f23p-1}, /* 109 */
   {-0x1.12634p-37, 0x1.530d880a28693p-2, 0x1.e31eae87308fbp-1}, /* 110 */
   {-0x1.1026cp-38, 0x1.5604012ee1bedp-2, 0x1.e298f443a370bp-1}, /* 111 */
   {-0x1.00ba4p-39, 0x1.58f9a75a8287ap-2, 0x1.e212104f70ecp-1}, /* 112 */
   {0x1.c6843p-36, 0x1.5bee78bc7ab57p-2, 0x1.e18a02fd4d22cp-1}, /* 113 */
   {0x1.21f708p-36, 0x1.5ee2737b964fep-2, 0x1.e100cca24a015p-1}, /* 114 */
   {0x1.b06eap-37, 0x1.61d595c9cad64p-2, 0x1.e0766d924647ap-1}, /* 115 */
   {-0x1.b333p-40, 0x1.64c7ddd3ca7p-2, 0x1.dfeae622e3542p-1}, /* 116 */
   {0x1.65b64cp-35, 0x1.67b949cef2662p-2, 0x1.df5e36a8f4f63p-1}, /* 117 */
   {0x1.b12cp-43, 0x1.6aa9d7dc7cda1p-2, 0x1.ded05f7de38cap-1}, /* 118 */
   {0x1.1480bp-36, 0x1.6d99863a367f9p-2, 0x1.de4160f68b4f4p-1}, /* 119 */
   {0x1.93ff9p-37, 0x1.70885310cc63bp-2, 0x1.ddb13b6c930aep-1}, /* 120 */
   {-0x1.40e94p-38, 0x1.73763c91eb992p-2, 0x1.dd1fef38bff12p-1}, /* 121 */
   {0x1.c2dbep-37, 0x1.766340f38b25fp-2, 0x1.dc8d7cb3cf6a8p-1}, /* 122 */
   {0x1.f1de8p-40, 0x1.794f5e616b6b5p-2, 0x1.dbf9e4394e586p-1}, /* 123 */
   {0x1.0ce78p-40, 0x1.7c3a9311f5519p-2, 0x1.db6526238522fp-1}, /* 124 */
   {0x1.39238p-37, 0x1.7f24dd3818314p-2, 0x1.dacf42ce3aa8dp-1}, /* 125 */
   {-0x1.163d6p-38, 0x1.820e3b0485789p-2, 0x1.da383a967d314p-1}, /* 126 */
   {0x1.0ebcep-37, 0x1.84f6aaaffdb71p-2, 0x1.d9a00dd88b71ep-1}, /* 127 */
   {-0x1.30f98p-35, 0x1.87de2a6775698p-2, 0x1.d906bcf3e027cp-1}, /* 128 */
   {0x1.7fcb88p-35, 0x1.8ac4b871b75bp-2, 0x1.d86c484371db1p-1}, /* 129 */
   {-0x1.3cd02p-37, 0x1.8daa52eba4ff7p-2, 0x1.d7d0b02bbf203p-1}, /* 130 */
   {-0x1.7b168p-38, 0x1.908ef81e6ebbbp-2, 0x1.d733f508ddfdbp-1}, /* 131 */
   {0x1.b274e4p-34, 0x1.9372a6459635fp-2, 0x1.d696173a84a7ap-1}, /* 132 */
   {0x1.99f79p-34, 0x1.96555b83f5b46p-2, 0x1.d5f717268995bp-1}, /* 133 */
   {0x1.3a5ep-38, 0x1.993716148d0a5p-2, 0x1.d556f52e7b404p-1}, /* 134 */
   {0x1.06e368p-34, 0x1.9c17d446c7bcep-2, 0x1.d4b5b1b03af49p-1}, /* 135 */
   {0x1.79c6ep-36, 0x1.9ef7943cad591p-2, 0x1.d4134d146456ep-1}, /* 136 */
   {-0x1.014bcp-38, 0x1.a1d6543af46d2p-2, 0x1.d36fc7bcd45bp-1}, /* 137 */
   {-0x1.38e858p-34, 0x1.a4b41276e9a72p-2, 0x1.d2cb220fa2d8cp-1}, /* 138 */
   {-0x1.6067cp-34, 0x1.a790cd35df451p-2, 0x1.d2255c70243fep-1}, /* 139 */
   {-0x1.47f4bp-36, 0x1.aa6c82b4ffa09p-2, 0x1.d17e77444ea0ep-1}, /* 140 */
   {-0x1.b3e01p-35, 0x1.ad473120f28b9p-2, 0x1.d0d672f6bc322p-1}, /* 141 */
   {0x1.e2bcp-39, 0x1.b020d6c849efp-2, 0x1.d02d4feb17d94p-1}, /* 142 */
   {-0x1.75a12p-35, 0x1.b2f971d70af6p-2, 0x1.cf830e8dddb4dp-1}, /* 143 */
   {-0x1.113d38p-34, 0x1.b5d1009805d13p-2, 0x1.ced7af453b7b4p-1}, /* 144 */
   {-0x1.5e06p-36, 0x1.b8a7814de51b7p-2, 0x1.ce2b327a1052dp-1}, /* 145 */
   {0x1.6a6bp-40, 0x1.bb7cf2306be21p-2, 0x1.cd7d9898ab7afp-1}, /* 146 */
   {-0x1.fcfp-38, 0x1.be51517f482d3p-2, 0x1.cccee20c59781p-1}, /* 147 */
   {0x1.3ff588p-35, 0x1.c1249d8399431p-2, 0x1.cc1f0f3ef354p-1}, /* 148 */
   {0x1.23cfdp-35, 0x1.c3f6d47599b2cp-2, 0x1.cb6e209f43591p-1}, /* 149 */
   {0x1.710cp-37, 0x1.c6c7f49a73b28p-2, 0x1.cabc1699cb32dp-1}, /* 150 */
   {0x1.1c364p-38, 0x1.c997fc38c9112p-2, 0x1.ca08f19b83549p-1}, /* 151 */
   {-0x1.ad9d4p-38, 0x1.cc66e9928593fp-2, 0x1.c954b213670c9p-1}, /* 152 */
   {-0x1.c8cb5p-36, 0x1.cf34baeb9ce5cp-2, 0x1.c89f5870cc0acp-1}, /* 153 */
   {-0x1.cc923p-36, 0x1.d2016e8c1980bp-2, 0x1.c7e8e522d86d7p-1}, /* 154 */
   {0x1.06d23p-36, 0x1.d4cd02bbf5125p-2, 0x1.c73158998c2d1p-1}, /* 155 */
   {-0x1.b56c88p-35, 0x1.d79775b3aa6b7p-2, 0x1.c678b349c3a9dp-1}, /* 156 */
   {0x1.7b8f5p-36, 0x1.da60c5d1b16dcp-2, 0x1.c5bef59f656bfp-1}, /* 157 */
   {-0x1.852adp-36, 0x1.dd28f145ffe49p-2, 0x1.c504201344fd9p-1}, /* 158 */
   {-0x1.11d7p-39, 0x1.dfeff66a649eap-2, 0x1.c44833142899fp-1}, /* 159 */
   {0x1.238dep-36, 0x1.e2b5d3820348fp-2, 0x1.c38b2f179fe99p-1}, /* 160 */
   {0x1.5de1ep-36, 0x1.e57a86d5b1691p-2, 0x1.c2cd14929bf71p-1}, /* 161 */
   {-0x1.9898ap-37, 0x1.e83e0eae6afbap-2, 0x1.c20de3fae39c4p-1}, /* 162 */
   {-0x1.a22398p-35, 0x1.eb00695aa49eap-2, 0x1.c14d9dc5a0d55p-1}, /* 163 */
   {0x1.e45e28p-35, 0x1.edc195342ca7fp-2, 0x1.c08c4265b676dp-1}, /* 164 */
   {0x1.9958dp-35, 0x1.f081907064381p-2, 0x1.bfc9d258e34fp-1}, /* 165 */
   {-0x1.420a8p-39, 0x1.f3405963c5d13p-2, 0x1.bf064e1546e7fp-1}, /* 166 */
   {0x1.33ceap-36, 0x1.f5fdee6712454p-2, 0x1.be41b6109ec8cp-1}, /* 167 */
   {-0x1.5d514p-37, 0x1.f8ba4dbe9af59p-2, 0x1.bd7c0ac73cefdp-1}, /* 168 */
   {0x1.eb1e2p-37, 0x1.fb7575c39c34fp-2, 0x1.bcb54cb0729f2p-1}, /* 169 */
   {-0x1.7091p-36, 0x1.fe2f64bc7b28ap-2, 0x1.bbed7c49c847ap-1}, /* 170 */
   {-0x1.902ap-41, 0x1.00740c82afadep-1, 0x1.bb249a0b712c3p-1}, /* 171 */
   {-0x1.c359p-35, 0x1.01cfc8725f626p-1, 0x1.ba5aa674be0c3p-1}, /* 172 */
   {-0x1.be042p-36, 0x1.032ae55dadbdep-1, 0x1.b98fa1fe42a72p-1}, /* 173 */
   {0x1.26d36p-36, 0x1.0485626ba9785p-1, 0x1.b8c38d26da7cp-1}, /* 174 */
   {-0x1.56462p-35, 0x1.05df3ec14d8c2p-1, 0x1.b7f6686f8c2b6p-1}, /* 175 */
   {-0x1.89d9ap-36, 0x1.0738799126acdp-1, 0x1.b728345235eacp-1}, /* 176 */
   {0x1.ddb76p-36, 0x1.089112046c423p-1, 0x1.b658f14f19e41p-1}, /* 177 */
   {-0x1.75d84p-36, 0x1.09e9074081746p-1, 0x1.b5889fe9b9be3p-1}, /* 178 */
   {0x1.a074p-41, 0x1.0b4058790116ap-1, 0x1.b4b7409de23c8p-1}, /* 179 */
   {0x1.6fb47p-35, 0x1.0c9704d7c454cp-1, 0x1.b3e4d3ee2671cp-1}, /* 180 */
   {-0x1.67247p-36, 0x1.0ded0b83cc9bap-1, 0x1.b3115a5fcc743p-1}, /* 181 */
   {0x1.40f08p-38, 0x1.0f426bb2de59ep-1, 0x1.b23cd46fdfd86p-1}, /* 182 */
   {-0x1.3e38cp-38, 0x1.1097248cd5b15p-1, 0x1.b16742a4eb736p-1}, /* 183 */
   {-0x1.4907ap-36, 0x1.11eb3540da5e9p-1, 0x1.b090a581da60fp-1}, /* 184 */
   {0x1.7e9cdp-36, 0x1.133e9cffdfb7fp-1, 0x1.afb8fd8953ef7p-1}, /* 185 */
   {-0x1.29ea8p-37, 0x1.14915af2d45a5p-1, 0x1.aee04b4400794p-1}, /* 186 */
   {-0x1.3f98dap-33, 0x1.15e36e4727973p-1, 0x1.ae068f38a0b31p-1}, /* 187 */
   {-0x1.f08a8p-37, 0x1.1734d63d4a424p-1, 0x1.ad2bc9e287a67p-1}, /* 188 */
   {0x1.f4d568p-35, 0x1.188591f6368cap-1, 0x1.ac4ffbd240a43p-1}, /* 189 */
   {0x1.4044fp-35, 0x1.19d5a0a0cf9bcp-1, 0x1.ab7325905720fp-1}, /* 190 */
   {-0x1.48756p-36, 0x1.1b250170604cfp-1, 0x1.aa9547a35a426p-1}, /* 191 */
   {-0x1.0b381p-36, 0x1.1c73b39a380bbp-1, 0x1.a9b662915eb3ap-1}, /* 192 */
   {-0x1.01fdp-39, 0x1.1dc1b64daf83p-1, 0x1.a8d676e553cffp-1}, /* 193 */
   {-0x1.621c1p-36, 0x1.1f0f08bae2169p-1, 0x1.a7f5852a9a57cp-1}, /* 194 */
   {0x1.7015d8p-34, 0x1.205baa1b1196bp-1, 0x1.a7138de74aca6p-1}, /* 195 */
   {0x1.b4708p-36, 0x1.21a79994595c6p-1, 0x1.a63091af6dc26p-1}, /* 196 */
   {-0x1.b9d398p-34, 0x1.22f2d65e4b185p-1, 0x1.a54c910c24163p-1}, /* 197 */
   {-0x1.d9b08p-38, 0x1.243d5fb93e635p-1, 0x1.a4678c814ec2fp-1}, /* 198 */
   {0x1.4f27d8p-35, 0x1.258734cd666bdp-1, 0x1.a38184a465eacp-1}, /* 199 */
   {-0x1.b2008p-35, 0x1.26d054cba3cefp-1, 0x1.a29a7a05eb03ap-1}, /* 200 */
   {0x1.73eb6p-37, 0x1.2818bef54af27p-1, 0x1.a1b26d2bb5ff1p-1}, /* 201 */
   {-0x1.5edp-39, 0x1.296072760dc14p-1, 0x1.a0c95eabc3946p-1}, /* 202 */
   {0x1.41f19cp-34, 0x1.2aa76e8ae43b8p-1, 0x1.9fdf4f10c6a62p-1}, /* 203 */
   {-0x1.1ae4d8p-35, 0x1.2bedb25e471a5p-1, 0x1.9ef43ef39f48bp-1}, /* 204 */
   {-0x1.27d7p-37, 0x1.2d333d348bc93p-1, 0x1.9e082edb869f4p-1}, /* 205 */
   {0x1.fa6p-38, 0x1.2e780e3ededa4p-1, 0x1.9d1b1f5e6d50ap-1}, /* 206 */
   {-0x1.1484ep-35, 0x1.2fbc24b2e355fp-1, 0x1.9c2d111009089p-1}, /* 207 */
   {0x1.51c2ep-35, 0x1.30ff7fcfc1287p-1, 0x1.9b3e047dfc671p-1}, /* 208 */
   {0x1.ba157p-36, 0x1.32421ec5b0a14p-1, 0x1.9a4dfa41e0baep-1}, /* 209 */
   {0x1.96c6ep-36, 0x1.338400d1c8557p-1, 0x1.995cf2ecc0ef5p-1}, /* 210 */
   {0x1.fd47fp-36, 0x1.34c5252d53ef4p-1, 0x1.986aef136621p-1}, /* 211 */
   {0x1.4f998p-35, 0x1.36058b1209277p-1, 0x1.9777ef4b3e40dp-1}, /* 212 */
   {0x1.4f0fep-37, 0x1.374531b880713p-1, 0x1.9683f42b87ffdp-1}, /* 213 */
   {0x1.cde83p-36, 0x1.3884185f1e0e9p-1, 0x1.958efe48096ep-1}, /* 214 */
   {-0x1.6a026p-34, 0x1.39c23e39e04aep-1, 0x1.94990e3d7d977p-1}, /* 215 */
   {-0x1.64963p-35, 0x1.3affa2904b79p-1, 0x1.93a2249a7edbp-1}, /* 216 */
   {0x1.4c462p-36, 0x1.3c3c4498e98ebp-1, 0x1.92aa41fbb951cp-1}, /* 217 */
   {-0x1.c74bep-36, 0x1.3d78238b3fa7cp-1, 0x1.91b166fe2793ep-1}, /* 218 */
   {-0x1.45fdcp-37, 0x1.3eb33eab7c36dp-1, 0x1.90b79435c59f5p-1}, /* 219 */
   {-0x1.1af724p-34, 0x1.3fed95319f619p-1, 0x1.8fbcca4124bacp-1}, /* 220 */
   {0x1.6b1a9p-36, 0x1.41272664af235p-1, 0x1.8ec109b3d3e33p-1}, /* 221 */
   {-0x1.25404p-38, 0x1.425ff178b9ff7p-1, 0x1.8dc453318dcdep-1}, /* 222 */
   {0x1.9f9ap-40, 0x1.4397f5b2b4074p-1, 0x1.8cc6a75177808p-1}, /* 223 */
   {0x1.acc6ap-35, 0x1.44cf32529a80ap-1, 0x1.8bc806afa62d5p-1}, /* 224 */
   {0x1.7a45788p-31, 0x1.4605a6af55a26p-1, 0x1.8ac871d63c272p-1}, /* 225 */
   {0x1.29a3bp-35, 0x1.473b51baeec8dp-1, 0x1.89c7e9a3b27b3p-1}, /* 226 */
   {-0x1.a889c4p-34, 0x1.4870330209f8cp-1, 0x1.88c66e77d9496p-1}, /* 227 */
   {0x1.32ecap-36, 0x1.49a449ba6906ap-1, 0x1.87c400fb07b89p-1}, /* 228 */
   {-0x1.705e1p-36, 0x1.4ad795159561cp-1, 0x1.86c0a1da650c2p-1}, /* 229 */
   {0x1.40e304p-34, 0x1.4c0a1461bf5dap-1, 0x1.85bc51ac07c91p-1}, /* 230 */
   {-0x1.1e5018p-35, 0x1.4d3bc6d43485cp-1, 0x1.84b7111c1cf5bp-1}, /* 231 */
   {-0x1.7028p-38, 0x1.4e6cabbe07a2p-1, 0x1.83b0e0c028ae3p-1}, /* 232 */
   {0x1.d4a58p-37, 0x1.4f9cc25d55449p-1, 0x1.82a9c13edbbd6p-1}, /* 233 */
   {-0x1.065fe4p-34, 0x1.50cc09f32d34ap-1, 0x1.81a1b33d75e34p-1}, /* 234 */
   {-0x1.8054ap-34, 0x1.51fa81ca0eb38p-1, 0x1.8098b75a0237p-1}, /* 235 */
   {0x1.f320cp-38, 0x1.5328292a7ec46p-1, 0x1.7f8ece35308bcp-1}, /* 236 */
   {-0x1.bf8aap-35, 0x1.5454ff4f4caa3p-1, 0x1.7e83f87cd6b28p-1}, /* 237 */
   {0x1.27e92p-36, 0x1.5581038a223bbp-1, 0x1.7d7836cb98d71p-1}, /* 238 */
   {0x1.4e179cp-34, 0x1.56ac351c8222ep-1, 0x1.7c6b89cb6ebbbp-1}, /* 239 */
   {-0x1.39219p-36, 0x1.57d693481890ap-1, 0x1.7b5df22750237p-1}, /* 240 */
   {-0x1.b78c8p-38, 0x1.59001d5f3278ep-1, 0x1.7a4f707c33a4dp-1}, /* 241 */
   {0x1.8a601p-35, 0x1.5a28d2a79f971p-1, 0x1.79400573528a5p-1}, /* 242 */
   {0x1.11bcb8p-35, 0x1.5b50b2663331dp-1, 0x1.782fb1b7e7877p-1}, /* 243 */
   {0x1.9e7978p-35, 0x1.5c77bbe82d185p-1, 0x1.771e75ee7c0a4p-1}, /* 244 */
   {-0x1.881cep-36, 0x1.5d9dee730248fp-1, 0x1.760c52c3d6c11p-1}, /* 245 */
   {-0x1.a0d318p-35, 0x1.5ec349565a11ap-1, 0x1.74f948dc4db73p-1}, /* 246 */
   {-0x1.17fa5p-36, 0x1.5fe7cbddb6e82p-1, 0x1.73e558e110b71p-1}, /* 247 */
   {0x1.08f7ap-37, 0x1.610b75521e2ap-1, 0x1.72d0837eb839dp-1}, /* 248 */
   {0x1.f9bp-43, 0x1.622e44fec46d9p-1, 0x1.71bac960e1f67p-1}, /* 249 */
   {-0x1.e61p-37, 0x1.63503a3138502p-1, 0x1.70a42b31fb4e1p-1}, /* 250 */
   {-0x1.4f2198p-35, 0x1.647154367b4e7p-1, 0x1.6f8ca99e04332p-1}, /* 251 */
   {0x1.c83e76p-33, 0x1.659192670b46ap-1, 0x1.6e744546d8881p-1}, /* 252 */
   {-0x1.91638p-39, 0x1.66b0f3f50f19cp-1, 0x1.6d5afef4c6982p-1}, /* 253 */
   {-0x1.2f1d8p-35, 0x1.67cf7847c834ap-1, 0x1.6c40d73d66c2p-1}, /* 254 */
   {0x1.21862p-36, 0x1.68ed1eaabb0f7p-1, 0x1.6b25ced25dddfp-1}, /* 255 */
};

// Multiply exactly a and b, such that *hi + *lo = a * b. 
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma (a, b, -*hi);
}

/* Multiply a double with a double double : a * (bh + bl)
   with error bounded by ulp(lo) */
static inline void s_mul (double *hi, double *lo, double a, double bh,
                          double bl) {
  a_mul (hi, lo, a, bh); /* exact */
  *lo = __builtin_fma (a, bl, *lo);
  /* the error is bounded by ulp(lo), where |lo| < |a*bl| + ulp(hi) */
}

// Returns (ah + al) * (bh + bl) - (al * bl)
// We can ignore al * bl when assuming al <= ulp(ah) and bl <= ulp(bh)
static inline void d_mul(double *hi, double *lo, double ah, double al,
                         double bh, double bl) {
  double s, t;

  a_mul(hi, &s, ah, bh);
  t = __builtin_fma(al, bh, s);
  *lo = __builtin_fma(ah, bl, t);
}

static inline void
fast_two_sum(double *hi, double *lo, double a, double b)
{
  double e;

  //  assert (a == 0 || __builtin_fabs (a) >= __builtin_fabs (b));
  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
}

/* Put in h+l an approximation of sin2pi(xh+xl), for |xh| < 2^-11 + 2^-30,
   and |xl| < 2^-52.36, with absolute error < 2^-77.09. */
static void
evalPSfast (double *h, double *l, double xh, double xl)
{
  double uh, ul, t;
  a_mul (&uh, &ul, xh, xh);
  ul = __builtin_fma (xh + xh, xl, ul);
  // uh+ul approximates (xh+xl)^2
  *h = PSfast[4]; // degree 7
  *h = __builtin_fma (*h, uh, PSfast[3]); // degree 5
  *h = __builtin_fma (*h, uh, PSfast[2]); // degree 3
  s_mul (h, l, *h, uh, ul);
  fast_two_sum (h, &t, PSfast[0], *h);
  *l += PSfast[1] + t;
  // multiply by xh+xl
  d_mul (h, l, *h, *l, xh, xl);
}

/* Put in h+l an approximation of cos2pi(xh+xl), for |xh| < 2^-11 + 2^-30,
   and |xl| < 2^-52.36, with relative error < 2^-69.96. */
static void
evalPCfast (double *h, double *l, double xh, double xl)
{
  double uh, ul, t;
  a_mul (&uh, &ul, xh, xh);
  ul = __builtin_fma (xh + xh, xl, ul);
  // uh+ul approximates (xh+xl)^2
  *h = PCfast[4]; // degree 6
  *h = __builtin_fma (*h, uh, PCfast[3]); // degree 4
  *h = __builtin_fma (*h, uh, PCfast[2]); // degree 2
  s_mul (h, l, *h, uh, ul);
  fast_two_sum (h, &t, PCfast[0], *h);
  *l += PCfast[1] + t;
}

/* Put in Y an approximation of sin2pi(X), for 0 <= X < 2^-11,
   where X2 approximates X^2.
   Absolute error bounded by 2^-132.999 with 0 <= Y < 0.003068
   (see evalPS() in sin.sage), and relative error bounded by
   2^-123.651 (see evalPSrel() in sin.sage). */
static void
evalPS (dint64_t *Y, dint64_t *X, dint64_t *X2)
{
  mul_dint_21 (Y, X2, PS+5); // degree 11
  add_dint (Y, Y, PS+4);     // degree 9
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PS+3);     // degree 7
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PS+2);     // degree 5
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PS+1);     // degree 3
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PS+0);     // degree 1
  mul_dint (Y, Y, X);        // multiply by X
}

/* Put in Y an approximation of cos2pi(X), for 0 <= X < 2^-11,
   where X2 approximates X^2.
   Absolute/relative error bounded by 2^-125.999 with 0.999995 < Y <= 1
   (see evalPC() in sin.sage). */
static void
evalPC (dint64_t *Y, dint64_t *X2)
{
  mul_dint_21 (Y, X2, PC+5); // degree 10
  add_dint (Y, Y, PC+4);     // degree 8
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PC+3);     // degree 6
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PC+2);     // degree 4
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PC+1);     // degree 2
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PC+0);     // degree 0
}

// normalize X such that X->hi has its most significant bit set (if X <> 0)
static void
normalize (dint64_t *X)
{
  int cnt;
  if (X->hi != 0)
  {
    cnt = __builtin_clzl (X->hi);
    if (cnt)
    {
      X->hi = (X->hi << cnt) | (X->lo >> (64 - cnt));
      X->lo = X->lo << cnt;
    }
    X->ex -= cnt;
  }
  else if (X->lo != 0)
  {
    cnt = __builtin_clzl (X->lo);
    X->hi = X->lo << cnt;
    X->lo = 0;
    X->ex -= 64 + cnt;
  }
}

/* Approximate X/(2pi) mod 1. If Xin is the input value, and Xout the
   output value, we have:
   |Xout - (Xin/(2pi) mod 1)| < 2^-126.67*|Xout|
   Assert X is normalized at input, and normalize X at output.
*/
static void
reduce (dint64_t *X)
{
  int e = X->ex;
  u128 u;

  if (e <= 1) // |X| < 2
  {
    /* multiply by T[0]/2^64 + T[1]/2^128, where
       |T[0]/2^64 + T[1]/2^128 - 1/(2pi)| < 2^-130.22 */
    u = (u128) X->hi * (u128) T[1];
    uint64_t tiny = u;
    X->lo = u >> 64;
    u = (u128) X->hi * (u128) T[0];
    X->lo += u;
    X->hi = (u >> 64) + (X->lo < (uint64_t) u);
    /* hi + lo/2^64 + tiny/2^128 = hi_in * (T[0]/2^64 + T[1]/2^128) thus
       |hi + lo/2^64 + tiny/2^128 - hi_in/(2*pi)| < hi_in * 2^-130.22
       Since X is normalized at input, hi_in >= 2^63, and since T[0] >= 2^61,
       we have hi >= 2^(63+61-64) = 2^60, thus the normalize() below
       perform a left shift by at most 3 bits */
    int e = X->ex;
    normalize (X);
    e = e - X->ex;
    assert (0 <= e && e <= 3);
    // put the upper e bits of tiny into X->lo
    if (e)
      X->lo |= tiny >> (64 - e);
    /* The error is bounded by 2^-130.22 (relative) + ulp(lo) (absolute).
       Since now X->hi >= 2^63, the absolute error of ulp(lo) converts into
       a relative error of less than 2^-127.
       This yields a maximal relative error of:
       (1 + 2^-130.22) * (1 + 2^-127) - 1 < 2^-126.852.
    */
    return;
  }

  // now 2 <= e <= 1024
  // assert (2 <= e && e <= 1024);

  /* The upper 64-bit word X->hi corresponds to hi/2^64*2^e, if multiplied by
     T[i]/2^((i+1)*64) it yields hi*T[i]/2^128 * 2^(e-i*64).
     If e-64i <= -128, it contributes to less than 2^-128;
     if e-64i >= 128, it yields an integer, which is 0 modulo 1.
     We thus only consider the values of i such that -127 <= e-64i <= 127,
     i.e., (-127+e)/64 <= i <= (127+e)/64.
     Up to 4 consecutive values of T[i] can contribute (only 3 when e is a
     multiple of 64). */
  int i = (e < 127) ? 0 : (e - 127 + 64 - 1) / 64; // ceil((e-127)/64)
  // 0 <= i <= 15
  uint64_t c[5];
  u = (u128) X->hi * (u128) T[i+3]; // i+3 <= 18
  c[0] = u;
  c[1] = u >> 64;
  u = (u128) X->hi * (u128) T[i+2];
  c[1] += u;
  c[2] = (u >> 64) + (c[1] < (uint64_t) u);
  u = (u128) X->hi * (u128) T[i+1];
  c[2] += u;
  c[3] = (u >> 64) + (c[2] < (uint64_t) u);
  u = (u128) X->hi * (u128) T[i];
  c[3] += u;
  c[4] = (u >> 64) + (c[3] < (uint64_t) u);

  /* up to here, the ignored part hi*(T[i+4]+T[i+5]+...) can contribute by
     less than 2^64 in c[0], thus less than 1 in c[1] */

  int f = e - 64 * i; // hi*T[i]/2^128 is multiplied by 2^f
  /* {c, 5} = hi*(T[i]+T[i+1]/2^64+T[i+2]/2^128+T[i+3]/2^192) */
  // assert (2 <= f && f <= 127);
  /* now shift c[0..4] by f bits to the left */
  uint64_t tiny;
  if (f < 64)
  {
    X->hi = (c[4] << f) | (c[3] >> (64 - f));
    X->lo = (c[3] << f) | (c[2] >> (64 - f));
    tiny = (c[2] << f) | (c[1] >> (64 - f));
    /* the ignored part was less than 1 in c[1],
       thus less than 2^(f-64) <= 1/2 in tiny */
  }
  else if (f == 64)
  {
    X->hi = c[3];
    X->lo = c[2];
    tiny = c[1];
    /* the ignored part was less than 1 in c[1],
       thus less than 1 in tiny */
  }
  else /* 65 <= f <= 127: this case can only occur when e >= 65 */
  {
    int g = f - 64; /* 1 <= g <= 63 */
    /* we compute an extra term */
    u = (u128) X->hi * (u128) T[i+4]; // i+4 <= 19
    u = u >> 64;
    c[0] += u;
    c[1] += (c[0] < u);
    c[2] += (c[0] < u) && c[1] == 0;
    c[3] += (c[0] < u) && c[1] == 0 && c[2] == 0;
    c[4] += (c[0] < u) && c[1] == 0 && c[2] == 0 && c[3] == 0;
    X->hi = (c[3] << g) | (c[2] >> (64 - g));
    X->lo = (c[2] << g) | (c[1] >> (64 - g));
    tiny = (c[1] << g) | (c[0] >> (64 - g));
    /* the ignored part was less than 1 in c[0],
       thus less than 1/2 in tiny */
  }
  /* The approximation error between X/in(2pi) mod 1 and
     X->hi/2^64 + X->lo/2^128 + tiny/2^192 is:
     (a) the ignored part in tiny, which is less than ulp(tiny),
         thus less than 1/2^192;
     (b) the ignored terms hi*T[i+4] + ... or hi*T[i+5] + ...,
         which accumulate to less than ulp(tiny) too, thus
         less than 1/2^192.
     Thus the approximation error is less than 2^-191 (absolute).
  */
  X->ex = 0;
  normalize (X);
  /* the worst case (for 2^25 <= x < 2^1024) is X->ex = -61, attained
     for |x| = 0x1.6ac5b262ca1ffp+851 */
  if (X->ex < 0) // put the upper -ex bits of tiny into low bits of lo
    X->lo |= tiny >> (64 + X->ex);
  /* Since X->ex >= -61, it means X >= 2^-62 before the normalization,
     thus the maximal absolute error of 2^-191 yields a relative error
     bounded by 2^-191/2^-62 = 2^-129.
     There is an additional truncation error (for tiny) of at most 1 ulp
     of X->lo, thus at most 2^-127.
     The relative error is thus bounded by 2^-126.67. */
}

/* Given Xin:=X with 0 <= Xin < 1, return i and modify X such that
   Xin = i/2^11 + Xout, with 0 <= Xout < 2^-11.
   This operation is exact. */
static int
reduce2 (dint64_t *X)
{
  // assert (X->ex <= 0);
  if (X->ex <= -11)
    return 0;
  int sh = 64 - 11 - X->ex;
  int i = X->hi >> sh;
  X->hi = X->hi & ((1ul << sh) - 1);
  normalize (X);
  return i;
}

/* h+l <- c1/2^64 + c0/2^128 */
static void
set_dd (double *h, double *l, uint64_t c1, uint64_t c0)
{
  uint64_t e, f, g;
  b64u64_u t;
  if (c1)
    {
      e = __builtin_clzl (c1);
      if (e)
        {
          c1 = (c1 << e) | (c0 >> (64 - e));
          c0 = c0 << e;
        }
      f = 0x3fe - e;
      t.u = (f << 52) | ((c1 << 1) >> 12);
      *h = t.f;
      c0 = (c1 << 53) | (c0 >> 11);
      if (c0)
        {
          g = __builtin_clzl (c0);
          if (g)
            c0 = c0 << g;
          t.u = ((f - 53 - g) << 52) | ((c0 << 1) >> 12);
          *l = t.f;
        }
      else
        *l = 0;
    }
  else if (c0)
    {
      e = __builtin_clzl (c0);
      f = 0x3fe - 64 - e;
      c0 = c0 << (e+1); // most significant bit shifted out
      /* put the upper 52 bits of c0 into h */
      t.u = (f << 52) | (c0 >> 12);
      *h = t.f;
      /* put the lower 12 bits of c0 into l */
      c0 = c0 << 52;
      if (c0)
        {
          int g = __builtin_clzl (c0);
          c0 = c0 << (g+1);
          t.u = ((f - 64 - g) << 52) | (c0 >> 12);
          *l = t.f;
        }
      else
        *l = 0;
    }
  else
    *h = *l = 0;
  /* Since we truncate from two 64-bit words to a double-double,
     we have another truncation error of less than 2^-106, thus
     the absolute error is bounded as follows:
     | h + l - frac(x/(2pi)) | < 2^-75.999 + 2^-106 < 2^-75.998 */
}

/* Assuming 0x1.7137449123ef6p-26 < x < +Inf,
   return i and set h,l such that i/2^11+h+l approximates frac(x/(2pi)).
   If x <= 0x1.921fb54442d18p+2:
   | i/2^11 + h + l - frac(x/(2pi)) | < 2^-104.116 * |i/2^11 + h + l|
   with |h| < 2^-11 and |l| < 2^-52.36.

   Otherwise only the absolute error is bounded:
   | i/2^11 + h + l - frac(x/(2pi)) | < 2^-75.998
   with 0 <= h < 2^-11 and |l| < 2^-53.

   In both cases we have |l| < 2^-51.64*|i/2^11 + h|.
*/
static int
reduce_fast (double *h, double *l, double x)
{
  if (x <= 0x1.921fb54442d17p+2) // x < 2*pi
    {
      /* | CH+CL - 1/(2pi) | < 2^-110.523 */
#define CH 0x1.45f306dc9c883p-3
#define CL -0x1.6b01ec5417056p-57
      a_mul (h, l, CH, x);            // exact
      *l = __builtin_fma (CL, x, *l);
      /* The error in the above fma() is at most ulp(l),
         where |l| <= CL*|x|+|l_in|.
         Assume 2^(e-1) <= x < 2^e.
         Then |h| < 2^(e-2) and |l_in| <= 1/2 ulp(2^(e-2)) = 2^(e-55),
         where l_in is the value of l after a_mul.
         Then |l| <= CL*x + 2^(e-55) <= 2^e*(CL+2-55) < 2^e * 2^-55.6.
         The rounding error of the fma() is bounded by
         ulp(l) <= 2^e * ulp(2^-55.6) = 2^(e-108).
         The error due to the approximation of 1/(2pi)
         is bounded by 2^-110.523*x <= 2^(e-110.523).
         Adding both errors yields:
         |h + l - x/(2pi)| < 2^e * (2^-108 + 2^-110.523) < 2^e * 2^-107.768.
         Since |x/(2pi)| > 2^(e-1)/(2pi), the relative error is bounded by:
         2^e * 2^-107.768 / (2^(e-1)/(2pi)) = 4pi * 2^-107.768 < 2^-104.116.

         Bound on l: since |h| < 1, we have after |l| <= ulp(h) <= 2^-53
         after a_mul(), and then |l| <= |CL|*0x1.921fb54442d17p+2 + 2^-53
         < 2^-52.36.

         Bound on l relative to h: after a_mul() we have |l| <= ulp(h)
         <= 2^-52*h. After fma() we have |l| <= CL*x + 2^-52*h
         <= 2^-53.84*CH*x + 2^-52*h <= (2^-53.84+2^-52)*h < 2^-51.64*h.
      */
    }
  else // x > 0x1.921fb54442d17p+2
    {
      b64u64_u t = {.f = x};
      int e = (t.u >> 52) & 0x7ff; /* 1025 <= e <= 2046 */
      if (x == TRACE) printf ("e=%d\n", e);
      /* We have 2^(e-1023) <= x < 2^(e-1022), thus
         ulp(x) is a multiple of 2^(e-1075), for example
         if x is just above 2*pi, e=1025, 2^2 <= x < 2^e,
         and ulp(x) is a multiple of 2^-50.
         On the other side 1/(2pi) ~ T[0]/2^64 + T[1]/2^128 + T[2]/2^192 + ...
         Let i be the smallest integer such that 2^(e-1075)/2^(64*(i+1))
         is not an integer, i.e., e - 1139 - 64i < 0, i.e.,
         i >= (e-1138)/64. */
      uint64_t m = (1ul << 52) | (t.u & 0xffffffffffffful);
      uint64_t c[3];
      u128 u;
      // x = m/2^53 * 2^(e-1022)
      if (e <= 1074) // 1025 <= e <= 1074: 2^2 <= x < 2^52
        {
          /* In that case the contribution of x*T[2]/2^192 is less than
             2^(52+64-192) <= 2^-76. */
          u = (u128) m * (u128) T[1];
          c[0] = u;
          c[1] = u >> 64;
          u = (u128) m * (u128) T[0];
          c[1] += u;
          c[2] = (u >> 64) + (c[1] < (uint64_t) u);
          /* | c[2]*2^128+c[1]*2^64+c[0] - m/(2pi)*2^128 | < m*T[2]/2^64 < 2^53
             thus:
             | (c[2]*2^128+c[1]*2^64+c[0])*2^(e-1203) - x/(2pi) | < 2^(e-1150)
             The low 1075-e bits of c[2] contribute to frac(x/(2pi)).
          */
          e = 1075 - e; // 1 <= e <= 50
          // e is the number of low bits of C[2] contributing to frac(x/(2pi))
          assert (1 <= e && e <= 50);
        }
      else // 1075 <= e <= 2046, 2^52 <= x < 2^1024
        {
          int i = (e - 1138 + 63) / 64; // i = ceil((e-1138)/64), 0 <= i <= 15
          /* m*T[i] contributes to f = 1139 + 64*i - e bits to frac(x/(2pi))
             with 1 <= f <= 64
             m*T[i+1] contributes a multiple of 2^(-f-64),
                      and at most to 2^(53-f)
             m*T[i+2] contributes a multiple of 2^(-f-128),
                      and at most to 2^(-11-f)
             m*T[i+3] contributes a multiple of 2^(-f-192),
                      and at most to 2^(-75-f) <= 2^-76
          */
          u = (u128) m * (u128) T[i+2];
          c[0] = u;
          c[1] = u >> 64;
          u = (u128) m * (u128) T[i+1];
          c[1] += u;
          c[2] = (u >> 64) + (c[1] < (uint64_t) u);
          u = (u128) m * (u128) T[i];
          c[2] += u;
          e = 1139 + (i<<6) - e; // 1 <= e <= 64
          // e is the number of low bits of C[2] contributing to frac(x/(2pi))
          assert (1 <= e && e <= 64);
        }
      if (e == 64)
        {
          c[0] = c[1];
          c[1] = c[2];
        }
      else
        {
          c[0] = (c[1] << (64 - e)) | c[0] >> e;
          c[1] = (c[2] << (64 - e)) | c[1] >> e;
        }
      if (x == TRACE) printf ("c[1]=%lu c[0]=%lu\n", c[1], c[0]);
      /* In all cases the ignored contribution from x*T[2] or x*T[i+3]
         is less than 2^-76,
         and the truncated part from the above shift is less than 2^-128 thus:
         | c[1]/2^64 + c[0]/2^128 - frac(x/(2pi)) | < 2^-76+2^-128 < 2^-75.999
      */
      set_dd (h, l, c[1], c[0]);
      /* set_dd() ensures |h| < 1 and |l| < ulp(h) <= 2^-53 */
    }
  if (x == TRACE) printf ("h=%la l=%la\n", *h, *l);

  double i = __builtin_floor (*h * 0x1p11);
  *h = __builtin_fma (i, -0x1p-11, *h);
  return i;
}

/* return the maximal relative error */
static double
sin_fast (double *h, double *l, double x)
{
  int neg = x < 0, is_sin = 1;
  double absx = neg ? -x : x;
  if (x == TRACE) printf ("absx=%la\n", absx);

  /* now x > 0x1.7137449123ef6p-26 */
  int i = reduce_fast (h, l, absx);
  assert (0 <= i && i < 2048);
  if (x == TRACE) printf ("i=%d h=%la l=%la\n", i, *h, *l);
  /* If x <= 0x1.921fb54442d18p+2:
     | i/2^11 + h + l - frac(x/(2pi)) | < 2^-104.116 * |h + l|
     otherwise | i/2^11 + h + l - frac(x/(2pi)) | < 2^-75.998. */

  // if i >= 2^10: 1/2 <= frac(x/(2pi)) < 1 thus pi <= x <= 2pi
  neg = neg ^ (i >> 10);
  i = i & 0x3ff;

  // now i < 2^10
  // if i >= 2^9: 1/4 <= frac(x/(2pi)) < 1/2 thus pi/2 <= x <= pi
  is_sin = is_sin ^ (i >> 9);
  i = i & 0x1ff;

  // now 0 <= i < 2^9
  // if i >= 2^8: 1/8 <= frac(x/(2pi)) < 1/4 thus pi/4 <= x <= pi/2
  if (i & 0x100) // case 
    {
      is_sin = !is_sin;
      i = 0x1ff - i;
      /* 0x1p-11 - h is exact below: indeed, reduce_fast first computes
         a first value of h (say h0, with 0 <= h0 < 1), then i = floor(h0*2^11)
         and h1 = h0 - 2^11*i with 0 <= h1 < 2^-11.
         If i >= 2^8 here, this implies h0 >= 1/2^3, thus ulp(h0) >= 2^-55:
         h0 and h1 are integer multiples of 2^-55.
         Thus h1 = k*2^-55 with 0 <= k < 2^44 (since 0 <= h1 < 2^-11).
         Then 0x1p-11 - h = (2^44-k)*2^-55 is exactly representable. */
      *h = 0x1p-11 - *h;
      *l = -*l;
    }

  if (x == TRACE) printf ("neg=%d is_sin=%d i=%d h=%la l=%la\n",
                          neg, is_sin, i, *h, *l);
  /* Now 0 <= i < 256 and 0 <= h+l < 2^-11
     with | i/2^11 + h + l - frac(x/(2pi)) | cmod 1/4 < eps
     with |eps| < 2^-104.116 * |h + l| for i=0
          |eps| < 2^-104.116 for i<>0 and x <= 0x1.921fb54442d18p+2
          |eps| < 2^-75.998 otherwise
     If is_sin=1, sin |x| = sin2pi (R + eps);
     if is_sin=0, sin |x| = cos2pi (R + eps).
     In both cases R = i/2^11 + h+l, 0 <= R < 1/4, and |eps| as above.
  */
  double sh, sl, ch, cl;
  /* since the SC[] table evaluates at i/2^11 + eps and not at i/2^11,
     we must subtract eps from h+l */
  /* Here h = k*2^-55 with 0 <= k < 2^44, and SC[i][0] is an integer
     multiple of 2^-63, with |SC[i][0]| < 2^-30, thus SC[i][0] = m*2^-63
     with |m| < 2^33. It follows h-SC[i][0] = (k*2^8 + m)*2^-63 with
     2^52 - 2^33 < k*2^8 + m < 2^52 + 2^33, thus h-SC[i][0] is exact.
     Now |h| < 2^-11 + 2^-30. */
  *h -= SC[i][0];
  // from reduce_fast() we have |l| < 2^-52.36
  static double maxl = 0;
  if (__builtin_fabs (*l) > maxl) printf ("maxl=%la\n", maxl = __builtin_fabs (*l));
  if (x == TRACE) printf ("modified h=%la\n", *h);
  evalPSfast (&sh, &sl, *h, *l); // absolute error < 2^-77.09
  if (x == TRACE) printf ("sh=%la sl=%la\n", sh, sl);
  evalPCfast (&ch, &cl, *h, *l); // relative error < 2^-69.96
  if (x == TRACE) printf ("ch=%la cl=%la\n", ch, cl);
  double err;
  if (is_sin)
    {
      s_mul (&sh, &sl, SC[i][2], sh, sl);
      s_mul (&ch, &cl, SC[i][1], ch, cl);
      fast_two_sum (h, l, ch, sh);
      *l += sl + cl;
      /* relative error bounded by 2^-63.48 < 0x1.70p-64 */
      err = 0x1.70p-64;
    }
  else
    {
      /* Let xi = i/2^11 + SC[i][0], then we have modulo 1/4:
         x/(2pi) = i/2^11 + h + l (original value of h)
                 = xi + h + l     (new value of h)
      */
      s_mul (&ch, &cl, SC[i][2], ch, cl);
      s_mul (&sh, &sl, SC[i][1], sh, sl);
      fast_two_sum (h, l, ch, -sh);
      *l += cl - sl;
      /* relative error bounded by 2^-64.61 < 0x1.50p-65 */
      err = 0x1.50p-65;
    }
  if (x == TRACE) printf ("h=%la l=%la neg=%d\n", *h, *l, neg);
  static double sgn[2] = {1.0, -1.0};
  *h *= sgn[neg];
  *l *= sgn[neg];
  if (x == TRACE) printf ("h=%la l=%la\n", *h, *l);
  return err;
}

static double
sin_accurate (double x)
{
  // printf ("enter sin_accurate, x=%la\n", x);
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
    // Taylor expansion of sin(x) is x - x^3/6 around zero
    return __builtin_fma (x, -0x1p-54, x);

  double absx = (x > 0) ? x : -x;

  /* now x > 0x1.7137449123ef6p-26 */
  dint64_t X[1];
  dint_fromd (X, absx);

  /* reduce argument */
  reduce (X);
  
  // now |X - x/(2pi) mod 1| < 2^-126.67*X, with 0 <= X < 1.

  int neg = x < 0, is_sin = 1;

  // Write X = i/2^11 + r with 0 <= r < 2^11.
  int i = reduce2 (X); // exact
  // if (x == TRACE) { printf ("i=%d Xred=", i); print_dint (X); }

  if (i & 0x400) // pi <= x < 2*pi: sin(x) = -sin(x-pi)
  {
    neg = !neg;
    i = i & 0x3ff;
  }

  // now i < 2^10

  if (i & 0x200) // pi/2 <= x < pi: sin(x) = cos(x-pi/2)
  {
    is_sin = 0;
    i = i & 0x1ff;
  }

  // now 0 <= i < 2^9

  if (i & 0x100)
    // pi/4 <= x < pi/2: sin(x) = cos(pi/2-x), cos(x) = sin(pi/2-x)
  {
    is_sin = !is_sin;
    // assert (X->sgn == 0);
    X->sgn = 1; // negate X
    // assert (X->ex < -1);
    add_dint (X, &MAGIC, X); // X -> 2^-11 - X
    // here: 256 <= i <= 511
    i = 0x1ff - i;
    // now 0 <= i < 256
  }
  // if (x == TRACE) { printf ("i=%d Xred=", i); print_dint (X); }

  // now 0 <= i < 256 and 0 <= X < 2^-11
  // assert (0 <= i && i < 256);

  /* If is_sin=1, sin |x| = sin2pi (R * (1 + eps))
        (cases 0 <= x < pi/4 and 3pi/4 <= x < pi)
     if is_sin=0, sin |x| = cos2pi (R * (1 + eps))
        (case pi/4 <= x < 3pi/4)
     In both cases R = i/2^11 + X, 0 <= R < 1/4, and |eps| < 2^-126.67.
  */

  dint64_t U[1], V[1], X2[1];
  mul_dint (X2, X, X);       // X2 approximates X^2
  evalPC (U, X2);    // cos2pi(X)
  /* since 0 <= X < 2^-11, we have 0.999 < U <= 1 */
  evalPS (V, X, X2); // sin2pi(X)
  /* since 0 <= X < 2^-11, we have 0 <= V < 0.0005 */
  if (is_sin)
  {
    // sin2pi(x) ~ sin2pi(i/2^11)*cos2pi(X)+cos2pi(i/2^11)*sin2pi(X)
    mul_dint (U, S+i, U);
    /* since 0 <= S[i] < 0.705 and 0.999 < Uin <= 1, we have
       0 <= U < 0.705 */
    mul_dint (V, C+i, V);
    /* For the error analysis, we distinguish the case i=0.
       For i=0, we have S[i]=0 and C[1]=1, thus V is the value computed
       by evalPS() above, with relative error < 2^-123.651.

       For 1 <= i < 256, analyze_sin_case1(rel=true) from sin.sage gives a
       relative error bound of -122.797 (obtained for i=1).
       In all cases, the relative error for the computation of
       sin2pi(i/2^11)*cos2pi(X)+cos2pi(i/2^11)*sin2pi(X) is bounded by -122.797
       not taking into account the approximation error in R:
       |U - sin2pi(R)| < |U| * 2^-122.797, with U the value computed
       after add_dint (U, U, V) below.

       For the approximation error in R, we have:
       sin |x| = sin2pi (R * (1 + eps))
       R = i/2^11 + X, 0 <= R < 1/4, and |eps| < 2^-126.67.
       Thus sin|x| = sin2pi(R+R*eps)
                   = sin2pi(R)+R*eps*2*pi*cos2pi(theta), theta in [R,R+R*eps]
       Since 2*pi*R/sin(2*pi*R) < pi/2 for R < 1/4, it follows:
       | sin|x| - sin2pi(R) | < pi/2*R*|sin(2*pi*R)|
       | sin|x| - sin2pi(R) | < 2^-126.018 * |sin2pi(R)|.

       Adding both errors we get:
       | sin|x| - U | < |U| * 2^-122.797 + 2^-126.018 * |sin2pi(R)|
                      < |U| * 2^-122.797 + 2^-126.018 * |U| * (1 + 2^-122.797)
                      < |U| * 2^-122.650.
    */
  }
  else
  {
    // if (x == TRACE) printf ("is_cos\n");
    // approximate sin2pi(x) by cos2pi(i/2^11)*cos2pi(X)-sin2pi(i/2^11)*sin2pi(X)
    mul_dint (U, C+i, U);
    // if (x == TRACE) { printf ("cos2pi(i/2^11)*cos2pi(X2)="); print_dint (U); }
    mul_dint (V, S+i, V);
    // if (x == TRACE) { printf ("sin2pi(i/2^11)*sin2pi(X2)="); print_dint (V); }
    V->sgn = 1 - V->sgn; // negate V
    /* For 0 <= i < 256, analyze_sin_case2(rel=true) from sin.sage gives a
       relative error bound of -123.540 (obtained for i=0):
       |U - cos2pi(R)| < |U| * 2^-123.540, with U the value computed
       after add_dint (U, U, V) below.

       For the approximation error in R, we have:
       cos |x| = cos2pi (R * (1 + eps))
       R = i/2^11 + X, 0 <= R < 1/4, and |eps| < 2^-126.67.
       Thus sin|x| = cos2pi(R+R*eps)
                   = cos2pi(R)-R*eps*2*pi*sin2pi(theta), theta in [R,R+R*eps]
       Since we have R < 1/4, we have cos2pi(R) >= sqrt(2)/2,
       and it follows:
       | sin|x|/cos2pi(R) - 1 | < 2*pi*R*eps/(sqrt(2)/2)
                                < pi/2*eps/sqrt(2)          [since R < 1/4]
                                < 2^-126.518.
       Adding both errors we get:
       | sin|x| - U | < |U| * 2^-123.540 + 2^-126.518 * |cos2pi(R)|
                      < |U| * 2^-123.540 + 2^-126.518 * |U| * (1 + 2^-123.540)
                      < |U| * 2^-123.367.
    */
  }
  add_dint (U, U, V);
  /* If is_sin=1:
     | sin|x| - U | < |U| * 2^-122.650
     If is_sin=0:
     | cos|x| - U | < |U| * 2^-123.367.
     In all cases the total error is bounded by |U| * 2^-122.650.
     The term |U| * 2^-122.650 contributes to at most 2^(128-122.650) < 41 ulps
     relatively to U->lo.
  */
  uint64_t err = 41;
  uint64_t hi0, hi1, lo0, lo1;
  lo0 = U->lo - err;
  hi0 = U->hi - (lo0 > U->lo);
  lo1 = U->lo + err;
  hi1 = U->hi + (lo1 < U->lo);
  /* check the upper 54 bits are equal */
  if ((hi0 >> 10) != (hi1 >> 10))
    {
      static const double exceptions[][3] = {
        {0x1.e0000000001c2p-20, 0x1.dfffffffff02ep-20, 0x1.dcba692492527p-146},
      };
      for (int i = 0; i < 1; i++)
        {
          if (__builtin_fabs (x) == exceptions[i][0])
            return (x > 0) ? exceptions[i][1] + exceptions[i][2]
              : -exceptions[i][1] - exceptions[i][2];
        }
      printf ("Rounding test of accurate path failed for sin(x)=%la\n", x);
      printf ("Please report the above to core-math@inria.fr\n");
      exit (1);
    }

  // if (x == TRACE) { printf ("U="); print_dint (U); }

  if (neg)
    U->sgn = 1 - U->sgn;

  double y = dint_tod (U);

  return y;
}

double
cr_sin (double x)
{
  b64u64_u t = {.f = x};
  int e = (t.u >> 52) & 0x7ff;

  if (e == 0x7ff) /* NaN, +Inf and -Inf. */
    {
      t.u = ~0ul;
      return t.f;
    }

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
    // Taylor expansion of sin(x) is x - x^3/6 around zero
    // for x=-0, fma (x, -0x1p-54, x) returns +0
    return (x == 0) ? x :__builtin_fma (x, -0x1p-54, x);

  double h, l, err;
  err = sin_fast (&h, &l, x);
  double left  = h + __builtin_fma (h, -err, l);
  double right = h + __builtin_fma (h,  err, l);
  if (left == right)
    {
      if (x == TRACE) printf ("fast path succeeded\n");
      return left;
    }
  if (x == TRACE) printf ("fast path failed\n");

  return sin_accurate (x);
}
