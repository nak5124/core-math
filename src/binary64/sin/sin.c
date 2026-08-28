/* Correctly-rounded sine function for binary64 value.

Copyright (c) 2022-2025 Paul Zimmermann and Tom Hubrecht

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
#include <inttypes.h>
#include <fenv.h> // for fegetround, FE_TONEAREST, FE_DOWNWARD, FE_UPWARD
#ifdef CORE_MATH_SUPPORT_ERRNO
#include <errno.h>
#endif

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

/******************** code copied from dint.h and pow.[ch] *******************/

#if (defined(__clang__) && __clang_major__ >= 14) || (defined(__GNUC__) && __GNUC__ >= 14 && __BITINT_MAXWIDTH__ && __BITINT_MAXWIDTH__ >= 128)
typedef unsigned _BitInt(128) u128;
#else
typedef unsigned __int128 u128;
#endif

/* The dint64_t structure represents a 128-bit number:
   (-1)^sgn*(hi/2^64+lo/2^128)*2^ex */
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
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
#else
typedef union {
  struct {
    u128 r;
    int64_t _ex;
    uint64_t _sgn;
  };
  struct {
    uint64_t hi;
    uint64_t lo;
    int64_t ex;
    uint64_t sgn;
  };
} dint64_t;
#endif


#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
typedef union {
  u128 r;
  struct {
    uint64_t l;
    uint64_t h;
  };
} uint128_t;
#else
typedef union {
  u128 r;
  struct {
    uint64_t h;
    uint64_t l;
  };
} uint128_t;
#endif

typedef union {
  double f;
  uint64_t u;
} f64_u;

// Extract both the mantissa and exponent of a double
static inline void fast_extract (int64_t *e, uint64_t *m, double x) {
  f64_u _x = {.f = x};

  *e = (_x.u >> 52) & 0x7ff;
  *m = (_x.u & (~0ull >> 12)) + (*e ? (1ull << 52) : 0);
  *e = *e - 0x3fe;
}

// Return non-zero if a = 0
static inline int
dint_zero_p (const dint64_t *a)
{
  return a->hi == 0;
}

#if 0
// Prints a dint64_t value for debugging purposes
static inline void print_dint(const dint64_t *a) {
  printf("{.hi=0x%"PRIx64", .lo=0x%"PRIx64", .ex=%"PRId64", .sgn=0x%"PRIx64"}\n", a->hi, a->lo, a->ex,
         a->sgn);
}
#endif

static inline int cmp(int64_t a, int64_t b) { return (a > b) - (a < b); }

static inline int cmpu128 (u128 a, u128 b) { return (a > b) - (a < b); }

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
       where it is said that k is interpreted modulo n. */
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
       thus __builtin_clzll(C) is well-defined below. */
    uint64_t ex = ch ? __builtin_clzll(ch) : 64 + __builtin_clzll(C);
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
      ex = __builtin_clzll (C >> 64);
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

  uint32_t t = __builtin_clzll (a->hi);

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
  uint64_t lo = (a->hi & (~0ull >> ex)) || a->lo;

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
    a->hi = (1ull << 63);
  }
}

// Convert a dint64_t value to a double
static inline double dint_tod(dint64_t *a) {
  subnormalize_dint (a);

  f64_u r = {.u = (a->hi >> 11) | (0x3ffll << 52)};

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

/* This table approximates 1/(2pi) downwards with precision 1280:
   1/(2*pi) ~ T[0]/2^0 + T[1]/2^64 + ... + T[i]/2^(i*64) + ...
   Computed with computeT() from sin.sage, and manually added entry 0. */
static const uint64_t _T[21] = {
  0,
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
   Each entry is to be interpreted as (hi/2^64+lo/2^128)*2^ex*(-1)^sgn.
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

/* Put in Y an approximation of sin2pi(X), for 0 <= X < 2^-11,
   where X2 approximates X^2.
   Absolute error bounded by 2^-132.999 with 0 <= Y < 0.003068
   (see evalPS() in sin.sage), and relative error bounded by
   2^-124.648 (see evalPSrel(K=8) in sin.sage). */
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
    cnt = __builtin_clzll (X->hi);
    if (cnt)
    {
      X->hi = (X->hi << cnt) | (X->lo >> (64 - cnt));
      X->lo = X->lo << cnt;
    }
    X->ex -= cnt;
  }
  else if (X->lo != 0)
  {
    cnt = __builtin_clzll (X->lo);
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
  static const uint64_t *T = _T+1;

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
    e = X->ex;
    normalize (X);
    e = e - X->ex;
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
  if (X->ex <= -11)
    return 0;
  int sh = 64 - 11 - X->ex;
  int i = X->hi >> sh;
  X->hi = X->hi & ((1ull << sh) - 1);
  normalize (X);
  return i;
}

// argument reduction for |x| >= 2^31
// return k and r such that
// x/(2pi) mod 1 = k/2^15 + r + s with 0 <= r < 2^-15 and 0 <= s < 2^-67.988
static uint64_t
reduce_large (double *r, double x)
{
  b64u64_u t = {.f = x};
  int e = (t.u >> 52) & 0x7ff; /* 1054 <= e <= 2046 */
  uint64_t m = (1ull << 52) | (t.u & 0xfffffffffffffull);
  // x = m * 2^(e-1075)
  /* _T[j] corresponds to _T[j]/2^(64*j) thus _T[j]*x corresponds to
     m*_T[j]*2^(e-1075-64*j). To get a non-zero fractional value,
     we need e-1075-64*j < 0, thus 64*j > e-1075 or 64*j >= e-1074. */
  int i = (e - 1011) / 64; // i = ceil((e-1074)/64), 0 <= i <= 16
  int f = (e - 1011) & 0x3f;
  /* the number of fractional bits from m*_T[i] is 64-f, thus we have to
     shift _T[i] by f bits to get 64 fractional bits */
  uint64_t U0, U1;
  if (f == 0) {
    U0 = _T[i];
    U1 = _T[i+1];
  } else {
    U0 = (_T[i] << f) | (_T[i+1] >> (64-f));
    U1 = (_T[i+1] << f) | (_T[i+2] >> (64-f));
  }
  /* Remark: computing directly u with 128-bit arithmetic from _T[i],
     _T[i+1] and _T[i+2] is slower (surely because 128-bit arithmetic is
     emulated.) */
  u128 u = (u128) U1 | (((u128) U0) << 64);
  u = (u128) m * u;
  // round r to nearest, where 0x810000000000000 = 2^59 + 2^52
  static const u128 magic = ((u128) 1 << 112) + 0x810000000000000ull;
  u += magic;
  t.f = (u << 15) >> 75; // next 53 bits of u after the first 15
  *r = t.f * 0x1p-68 - 0x1p-16;
  return u >> 113;
  // since we return 15 bits in i and 53 in h, the accuracy is at most 2^-68
}

/* Assume x is a regular number, and |x| > 0x1.7137449123ef6p-26. */
__attribute__((cold))
static double
sin_accurate (double x)
{
  double absx = (x > 0) ? x : -x;

  dint64_t X[1];
  dint_fromd (X, absx);

  /* reduce argument */
  reduce (X);

  // now |X - x/(2pi) mod 1| < 2^-126.67*X, with 0 <= X < 1.

  int neg = x < 0, is_sin = 1;

  // Write X = i/2^11 + r with 0 <= r < 2^11.
  int i = reduce2 (X); // exact

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
    X->sgn = 1; // negate X
    add_dint (X, &MAGIC, X); // X -> 2^-11 - X
    // here: 256 <= i <= 511
    i = 0x1ff - i;
    // now 0 <= i < 256
  }

  // now 0 <= i < 256 and 0 <= X < 2^-11

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
    // sin2pi(R) ~ sin2pi(i/2^11)*cos2pi(X)+cos2pi(i/2^11)*sin2pi(X)
    mul_dint (U, S+i, U);
    /* since 0 <= S[i] < 0.705 and 0.999 < Uin <= 1, we have
       0 <= U < 0.705 */
    mul_dint (V, C+i, V);
    /* For the error analysis, we distinguish the case i=0.
       For i=0, we have S[i]=0 and C[1]=1, thus V is the value computed
       by evalPS() above, with relative error < 2^-124.648.

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
    // cos2pi(R) ~ cos2pi(i/2^11)*cos2pi(X)-sin2pi(i/2^11)*sin2pi(X)
    mul_dint (U, C+i, U);
    mul_dint (V, S+i, V);
    V->sgn = 1 - V->sgn; // negate V
    /* For 0 <= i < 256, analyze_sin_case2(rel=true) from sin.sage gives a
       relative error bound of -123.540 (obtained for i=0):
       |U - cos2pi(R)| < |U| * 2^-123.540, with U the value computed
       after add_dint (U, U, V) below.

       For the approximation error in R, we have:
       sin |x| = cos2pi (R * (1 + eps))
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
        /* the following worst case was reported by Erik E., it has 68
           identical bits after the round bit */
        {0x1.6ac5b262ca1ffp+849, 0x1p+0, -0x1.2b089ea1e692bp-123},
      };
      for (int j = 0; j < 2; j++)
        {
          if (__builtin_fabs (x) == exceptions[j][0])
            return (x > 0) ? exceptions[j][1] + exceptions[j][2]
              : -exceptions[j][1] - exceptions[j][2];
        }
      /* if we go here, we have a hard-to-round case, but since all hard-to-round
         cases are known and pass all tests, we are ok */
    }

  if (neg)
    U->sgn = 1 - U->sgn;

  double y = dint_tod (U);

  return y;
}

#if 0
/* for each j, 0 <= j <= 128, T1[j] contains sh, sl, ch, cl where
   sh+sl is a double-double approximation of sin(j/2^5) and
   ch+cl is a double-double approximation of cos(j/2^5) */
static const double T1[129][4] = {
  {0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0},
  {0x1.ffeaaaeeee86fp-6, -0x1.cd406fb224ae2p-60, 0x1.ffc00155527d3p-1, -0x1.3b54492d89b5bp-55},
  {0x1.ffaaaeeed4edbp-5, -0x1.2d16d32684b69p-59, 0x1.ff0015549f4d3p-1, 0x1.328387b99426fp-55},
  {0x1.7f701032550e4p-4, 0x1.afc2d1800501ap-60, 0x1.fdc06bf7e6b9bp-1, 0x1.31902b535f8dbp-55},
  {0x1.feaaeee86ee36p-4, -0x1.afcb2bcc6f03bp-59, 0x1.fc015527d5bd3p-1, 0x1.b68f35094efb8p-55},
  {0x1.3eb312c5d66cbp-3, 0x1.47d666b66cb91p-57, 0x1.f9c340a7cc428p-1, 0x1.c5b6b063b7462p-55},
  {0x1.7dc102fbaf2b5p-3, 0x1.5ab50e23c97c3p-59, 0x1.f706bdf9ece1cp-1, -0x1.698c80c36dcb4p-55},
  {0x1.bc6f84edc6199p-3, 0x1.9c1a56a7b0cabp-57, 0x1.f3cc7c3b3d16ep-1, -0x1.21a3ad28a3494p-57},
  {0x1.faaeed4f31577p-3, -0x1.15d88508e32b8p-57, 0x1.f01549f7deea1p-1, 0x1.d3c1e99e5cafdp-55},
  {0x1.1c37d64c6b876p-2, 0x1.46076fe0dcff4p-56, 0x1.ebe214f76efa8p-1, -0x1.02f9f12ba543ep-55},
  {0x1.3ad129769d3d8p-2, 0x1.03d550487839ap-63, 0x1.e733ea0193d4p-1, -0x1.6428b3546ce13p-55},
  {0x1.591bc9fa2f597p-2, 0x1.7c74bac3fe0cbp-57, 0x1.e20bf49acd6c1p-1, -0x1.660aec7ef636bp-58},
  {0x1.7710255764214p-2, -0x1.6ead7314bb6cep-57, 0x1.dc6b7eb995912p-1, 0x1.4b364776dcd35p-58},
  {0x1.94a6be9f546c5p-2, -0x1.69ce13e683f58p-56, 0x1.d653f073e404p-1, -0x1.76236434bec37p-55},
  {0x1.b1d8305321617p-2, -0x1.ae242cb99f519p-56, 0x1.cfc6cfa52ad9fp-1, 0x1.8b5b5508f2a0dp-55},
  {0x1.ce9d2e3d4a51fp-2, -0x1.2fc8a12dae298p-57, 0x1.c8c5bf8ce1a84p-1, 0x1.ab3d1a1590123p-56},
  {0x1.eaee8744b05fp-2, -0x1.789b43c9b027dp-58, 0x1.c1528065b7d5p-1, -0x1.892111312e828p-55},
  {0x1.0362939c69955p-1, -0x1.2d8cd78397b01p-55, 0x1.b96eeef58840ep-1, 0x1.45a3cc78fadep-58},
  {0x1.110d0c4b69c3bp-1, 0x1.d918998809981p-55, 0x1.b11d04162a4c6p-1, 0x1.1dd561efbc0c2p-56},
  {0x1.1e7343236574cp-1, 0x1.22a3fa4f41d5ap-56, 0x1.a85ed4373e02dp-1, 0x1.9be06385ec792p-57},
  {0x1.2b91dea88421ep-1, -0x1.fa371db216abp-55, 0x1.9f368ed912f85p-1, -0x1.1d200c5791606p-55},
  {0x1.386597456282bp-1, -0x1.10fada93b07a8p-56, 0x1.95a67e00cb1fdp-1, -0x1.0befda21f862dp-55},
  {0x1.44eb381cf386bp-1, -0x1.3ed6c1e6a5505p-55, 0x1.8bb105a5dc9p-1, 0x1.863e03e9474c1p-55},
  {0x1.511f9fd7b351cp-1, -0x1.5c0e861c48831p-55, 0x1.8158a31916d5dp-1, -0x1.de8b90b8228dep-57},
  {0x1.5cffc16bf8f0dp-1, 0x1.96cb370eb578ap-55, 0x1.769fec655211fp-1, -0x1.827d5cf8c68c5p-57},
  {0x1.6888a4e134b2fp-1, -0x1.6b7d37644d5e6p-55, 0x1.6b898fa9efb5dp-1, 0x1.15ac786ccf4b2p-56},
  {0x1.73b7680dea578p-1, -0x1.2248306dc12a2p-56, 0x1.6018526f563dfp-1, 0x1.46ca5e0e432dp-55},
  {0x1.7e893f5037959p-1, 0x1.0eefbaa650c4cp-55, 0x1.544f10f592ca5p-1, -0x1.e7ae8e6c7a62fp-55},
  {0x1.88fb7640b8da2p-1, -0x1.49987c11efaa3p-55, 0x1.4830bd7d4ceb3p-1, 0x1.df77ff20d5448p-55},
  {0x1.930b705f9f85ap-1, -0x1.09ae60f413f4p-61, 0x1.3bc05f8b3a656p-1, 0x1.dab7124aa8c6dp-55},
  {0x1.9cb6a9bbce64bp-1, -0x1.4f3e7a32f8d0cp-56, 0x1.2f011326420e4p-1, 0x1.8e30efe9e96c2p-56},
  {0x1.a5fab793d29c8p-1, 0x1.7482b1e8e6d85p-55, 0x1.21f608107e37ap-1, -0x1.0a3f22ad6358p-55},
  {0x1.aed548f090ceep-1, 0x1.06374f484e288p-59, 0x1.14a280fb5068cp-1, -0x1.b71edcc9344bcp-55},
  {0x1.b74427397fca2p-1, 0x1.da351af253ee4p-55, 0x1.0709d2b6b95eep-1, -0x1.71cc4ee678c32p-55},
  {0x1.bf4536c24bb85p-1, 0x1.97632053703fp-55, 0x1.f25ec6b852fc2p-2, 0x1.445cbca9a80a8p-56},
  {0x1.c6d67751be646p-1, 0x1.d163b7b4fe389p-56, 0x1.d62d52e9fdfa9p-2, 0x1.f6eae4ae67d35p-58},
  {0x1.cdf604a1cadcep-1, -0x1.6b50757f2fa4p-56, 0x1.b9865639d0596p-2, -0x1.931bd06786cb9p-56},
  {0x1.d4a216d89c717p-1, 0x1.d4810b29c8736p-55, 0x1.9c70fa40c279dp-2, -0x1.6346cef9b5fa7p-58},
  {0x1.dad902fa8ac87p-1, 0x1.ea5e370875907p-58, 0x1.7ef4842f0bccdp-2, 0x1.83529407722f1p-56},
  {0x1.e0993b54d68f6p-1, -0x1.f26cc0d6a7cecp-58, 0x1.611852fae0769p-2, -0x1.71272938d7ae8p-57},
  {0x1.e5e14fe11418cp-1, 0x1.f26492c1c25ap-57, 0x1.42e3dd88bd952p-2, -0x1.353a9f74bf255p-57},
  {0x1.eaafeea12b0c4p-1, 0x1.d7af5fa4a5c74p-57, 0x1.245eb0cdba154p-2, -0x1.c4555428fdfb4p-57},
  {0x1.ef03e3f3d42a2p-1, 0x1.0572b0573c404p-59, 0x1.05906dec537dap-2, 0x1.12c3f77448473p-61},
  {0x1.f2dc1ae18002ep-1, -0x1.be7521dc7c74p-58, 0x1.cd0190985ef77p-3, -0x1.11be2ffbeed45p-58},
  {0x1.f6379d619369dp-1, 0x1.6b296ac1928abp-55, 0x1.8e6f075a987d6p-3, 0x1.a57e7fd1918d8p-62},
  {0x1.f9159497e853fp-1, 0x1.66c77a4219a37p-56, 0x1.4f78e46e35a46p-3, -0x1.82bbe6c49f2bp-59},
  {0x1.fb75490a83c2cp-1, 0x1.d9fbeed39ae46p-55, 0x1.102ee507ff5fp-3, -0x1.77ec7eee89a9bp-57},
  {0x1.fd5622cf734eap-1, 0x1.576f5c33de713p-55, 0x1.a141b6a6da89dp-4, 0x1.dd0de04944ab6p-58},
  {0x1.feb7a9b2c6d8bp-1, -0x1.0c8f40129a886p-56, 0x1.21bd54fc5f9a7p-4, 0x1.0fcb936b1ce7ep-58},
  {0x1.ff9985549ce69p-1, 0x1.57aa6cfbfc93dp-55, 0x1.43e10afde8436p-5, -0x1.fc499d21a932p-60},
  {0x1.fffb7d3f3a253p-1, -0x1.2d4934e6c1f3dp-56, 0x1.0fd9d5c093df5p-7, -0x1.50076d7383a18p-64},
  {0x1.ffdd78f5268bfp-1, 0x1.f41fc70ae37ddp-56, -0x1.780a3ac0ba58bp-6, 0x1.d5e43e408abb2p-63},
  {0x1.ff3f7ff74c9a7p-1, -0x1.10dae3aca52fep-55, -0x1.bbd1afe4369efp-5, 0x1.50fbc01ce6562p-59},
  {0x1.fe21b9c319278p-1, 0x1.8ac14da77e504p-59, -0x1.5d97a825ea2aap-4, -0x1.72c8c2a1b0d92p-58},
  {0x1.fc846dc89c3afp-1, 0x1.75931f07e378ap-55, -0x1.dcef1441cb33cp-4, -0x1.f2bc7445c5208p-58},
  {0x1.fa680358ad68ap-1, 0x1.89f16c1748c9ap-55, -0x1.2de7a38a3ff6fp-3, 0x1.054bfdacd158ep-59},
  {0x1.f7cd018b18246p-1, -0x1.c06b85582fc39p-56, -0x1.6d0c449d3e98ap-3, -0x1.623c28c417034p-58},
  {0x1.f4b40f1cd6831p-1, 0x1.98c5d3c1c9353p-55, -0x1.abd5a485cce28p-3, -0x1.ebfb11995e71ep-62},
  {0x1.f11df24662dadp-1, -0x1.09b7c1ab8f94bp-56, -0x1.ea34113fa728fp-3, 0x1.abd498353e0e9p-57},
  {0x1.ed0b908a2aac3p-1, -0x1.4ece5211b2c6ap-56, -0x1.140bf9c1636a7p-2, 0x1.4fbce747bfd47p-58},
  {0x1.e87dee7b2f393p-1, -0x1.06241f0ee831p-59, -0x1.32b8e9548fce1p-2, 0x1.3fc0930cc38b6p-56},
  {0x1.e3762f7be2204p-1, -0x1.0272412ab7375p-55, -0x1.51192c465a31bp-2, -0x1.053ee416dfe5ap-56},
  {0x1.ddf595754e444p-1, -0x1.4ce8990cb150ep-56, -0x1.6f252aae8625bp-2, 0x1.ae75f52c15a19p-57},
  {0x1.d7fd80869f372p-1, -0x1.c342d6d256f85p-57, -0x1.8cd561b589476p-2, -0x1.acf78510604dap-59},
  {0x1.d18f6ead1b446p-1, -0x1.02a3dbf3bffb2p-56, -0x1.aa22657537205p-2, 0x1.6f3341d4d1235p-56},
  {0x1.caacfb64a61cdp-1, -0x1.fbf52442206c4p-56, -0x1.c704e2d3b0cbfp-2, 0x1.0908c2140ecf5p-60},
  {0x1.c357df40e4024p-1, -0x1.f162bd32468fep-56, -0x1.e375a15821ab9p-2, -0x1.a0e030d758208p-59},
  {0x1.bb91ef7f1729ep-1, 0x1.ba36b4a8034e5p-59, -0x1.ff6d84f8d3facp-2, -0x1.b3aa6bb754ef4p-59},
  {0x1.b35d1d90d2dd6p-1, -0x1.d3d716afba31dp-57, -0x1.0d72c7f114e12p-1, 0x1.6788abb417645p-55},
  {0x1.aabb769fa1ad3p-1, 0x1.ead5c74acefc3p-55, -0x1.1aeb721b04367p-1, -0x1.4ee940f7119e4p-56},
  {0x1.a1af2309bdca6p-1, -0x1.8b169e843eaf8p-55, -0x1.281d62e1a3938p-1, 0x1.6a2cae7608016p-55},
  {0x1.983a65d7fc58p-1, 0x1.d8dba65860c9p-55, -0x1.35054dda59168p-1, -0x1.664c0a672acb8p-55},
  {0x1.8e5f9c2d0e3a9p-1, 0x1.5dc0da4ffdf4ep-55, -0x1.419ff91b9ba6dp-1, 0x1.9a10a4b5cbe7ep-55},
  {0x1.84213cae3a92p-1, 0x1.298047b6629bap-55, -0x1.4dea3e0b69097p-1, -0x1.2bc301ec35804p-55},
  {0x1.7981d6e5b8b11p-1, -0x1.9fcdb3acf5b7p-57, -0x1.59e10a28e82edp-1, 0x1.f53d598593a6cp-57},
  {0x1.6e84129ed0f95p-1, 0x1.a56bab25774afp-55, -0x1.65815fd1054fdp-1, -0x1.a156030f696b6p-55},
  {0x1.632aaf3bed93bp-1, 0x1.0637f900540a7p-60, -0x1.70c856fdd6b67p-1, 0x1.a18459c4d6abdp-55},
  {0x1.57788306c57f6p-1, 0x1.a7131e3be9006p-56, -0x1.7bb31e009a57bp-1, 0x1.541fc31d208bdp-55},
  {0x1.4b707a7acdecdp-1, -0x1.ef71ae7061d34p-55, -0x1.863efa361dc25p-1, -0x1.5e50f57769cbap-56},
  {0x1.3f15978a1f45fp-1, -0x1.be1f86c7149adp-56, -0x1.906948b56347dp-1, 0x1.26b777679a478p-57},
  {0x1.326af0dcfcab1p-1, -0x1.fd42734161659p-55, -0x1.9a2f7ef858b7dp-1, -0x1.587cfaa17e973p-56},
  {0x1.2573b10c2dffep-1, 0x1.0cb85186507c5p-56, -0x1.a38f2b7e75819p-1, 0x1.bd5e7c6d218f8p-57},
  {0x1.183315d65df2ap-1, -0x1.41089cbc8c0afp-55, -0x1.ac85f6691793ep-1, 0x1.eb962bc7b74ap-55},
  {0x1.0aac6f50aea35p-1, -0x1.49fd3bc15c939p-55, -0x1.b511a21177e5ep-1, -0x1.75f0809e1e829p-55},
  {0x1.f9c63e25718c7p-2, -0x1.da7d3b28b8de6p-58, -0x1.bd300b98112c3p-1, -0x1.0e2cbb26ca4edp-55},
  {0x1.ddb52ebc547f7p-2, 0x1.8b4ca4f49f731p-56, -0x1.c4df2b6d54e0cp-1, 0x1.f42713219f479p-55},
  {0x1.c12cb48474a24p-2, -0x1.7eea8e847d17dp-56, -0x1.cc1d15d38c71cp-1, -0x1.6b76b64db6c33p-55},
  {0x1.a433f17654f04p-2, -0x1.8273ee47f959dp-56, -0x1.d2e7fb59c6201p-1, -0x1.106e2c45a122ep-56},
  {0x1.86d2239c183fbp-2, 0x1.f838db9ee6256p-56, -0x1.d93e294faed14p-1, 0x1.421d74d654ed8p-56},
  {0x1.690ea3420861p-2, -0x1.5c3804d08d097p-56, -0x1.df1e0a323be1p-1, -0x1.f8360382131eep-55},
  {0x1.4af0e1208cd6dp-2, 0x1.4923b3ae7090ap-56, -0x1.e486261109c75p-1, -0x1.e72962145517bp-59},
  {0x1.2c80648006a85p-2, 0x1.c9458401665b5p-58, -0x1.e97522ec563bcp-1, 0x1.35dac6006c32ap-55},
  {0x1.0dc4c95708521p-2, 0x1.4fefad09e5717p-60, -0x1.ede9c50b7e58fp-1, -0x1.739952d0f281fp-57},
  {0x1.dd8b7cc6c48dbp-3, 0x1.20505b9f3773bp-57, -0x1.f1e2ef4beb207p-1, 0x1.b44f6d483c9bcp-55},
  {0x1.9f16067cfb738p-3, 0x1.4786db3b8ead4p-57, -0x1.f55fa36858a4p-1, 0x1.b5642982a1298p-55},
  {0x1.6038ccdb01312p-3, -0x1.fe5f02cef39abp-60, -0x1.f85f02386603dp-1, -0x1.178460cf1ed29p-58},
  {0x1.210386db6d55bp-3, 0x1.3c7205d08d063p-57, -0x1.fae04be85e5d2p-1, -0x1.83effc17efb54p-55},
  {0x1.c30c02f6f2e41p-4, 0x1.27df80431e208p-61, -0x1.fce2e0292cb7bp-1, 0x1.08f56002d0a5ep-56},
  {0x1.43a0378fadb65p-4, 0x1.7317f6e0fc189p-59, -0x1.fe663e586ef52p-1, 0x1.44a72b25b459cp-55},
  {0x1.87c70b94029d7p-5, -0x1.fcdc8b319b851p-62, -0x1.ff6a05a09dbe2p-1, -0x1.0dbce2e0658e1p-55},
  {0x1.0fd770a03e5aap-6, -0x1.96353881cf537p-60, -0x1.ffedf51141634p-1, 0x1.e060226d9f29ep-59},
  {-0x1.e04654b27e08ap-7, 0x1.a30a09ec6a024p-66, -0x1.fff1ebaf2da3fp-1, -0x1.f5e622c0e6966p-55},
  {-0x1.77f0dee42925cp-5, -0x1.cc6e70c125987p-59, -0x1.ff75e87cc04ep-1, -0x1.1093c3d953238p-55},
  {-0x1.3bb9172c9b5d8p-4, 0x1.74e861f4eff6cp-59, -0x1.fe7a0a7a20a48p-1, -0x1.385c8f10b6ed5p-56},
  {-0x1.bb2ad2464a48cp-4, -0x1.62baeb29e6797p-58, -0x1.fcfe909d7f7f8p-1, 0x1.3f803163b746p-55},
  {-0x1.1d16e27d233cp-3, 0x1.af43a26adac33p-58, -0x1.fb03d9c35a13ap-1, 0x1.1ed7613b89931p-56},
  {-0x1.5c51179a9d633p-3, -0x1.bd29dae986182p-60, -0x1.f88a6496c3517p-1, 0x1.3d43cd3a4b0f7p-57},
  {-0x1.9b343a429923cp-3, -0x1.418d85d629245p-57, -0x1.f592cf71b9c97p-1, 0x1.8b7251336dee6p-56},
  {-0x1.d9b09200454f7p-3, -0x1.a6111f33eb61cp-58, -0x1.f21dd83591ff9p-1, 0x1.494aa3fd99a7ap-57},
  {-0x1.0bdb4008811f4p-2, 0x1.c7fadfa47ac74p-56, -0x1.ee2c5c1b7f135p-1, 0x1.1833162040694p-55},
  {-0x1.2a9b41a5fed1fp-2, 0x1.5ee1f3a1c3d2cp-57, -0x1.e9bf577d4599dp-1, 0x1.d6e9391acc89ap-55},
  {-0x1.49109e01340b2p-2, 0x1.6999aa9e231c6p-56, -0x1.e4d7e596267d2p-1, -0x1.fbf25c49fbb5ap-55},
  {-0x1.6733b7eba621fp-2, -0x1.ae055844cf8c8p-57, -0x1.df77403c11a5fp-1, 0x1.094dd04296f85p-58},
  {-0x1.84fd06c708f17p-2, -0x1.301034b191101p-57, -0x1.d99ebf9132218p-1, 0x1.f7d3698e47beap-57},
  {-0x1.a26518675c6p-2, 0x1.aba1272dd6db8p-56, -0x1.d34fd9ade7622p-1, 0x1.e93a474b00113p-56},
  {-0x1.bf6492ef6d71ep-2, 0x1.3cf40fc46e3b1p-57, -0x1.cc8c22434119ep-1, -0x1.b225bfa32bdaap-55},
  {-0x1.dbf436a743c91p-2, -0x1.28c5b433b8062p-56, -0x1.c5554a3615112p-1, 0x1.39d87639a31a8p-58},
  {-0x1.f80cdfcc05f91p-2, -0x1.874f2a047d009p-56, -0x1.bdad1f32c831ep-1, -0x1.9c2ff5dd5723bp-55},
  {-0x1.09d3c42c705c2p-1, -0x1.b758d2b662c18p-56, -0x1.b5958b39e5d69p-1, 0x1.e9b5878af2346p-56},
  {-0x1.175ea4e43f5bcp-1, -0x1.ee921ee02c7fdp-57, -0x1.ad109425a2341p-1, 0x1.5779694d7752ap-55},
  {-0x1.24a3af6750621p-1, -0x1.a3d145c0f88eap-55, -0x1.a4205b28667f7p-1, 0x1.431eff5650152p-55},
  {-0x1.319f9284b3e88p-1, 0x1.57c0837231309p-55, -0x1.9ac71c44872a6p-1, 0x1.5bb621991534fp-61},
  {-0x1.3e4f0f54f24aap-1, 0x1.58c931643b365p-55, -0x1.91072dbd4648dp-1, -0x1.ea865c0f57ea9p-58},
  {-0x1.4aaefa09c1509p-1, 0x1.b53f2c50903ep-55, -0x1.86e2ff8145dep-1, 0x1.39f6c59fb1358p-56},
  {-0x1.56bc3ab8f386fp-1, 0x1.1ff11e3bc3a75p-56, -0x1.7c5d1a8e8f73ep-1, 0x1.b79e386300bd8p-57},
  {-0x1.6273ce226eaa8p-1, -0x1.d62d9b9afe29cp-55, -0x1.7178205057fa1p-1, 0x1.0b622467e57cdp-56},
  {-0x1.6dd2c670f7aa7p-1, -0x1.8b88cd0857facp-59, -0x1.6636c9f6a87a7p-1, 0x1.72233bab9ac71p-55},
  {-0x1.78d64bf5a40f2p-1, -0x1.14c27294a56bbp-55, -0x1.5a9be7c815b86p-1, 0x1.9f010f4932d68p-58},
  {-0x1.837b9dddc1eaep-1, -0x1.c33a601568391p-55, -0x1.4eaa606db24c1p-1, 0x1.dcc92f1e91c23p-56},
};

/* for each j, 0 <= j < 128, T2[j] contains sh, sl, ch, cl where
   sh+sl is a double-double approximation of sin(j/2^12) and
   ch+cl is a double-double approximation of cos(j/2^12) */
static const double T2[128][4] = {
  {0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0},
  {0x1.ffffffaaaaaabp-13, -0x1.1111112b12b13p-69, 0x1.ffffff0000001p-1, 0x1.55555527d27d3p-55},
  {0x1.fffffeaaaaaafp-12, -0x1.1111179179173p-68, 0x1.fffffc0000015p-1, 0x1.555549f49f4acp-55},
  {0x1.7ffffdc00001p-11, 0x1.99997dd41d455p-66, 0x1.fffff7000006cp-1, -0x1.033333098af8bp-72},
  {0x1.fffffaaaaaaefp-11, -0x1.1112b12b1254bp-67, 0x1.fffff00000155p-1, 0x1.55527d27d34d3p-55},
  {0x1.3ffffacaaab13p-10, -0x1.5557455d752b2p-65, 0x1.ffffe70000341p-1, 0x1.554a7b8e3dbbap-55},
  {0x1.7ffff70000103p-10, 0x1.9992a83a8720fp-65, 0x1.ffffdc00006cp-1, -0x1.0333328c92496p-66},
  {0x1.bffff1b555786p-10, -0x1.bbc5f227cb89ep-64, 0x1.ffffcf0000c81p-1, 0x1.5503a1f4e6c6fp-55},
  {0x1.ffffeaaaaaeefp-10, -0x1.117917911cap-66, 0x1.ffffc00001555p-1, 0x1.549f49f56f56fp-55},
  {0x1.1ffff0d0003d8p-9, 0x1.32f7e32c25783p-64, 0x1.ffffaf000222cp-1, -0x1.710e645096267p-63},
  {0x1.3fffeb2aab12dp-9, 0x1.551754519b323p-63, 0x1.ffff9c0003415p-1, 0x1.529ee39310f7ep-55},
  {0x1.5fffe44555fd2p-9, -0x1.ef67c2f29c278p-63, 0x1.ffff870004c41p-1, 0x1.5087153234b58p-55},
  {0x1.7fffdc0001033p-9, 0x1.97dd41d795f16p-64, 0x1.ffff700006cp-1, -0x1.03333098af8f2p-60},
  {0x1.9fffd23aac2d7p-9, -0x1.e3f22117ed7e3p-65, 0x1.ffff5700094c1p-1, 0x1.483d621c22ffp-55},
  {0x1.bfffc6d557859p-9, 0x1.06daa5161cae9p-65, 0x1.ffff3c000c815p-1, 0x1.40e87d6f4f713p-55},
  {0x1.dfffb9b00317p-9, 0x1.f7b9356301cb9p-64, 0x1.ffff1f00107acp-1, -0x1.ee62783da208p-59},
  {0x1.ffffaaaaaeeefp-9, -0x1.2b12b0ce9b237p-65, 0x1.ffff000015555p-1, 0x1.27d27df7df7bbp-55},
  {0x1.0fffccd2ad8e3p-8, -0x1.9278cb9abfc78p-63, 0x1.fffedf001b301p-1, 0x1.13db234688337p-55},
  {0x1.1fffc34003d82p-8, 0x1.922f98d0e7fc7p-62, 0x1.fffebc00222cp-1, -0x1.710e5e0f257d2p-57},
  {0x1.2fffb88d5a5efp-8, 0x1.eca449a9e19c3p-62, 0x1.fffe97002a6c1p-1, 0x1.ab6d30bd0832dp-56},
  {0x1.3fffacaab12d5p-8, 0x1.45d514a762ee1p-62, 0x1.fffe700034155p-1, 0x1.4f71d0cc9a3e9p-56},
  {0x1.4fff9f88084f2p-8, 0x1.3ac6f0be63e61p-63, 0x1.fffe47003f4ecp-1, -0x1.d14f8e7c7a202p-56},
  {0x1.5fff91155fd18p-8, 0x1.e5b8217c1c4b8p-63, 0x1.fffe1c004c415p-1, 0x1.0e2aa2b6bba3ap-58},
  {0x1.6fff8142b7c2fp-8, -0x1.7499649a6260cp-63, 0x1.fffdef005b181p-1, -0x1.e1ea79cf95755p-58},
  {0x1.7fff700010333p-8, 0x1.2a83abb33321p-63, 0x1.fffdc0006bfffp-1, 0x1.f999ae6db6562p-55},
  {0x1.8fff5d3d69339p-8, 0x1.e0c0912734b74p-62, 0x1.fffd8f007f281p-1, -0x1.40f0a79254f2fp-55},
  {0x1.9fff48eac2d6ep-8, 0x1.3a22bca3cb044p-65, 0x1.fffd5c0094c15p-1, -0x1.f0a75b5479558p-55},
  {0x1.afff32f81d316p-8, 0x1.9ee4ff27a6b5fp-64, 0x1.fffd2700acfebp-1, -0x1.af1cc288ce35dp-59},
  {0x1.bfff1b5578591p-8, -0x1.7c89d8d240eeap-64, 0x1.fffcf000c8154p-1, 0x1.d0fc8b8c8bd4p-58},
  {0x1.cfff01f2d4659p-8, -0x1.2ec3de5bb54ffp-62, 0x1.fffcb700e63cp-1, -0x1.f074185fe8d96p-56},
  {0x1.dffee6c031704p-8, -0x1.08d949ecd162fp-62, 0x1.fffc7c0107abep-1, 0x1.19d9f0976f758p-57},
  {0x1.effec9ad8f945p-8, 0x1.10af5f369e5c9p-62, 0x1.fffc3f012c9ffp-1, -0x1.22c0975a67b07p-59},
  {0x1.fffeaaaaeeeefp-8, -0x1.e45e2ec67b77cp-62, 0x1.fffc000155552p-1, 0x1.f4a01a0196daep-55},
  {0x1.07ff44d427cf8p-7, -0x1.d507d1ffd9d83p-63, 0x1.fffbbf01820a9p-1, -0x1.af53b77a5f1d7p-55},
  {0x1.0fff334ad8e2dp-7, -0x1.824c979587b6cp-61, 0x1.fffb7c01b3011p-1, 0x1.ed939e215e812p-56},
  {0x1.17ff20b18ac2ep-7, 0x1.9459e62d75e19p-63, 0x1.fffb3701e87bcp-1, 0x1.dab9a5a85447cp-55},
  {0x1.1fff0d003d826p-7, -0x1.0399fe343032fp-63, 0x1.fffaf00222bfap-1, 0x1.de375ed37acdap-56},
  {0x1.27fef82ef134fp-7, 0x1.09b6df6985f25p-61, 0x1.fffaa7026213bp-1, -0x1.daa28568eb4acp-55},
  {0x1.2ffee235a5ef7p-7, 0x1.8558598d90cd7p-62, 0x1.fffa5c02a6c0dp-1, 0x1.6da880a614dfbp-55},
  {0x1.37fecb0c5bc7dp-7, 0x1.6774fe4973282p-61, 0x1.fffa0f02f1123p-1, -0x1.493019b01ecc7p-55},
  {0x1.3ffeb2ab12d54p-7, 0x1.75456a6f1b29dp-61, 0x1.fff9c0034154ap-1, 0x1.ee3dbba2340adp-55},
  {0x1.47fe9909cb302p-7, 0x1.69043dcc5434fp-61, 0x1.fff96f0397d75p-1, -0x1.002043ee09aa1p-55},
  {0x1.4ffe7e2084f21p-7, 0x1.bf44e2651ff34p-61, 0x1.fff91c03f4eb1p-1, 0x1.d6138e8f53867p-55},
  {0x1.57fe61e74036p-7, 0x1.e881ee0b60c64p-62, 0x1.fff8c70458e31p-1, -0x1.a66d7a64a2689p-55},
  {0x1.5ffe4455fd182p-7, 0x1.83d1949b841dep-61, 0x1.fff87004c4142p-1, 0x1.c5737d7d55d85p-57},
  {0x1.67fe2564bbb61p-7, -0x1.697c0711cd593p-65, 0x1.fff8170536d56p-1, 0x1.5242e2d644432p-62},
  {0x1.6ffe050b7c2ebp-7, 0x1.400f063986d84p-62, 0x1.fff7bc05b17fcp-1, 0x1.e16e9d29a9286p-56},
  {0x1.77fde3423ea27p-7, -0x1.ddc4cb4ee3ee3p-61, 0x1.fff75f06346e5p-1, -0x1.c62e88e8187bp-56},
  {0x1.7ffdc0010333p-7, -0x1.15efa2be503dbp-61, 0x1.fff70006bffep-1, -0x1.9984c57e6cfb8p-55},
  {0x1.87fd9b3fca03bp-7, -0x1.1de672ad54c63p-62, 0x1.fff69f07548ddp-1, -0x1.55771a7ba26adp-55},
  {0x1.8ffd74f693394p-7, 0x1.8131581e4c084p-64, 0x1.fff63c07f27ecp-1, -0x1.e0a1e80fb9327p-58},
  {0x1.97fd4d1d5efap-7, -0x1.811befd7ce4a2p-63, 0x1.fff5d7089a33dp-1, 0x1.8b2585d87c60ap-55},
  {0x1.9ffd23ac2d6dcp-7, 0x1.bbf7c37af1411p-64, 0x1.fff570094c121p-1, -0x1.4dc992d56e09p-58},
  {0x1.a7fcf89afebep-7, -0x1.917b1ac1263ebp-61, 0x1.fff5070a08807p-1, -0x1.872a193837753p-55},
  {0x1.affccbe1d315cp-7, -0x1.44df7ec641172p-61, 0x1.fff49c0acfe7ep-1, 0x1.43b50aa13d353p-55},
  {0x1.b7fc9d78aaa1cp-7, -0x1.b4699f263634dp-63, 0x1.fff42f0ba2b38p-1, 0x1.4207f0137ff33p-61},
  {0x1.bffc6d5785907p-7, -0x1.2aca46ff7e621p-62, 0x1.fff3c00c81504p-1, -0x1.77e605ed492eap-55},
  {0x1.c7fc3b766411fp-7, -0x1.49d1b3c8125c2p-61, 0x1.fff34f0d6c2d1p-1, 0x1.2a823634124e8p-56},
  {0x1.cffc07cd46581p-7, 0x1.e8a928b51c269p-61, 0x1.fff2dc0e63bbp-1, 0x1.f1c3f22357e3bp-55},
  {0x1.d7fbd2542c96ap-7, -0x1.17b9042ec78eep-62, 0x1.fff2670f686d2p-1, -0x1.a043637e0cc6ep-55},
  {0x1.dffb9b031702fp-7, 0x1.c9b737bfa8c4dp-61, 0x1.fff1f0107ab84p-1, 0x1.9dfc25cd50f26p-55},
  {0x1.e7fb61d205d47p-7, 0x1.5972c38824468p-61, 0x1.fff177119b139p-1, -0x1.30256fbb84521p-56},
  {0x1.effb26b8f9445p-7, -0x1.6db4c4aa1f282p-61, 0x1.fff0fc12c9f7fp-1, -0x1.1512a7311528fp-56},
  {0x1.f7fae9aff18d9p-7, -0x1.df1e696fe6169p-65, 0x1.fff07f1407e06p-1, 0x1.7adcaa95bdce5p-55},
  {0x1.fffaaaaeeeed5p-7, -0x1.2ab639a9f0776p-63, 0x1.fff000155549fp-1, 0x1.28a28a03a5ef3p-55},
  {0x1.03fd34d6f8d14p-6, 0x1.3bb4d755ff497p-60, 0x1.ffef7f16b2b3ap-1, -0x1.d494be4afb1dp-55},
  {0x1.07fd13527cf72p-6, 0x1.249fa3c74f086p-62, 0x1.ffeefc18209e5p-1, 0x1.5ecdc572cc8d7p-58},
  {0x1.0bfcf0c60409cp-6, 0x1.ba6cc916b8704p-63, 0x1.ffee77199f8d2p-1, -0x1.318f79b3d580ap-55},
  {0x1.0ffccd2d8e2bbp-6, 0x1.cdaee867eed8ep-63, 0x1.ffedf01b3004fp-1, 0x1.b371329e85efbp-55},
  {0x1.13fca8851b809p-6, -0x1.7cd32fd54abffp-60, 0x1.ffed671cd28cep-1, 0x1.daebdb230c8bap-57},
  {0x1.17fc82c8ac2cfp-6, 0x1.45b1ce04c2096p-60, 0x1.ffecdc1e87adep-1, -0x1.5057029a5a134p-55},
  {0x1.1bfc5bf44056bp-6, -0x1.3646e0544ccf7p-62, 0x1.ffec4f204ff2ep-1, -0x1.e1fba1ff33724p-60},
  {0x1.1ffc3403d8249p-6, -0x1.0653aa527c178p-60, 0x1.ffebc0222be8fp-1, -0x1.bc1e4e7f0892ap-58},
  {0x1.23fc0af373be8p-6, -0x1.6a65ff5086289p-61, 0x1.ffeb2f241c1fp-1, 0x1.bf4cb8ada359ap-55},
  {0x1.27fbe0bf134d9p-6, 0x1.a87c8bd5fe0d7p-61, 0x1.ffea9c2621262p-1, 0x1.59511f26c46p-55},
  {0x1.2bfbb562b6fcp-6, 0x1.86e1aa36412bcp-61, 0x1.ffea07283b915p-1, -0x1.eabb6df3bcd03p-55},
  {0x1.2ffb88da5ef53p-6, -0x1.bb327ff45e021p-60, 0x1.ffe9702a6bf57p-1, -0x1.26eda3fa8d528p-56},
  {0x1.33fb5b220b659p-6, -0x1.f1566bc83e5p-62, 0x1.ffe8d72cb2e99p-1, 0x1.a60588d2737e2p-56},
  {0x1.37fb2c35bc7afp-6, -0x1.55c0339f44b0cp-60, 0x1.ffe83c2f1106bp-1, 0x1.b6f11be57003bp-55},
  {0x1.3bfafc1172643p-6, -0x1.1e63460f9850dp-60, 0x1.ffe79f3186e7dp-1, 0x1.80c9b0f563d98p-55},
  {0x1.3ffacab12d517p-6, 0x1.519b3218acccfp-60, 0x1.ffe700341529fp-1, -0x1.b3bc25e3e4cb3p-57},
  {0x1.43fa9810ed743p-6, -0x1.df2bcf989a364p-60, 0x1.ffe65f36bc6cp-1, -0x1.6b99d2d499d3dp-56},
  {0x1.47fa642cb2feep-6, 0x1.74934f60a45e3p-60, 0x1.ffe5bc397d4fp-1, -0x1.d1e6991ea226bp-62},
  {0x1.4bfa2f007e259p-6, 0x1.c328811043de7p-62, 0x1.ffe5173c5875fp-1, 0x1.f2a04a6ccd9dap-56},
  {0x1.4ff9f8884f1d6p-6, -0x1.c7fcce2c5b193p-60, 0x1.ffe4703f4e85dp-1, 0x1.8a41c0e1e6f4dp-55},
  {0x1.53f9c0c0261cbp-6, 0x1.d6a98a37972edp-61, 0x1.ffe3c7426025ap-1, 0x1.e3532548e408bp-56},
  {0x1.57f987a4035b7p-6, -0x1.55acd0c307e42p-60, 0x1.ffe31c458dfe6p-1, -0x1.94e3de2e570cap-55},
  {0x1.5bf94d2fe712ap-6, -0x1.6cfd18c943537p-60, 0x1.ffe26f48d8bafp-1, 0x1.2b1fb1010dd4fp-55},
  {0x1.5ff9115fd17cbp-6, 0x1.c1ca3eb4720c4p-60, 0x1.ffe1c04c41087p-1, 0x1.7c05fb660499bp-57},
  {0x1.63f8d42fc2d59p-6, 0x1.458b136df462cp-62, 0x1.ffe10f4fc795dp-1, -0x1.a6aaf2380fa88p-56},
  {0x1.67f8959bbb5a6p-6, -0x1.a535ade4bedbp-60, 0x1.ffe05c536d14p-1, 0x1.64e784cb2f31cp-56},
  {0x1.6bf8559fbb49ap-6, 0x1.ea6777fd9c253p-61, 0x1.ffdfa75732361p-1, 0x1.be7a505792261p-61},
  {0x1.6ff81437c2e37p-6, -0x1.94119fff4d189p-61, 0x1.ffdef05b17b0fp-1, 0x1.c785970c61f95p-58},
  {0x1.73f7d15fd2692p-6, -0x1.9efb9742db439p-61, 0x1.ffde375f1e3bap-1, 0x1.3bd65ca3cae0bp-57},
  {0x1.77f78d13ea1d9p-6, -0x1.09b8b20c20137p-60, 0x1.ffdd7c63468f2p-1, -0x1.713be1274b02ep-56},
  {0x1.7bf747500a45p-6, 0x1.96ad4a0f9ac82p-60, 0x1.ffdcbf6791666p-1, 0x1.81d57742e4de6p-59},
  {0x1.7ff7001033255p-6, 0x1.efe2b51527336p-64, 0x1.ffdc006bff7e6p-1, 0x1.ae6dae86977bdp-55},
  {0x1.83f6b7506505bp-6, -0x1.65aa4779d09b6p-60, 0x1.ffdb3f7091963p-1, -0x1.1136e785c4576p-55},
  {0x1.87f66d0ca02edp-6, 0x1.118fac888367fp-60, 0x1.ffda7c75486ebp-1, -0x1.4b5a612e3549dp-55},
  {0x1.8bf62140e4eb1p-6, 0x1.166b0c4334226p-61, 0x1.ffd9b77a24caep-1, -0x1.6db770fb99e48p-58},
  {0x1.8ff5d3e933863p-6, 0x1.6d7eddbb6e9a7p-65, 0x1.ffd8f07f276fcp-1, 0x1.109849f8206c1p-55},
  {0x1.93f585018c4d8p-6, 0x1.93dff04353b51p-60, 0x1.ffd8278451245p-1, 0x1.2ee5e89f9a0aep-55},
  {0x1.97f53485ef9p-6, -0x1.b89965ca9b50cp-61, 0x1.ffd75c89a2b19p-1, -0x1.1d3f662bcac2dp-55},
  {0x1.9bf4e2725d9e1p-6, -0x1.225ba014ad163p-62, 0x1.ffd68f8f1ce26p-1, 0x1.af81aaa794ce2p-56},
  {0x1.9ff48ec2d6c9dp-6, 0x1.2345af9bbd83ap-62, 0x1.ffd5c094c083dp-1, 0x1.af561ff13a75fp-55},
  {0x1.a3f439735b66fp-6, 0x1.9de9de9ed472ap-61, 0x1.ffd4ef9a8e64ep-1, 0x1.2799e607422ddp-64},
  {0x1.a7f3e27febcacp-6, 0x1.3dffe805845f6p-60, 0x1.ffd41ca087568p-1, -0x1.a802c083cd99dp-55},
  {0x1.abf389e4884c4p-6, -0x1.78d318ce2bb7dp-60, 0x1.ffd347a6ac2bap-1, -0x1.915ca8003afebp-56},
  {0x1.aff32f9d3143fp-6, -0x1.9d43cec919526p-62, 0x1.ffd270acfdb94p-1, 0x1.155759bf841abp-55},
  {0x1.b3f2d3a5e70c3p-6, -0x1.a26421e02919p-60, 0x1.ffd197b37cd67p-1, -0x1.a2599babe60a8p-55},
  {0x1.b7f275faaa00ep-6, 0x1.7ba9717750096p-61, 0x1.ffd0bcba2a5cp-1, 0x1.7072fd42a5c89p-55},
  {0x1.bbf216977a7fcp-6, 0x1.9defcf844df3p-60, 0x1.ffcfdfc107251p-1, 0x1.44ae5b37f12a4p-56},
  {0x1.bff1b57858e83p-6, 0x1.df20c232a4f03p-60, 0x1.ffcf00c8140e9p-1, -0x1.c3e3e857e2d3p-55},
  {0x1.c3f15299459b5p-6, 0x1.5f20c3debd91bp-60, 0x1.ffce1fcf51f76p-1, 0x1.ea0ae2764eac2p-57},
  {0x1.c7f0edf640fcp-6, -0x1.fa6b3ccbef77ap-66, 0x1.ffcd3cd6c1c09p-1, 0x1.8e0c2a64951abp-55},
  {0x1.cbf0878b4b6eep-6, -0x1.7f837e423c83bp-60, 0x1.ffcc57de644d2p-1, -0x1.94fa6a65a3c6ep-57},
  {0x1.cff01f54655a5p-6, -0x1.a7b19be514189p-63, 0x1.ffcb70e63a81fp-1, 0x1.6ff504d48ec22p-56},
  {0x1.d3efb54d8f269p-6, 0x1.a6967f983a92dp-60, 0x1.ffca87ee45461p-1, -0x1.0fda8383c71c8p-55},
  {0x1.d7ef4972c93dbp-6, 0x1.de02a8fba1a93p-60, 0x1.ffc99cf685826p-1, 0x1.0231cb91e8c08p-57},
  {0x1.dbeedbc0140b9p-6, -0x1.8fc8d2dd5dcc8p-61, 0x1.ffc8affefc21fp-1, -0x1.66e494d6af995p-55},
  {0x1.dfee6c316ffddp-6, -0x1.887f20baa1b72p-60, 0x1.ffc7c107aa11ap-1, -0x1.1ed367057f9fdp-58},
  {0x1.e3edfac2dd84p-6, -0x1.4aee139779fb7p-61, 0x1.ffc6d01090407p-1, 0x1.de2a03070ecbbp-55},
  {0x1.e7ed87705d0f9p-6, 0x1.9b2546d734038p-60, 0x1.ffc5dd19af9f7p-1, -0x1.9a68004756a5ap-55},
  {0x1.ebed1235ef13ep-6, 0x1.eec5d934dc983p-60, 0x1.ffc4e82309217p-1, -0x1.aada86454d932p-56},
  {0x1.efec9b0f94063p-6, -0x1.f38705d2c8b83p-61, 0x1.ffc3f12c9dbb7p-1, 0x1.d6b49f375683ap-55},
  {0x1.f3ec21f94c5d9p-6, -0x1.726084bda1122p-60, 0x1.ffc2f8366e648p-1, 0x1.347fbe10b0405p-61},
  {0x1.f7eba6ef18932p-6, -0x1.d3d43446c6013p-60, 0x1.ffc1fd407c158p-1, -0x1.7e892bb34cd41p-56},
  {0x1.fbeb29ecf921ep-6, 0x1.084d9ff5acb2p-61, 0x1.ffc1004ac7c96p-1, 0x1.06ff6710a573ep-55},
};
#endif

static inline double fasttwosum(double x, double y, double *e){
  double s = x + y, z = s - x;
  *e = y - z;
  return s;
}

// from acos.c (see comments there)
static inline double twosum(double a, double b, double *t){
  double s = a + b;
  double a_prime = s - b;
  double b_prime = s - a_prime;
  double delta_a = a - a_prime;
  double delta_b = b - b_prime;
  *t = delta_a + delta_b;
  return s;
}

static inline double fastsum(double xh, double xl, double yh, double yl, double *e){
  double sl, sh = fasttwosum(xh, yh, &sl);
  *e = (xl + yl) + sl;
  return sh;
}

static inline double muldd(double xh, double xl, double ch, double cl, double *l){
  double ahhh = xh*ch;
  *l = (xh*cl + xl*ch) + __builtin_fma(xh, ch, -ahhh);
  return ahhh;
}

// accurate version, with an extra normalization step
static inline double muldd_acc(double xh, double xl, double ch, double cl, double *l){
  double ahhh = xh*ch;
  *l = (xh*cl + xl*ch) + __builtin_fma(xh, ch, -ahhh);
  return fasttwosum (ahhh, *l, l);
}

static inline double mulddd(double xh, double xl, double ch, double *l){
  double ahhh = xh*ch;
  *l = xl*ch + __builtin_fma(xh, ch, -ahhh);
  return ahhh;
}

// accurate version, with an extra normalization step
static inline double mulddd_acc(double xh, double xl, double ch, double *l){
  double ahhh = xh*ch;
  *l = xl*ch + __builtin_fma(xh, ch, -ahhh);
  return fasttwosum (ahhh, *l, l);
}

#if 0
static double __attribute__((noinline))
as_sin_fast_acc (double x)
{
  b64u64_u t = {.f = x};
  int sgn = t.u >> 63; // save sign
  t.u &= 0x7fffffffffffffffull; // t.u is now the encoding of |x|
  double ax = __builtin_fabs(x), s = 0x1p+12 * ax;
  double jd = __builtin_roundeven (s); // 0 <= j <= 2^14
  double r = (s - jd) * 0x1p-12; // |r| <= 2^-13
  double r2h = r * r, r2l = __builtin_fma(r,r,-r2h);
  int j = jd, i1 = j >> 7, i2 = j & 0x7f;
  double s1h, s1l, s2h, s2l;
  s1h = muldd (T1[i1][0], T1[i1][1], T2[i2][2], T2[i2][3], &s1l);
  s2h = muldd (T2[i2][0] , T2[i2][1], T1[i1][2], T1[i1][3], &s2l);
  double Sh, Sl;
  if (__builtin_expect (i1 == 101 && i2 >= 61, 0))
    Sh = fastsum (s2h, s2l, s1h, s1l, &Sl);
  else
    Sh = fastsum (s1h, s1l, s2h, s2l, &Sl);
  double c1h, c1l, c2h, c2l;
  c1h = muldd (T1[i1][2], T1[i1][3], T2[i2][2], T2[i2][3], &c1l);
  c2h = muldd (T1[i1][0], T1[i1][1], T2[i2][0], T2[i2][1], &c2l);
  double Ch, Cl;
  if(__builtin_expect(i1==50&&i2>64,0))
    Ch = fastsum (-c2h, -c2l, c1h, c1l, &Cl);
  else
    Ch = fastsum (c1h, c1l, -c2h, -c2l, &Cl);
  double rCh = r*Ch, rCl = Cl*r + __builtin_fma(Ch,r,-rCh);
  static const double cs[] = {0x1.1111111111111p-7, -0x1.a01a019d36453p-13};
  static const double cc[] = {0x1.5555555555555p-5, -0x1.6c16c168d68d7p-10};
  double rC3l, rC3h = muldd(rCh,rCl, 0x1.5555555555555p-2, 0x1.5555555555555p-56, &rC3l);
  double fl, fh;
  if(__builtin_expect(j==12868,0)){
    if(__builtin_fabs(Sh) > __builtin_fabs(rCh))
      fh = fastsum(Sh,Sl, rCh,rCl, &fl);
    else
      fh = fastsum(rCh,rCl, Sh,Sl, &fl);
    if(__builtin_fabs(Sh) > __builtin_fabs(rC3h))
      rC3h = fastsum(Sh,Sl, rC3h,rC3l, &rC3l);
    else
      rC3h = fastsum(rC3h,rC3l, Sh,Sl, &rC3l);
  } else {
    fh = fastsum(Sh,Sl, rCh,rCl, &fl);
    rC3h = fastsum(Sh,Sl, rC3h,rC3l, &rC3l);
  }
  rC3h *= 0.5;
  rC3l *= 0.5;

  double tt = r2h*(Sh*(cc[0]+r2h*cc[1]) + rCh*(cs[0]+r2h*cs[1])), e;

  rC3h = fasttwosum(rC3h, -tt, &e);
  rC3l += e;

  rC3h = muldd(rC3h,rC3l, r2h,r2l, &rC3l);

  fh = fastsum(fh,fl, -rC3h,-rC3l, &fl);
  fh = fasttwosum(fh, fl, &fl);

  b64u64_u  rl = {.f = fl};
  uint64_t d = (rl.u + 16)&(~(uint64_t)0>>12);
  if(__builtin_expect(d<=16, 0)){
    static const double wow[] = {
      0x1.2359262c76506p-13, 0x1.3f69df45a2f3bp-13, 0x1.55fd6fc227d9dp-12,
      0x1.96350d587e672p-12, 0x1.1e2e72ca9e866p-11, 0x1.3bc6ca12143b6p-11,
      0x1.bac75a647203p-11, 0x1.ce8b994974d0bp-11, 0x1.7f9ec1226e157p-10,
      0x1.b960ebdc4ec13p-10, 0x1.933fb67c4d0afp-9, 0x1.ab56a7ae04d3bp-9,
      0x1.b3f2ba40dbc66p-9, 0x1.e77a1b55ccf96p-9, 0x1.efe186fe553d9p-9,
      0x1.46a9ab5a2c58ep-8, 0x1.8e113d622c77ap-8, 0x1.dafa3c69f3426p-8,
      0x1.1c730dcd71abap-7, 0x1.5dc43f86236ccp-7, 0x1.e17faefac7797p-7,
      0x1.41db571d96126p-6, 0x1.9c412d62c144p-6, 0x1.275a3d78c01ecp-5,
      0x1.4cd45ddee2881p-5, 0x1.69949b3d51fb1p-5, 0x1.9283586503fep-5,
      0x1.d7bdcd778049fp-5, 0x1.fccdc252cad1fp-5, 0x1.0023629fc9899p-4,
      0x1.21857ad584f7fp-4, 0x1.2e36813a9874p-4, 0x1.456ac98461b72p-4,
      0x1.5231b416ba885p-4, 0x1.7f4ea0f3bbc6fp-4, 0x1.9c5c0d685abd4p-4,
      0x1.a202b3fb84788p-4, 0x1.45c341fb80643p-3, 0x1.6f9a0f284d491p-3,
      0x1.967cda38032b3p-3, 0x1.c49ac7cde7b4cp-3, 0x1.d5064e6fe82c5p-3,
      0x1.dd04b12d498c6p-3, 0x1.e3095cae52dd7p-3, 0x1.223c48cd64801p-2,
      0x1.50954b7bbf87bp-2, 0x1.6b30c65ac788ap-2, 0x1.bdc8830ddf4e6p-2,
      0x1.c881b16b684b1p-2, 0x1.e05b0e0a809bcp-2, 0x1.ed25c5eb8c916p-2,
      0x1.fe767739d0f6dp-2, 0x1.22f26f182fabdp-1, 0x1.3eb4df7c5532ap-1,
      0x1.41516c909749cp-1, 0x1.4566e96eb9313p-1, 0x1.540e24e5f33f3p-1,
      0x1.d98c4c612718dp-1, 0x1.ee539c9654a36p-1, 0x1.02c2f02bd16d5p+0,
      0x1.3aa301f6ebb1ep+0, 0x1.640ac66708cp+0, 0x1.8272d4fd7730bp+0,
      0x1.921fb54442d16p+0, 0x1.921fb54442d17p+0, 0x1.921fb54442d18p+0,
      0x1.921fb54442d19p+0, 0x1.921fb54442d1ap+0, 0x1.bbfa05708792dp+0,
      0x1.e2fae1619a6afp+0, 0x1.4dbe000d5c1d2p+1, 0x1.6756745770a51p+1,
      0x1.6e6198df13b76p+1, 0x1.920745cc24d5ep+1, 0x1.9255291b529ecp+1,
      0x1.aa0b46aa9cc59p+1};
    static const unsigned char tls[] = {
      7, 0, 7, 6, 4, 1, 4, 3, 0, 7, 4, 0, 6, 2, 1, 4, 2, 4, 5, 7, 4, 2, 5, 3,
      5, 0, 5, 5, 4, 5, 1, 6, 7, 4, 3, 4, 7, 6, 3, 6, 0, 2, 0, 3, 2, 6, 3, 3,
      1, 3, 2, 5, 4, 3, 5, 7, 0, 7, 3, 4, 2, 6, 6, 1, 1, 1, 1, 1, 0, 3, 4, 2,
      3, 6, 2, 6 };
    const uint64_t *db = (const uint64_t*)wow;
    int a = 0, b = sizeof(wow)/sizeof(wow[0]) - 1, m = (a + b)/2;
    while (a <= b) {
      if (db[m] < t.u){
	a = m + 1;
      } else if (__builtin_expect(db[m] == t.u, 0)) {
	b64u64_u jf = {.f = fh},
	  dr = {.u = ((jf.u&(0x7fful<<52)) - (54ul<<52))|((tls[m]&1ul)<<63)};
	uint64_t t0 = tls[m]>>1;
	for(int k = -1; k<=1; k++){
	  b64u64_u jk = {.u = jf.u + k};
	  if((jk.u&3) == t0){
	    fh = jk.f;
	    fl = dr.f;
	    break;
	  }
	}
	break;
      } else {
	b = m - 1;
      }
      m = (a + b)>>1;
    }
  }

  static const double Sgn[] = {1.0, -1.0};
  fh *= Sgn[sgn];
  fl *= Sgn[sgn];
  return fh + fl;
}
#endif

/* for each j, 0 <= j < 128, U1[j] contains sh, sl, ch, cl where
   sh+sl is a double-double approximation of sin(j*pi/2^7) and
   ch+cl is a double-double approximation of cos(j*pi/2^7).
   Generated by U1() from sin.sage */
static const double U1[128][4] = {
  {0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0},
  {0x1.92155f7a3667ep-6, -0x1.b1d63091a013p-64, 0x1.ffd886084cd0dp-1, -0x1.1354d4556e4cbp-55},
  {0x1.91f65f10dd814p-5, -0x1.912bd0d569a9p-61, 0x1.ff621e3796d7ep-1, -0x1.c57bc2e24aa15p-57},
  {0x1.2d52092ce19f6p-4, -0x1.9a088a8bf6b2cp-59, 0x1.fe9cdad01883ap-1, 0x1.521ecd0c67e35p-57},
  {0x1.917a6bc29b42cp-4, -0x1.e2718d26ed688p-60, 0x1.fd88da3d12526p-1, -0x1.87df6378811c7p-55},
  {0x1.f564e56a9730ep-4, 0x1.a2704729ae56dp-59, 0x1.fc26470e19fd3p-1, 0x1.1ec8668ecaceep-55},
  {0x1.2c8106e8e613ap-3, 0x1.13000a89a11ep-58, 0x1.fa7557f08a517p-1, -0x1.7a0a8ca13571fp-55},
  {0x1.5e214448b3fc6p-3, 0x1.531ff779ddac6p-57, 0x1.f8764fa714ba9p-1, 0x1.ab256778ffcb6p-56},
  {0x1.8f8b83c69a60bp-3, -0x1.26d19b9ff8d82p-57, 0x1.f6297cff75cbp-1, 0x1.562172a361fd3p-56},
  {0x1.c0b826a7e4f63p-3, -0x1.af1439e521935p-62, 0x1.f38f3ac64e589p-1, -0x1.d7bafb51f72e6p-56},
  {0x1.f19f97b215f1bp-3, -0x1.42deef11da2c4p-57, 0x1.f0a7efb9230d7p-1, 0x1.52c7adc6b4989p-56},
  {0x1.111d262b1f677p-2, 0x1.824c20ab7aa9ap-56, 0x1.ed740e7684963p-1, 0x1.e82c791f59cc2p-56},
  {0x1.294062ed59f06p-2, -0x1.5d28da2c4612dp-56, 0x1.e9f4156c62ddap-1, 0x1.760b1e2e3f81ep-55},
  {0x1.4135c94176601p-2, 0x1.0c97c4afa2518p-56, 0x1.e6288ec48e112p-1, -0x1.16b56f2847754p-57},
  {0x1.58f9a75ab1fddp-2, -0x1.efdc0d58cf62p-62, 0x1.e212104f686e5p-1, -0x1.014c76c126527p-55},
  {0x1.7088530fa459fp-2, -0x1.44b19e0864c5dp-56, 0x1.ddb13b6ccc23cp-1, 0x1.83c37c6107db3p-55},
  {0x1.87de2a6aea963p-2, -0x1.72cedd3d5a61p-57, 0x1.d906bcf328d46p-1, 0x1.457e610231ac2p-56},
  {0x1.9ef7943a8ed8ap-2, 0x1.6da81290bdbabp-57, 0x1.d4134d14dc93ap-1, -0x1.4ef5295d25af2p-55},
  {0x1.b5d1009e15ccp-2, 0x1.5b362cb974183p-57, 0x1.ced7af43cc773p-1, -0x1.e7b6bb5ab58aep-58},
  {0x1.cc66e9931c45ep-2, 0x1.6850e59c37f8fp-58, 0x1.c954b213411f5p-1, -0x1.2fb761e946603p-58},
  {0x1.e2b5d3806f63bp-2, 0x1.e0d891d3c6841p-58, 0x1.c38b2f180bdb1p-1, -0x1.6e0b1757c8d07p-56},
  {0x1.f8ba4dbf89abap-2, -0x1.2ec1fc1b776b8p-60, 0x1.bd7c0ac6f952ap-1, -0x1.825a732ac700ap-55},
  {0x1.073879922ffeep-1, -0x1.a5a014347406cp-55, 0x1.b728345196e3ep-1, -0x1.bc69f324e6d61p-55},
  {0x1.11eb3541b4b23p-1, -0x1.ef23b69abe4f1p-55, 0x1.b090a581502p-1, -0x1.926da300ffccep-55},
  {0x1.1c73b39ae68c8p-1, 0x1.b25dd267f66p-55, 0x1.a9b66290ea1a3p-1, 0x1.9f630e8b6dac8p-60},
  {0x1.26d054cdd12dfp-1, -0x1.5da743ef3770cp-55, 0x1.a29a7a0462782p-1, -0x1.128bb015df175p-56},
  {0x1.30ff7fce17035p-1, -0x1.efcc626f74a6fp-57, 0x1.9b3e047f38741p-1, -0x1.30ee286712474p-55},
  {0x1.3affa292050b9p-1, 0x1.e3e25e3954964p-56, 0x1.93a22499263fbp-1, 0x1.3d419a920df0bp-55},
  {0x1.44cf325091dd6p-1, 0x1.8076a2cfdc6b3p-57, 0x1.8bc806b151741p-1, -0x1.2c5e12ed1336dp-55},
  {0x1.4e6cabbe3e5e9p-1, 0x1.3c293edceb327p-57, 0x1.83b0e0bff976ep-1, -0x1.6f420f8ea3475p-56},
  {0x1.57d69348cecap-1, -0x1.75720992bfbb2p-55, 0x1.7b5df226aafafp-1, -0x1.0f537acdf0ad7p-56},
  {0x1.610b7551d2cdfp-1, -0x1.251b352ff2a37p-56, 0x1.72d0837efff96p-1, 0x1.0d4ef0f1d915cp-55},
  {0x1.6a09e667f3bcdp-1, -0x1.bdd3413b26456p-55, 0x1.6a09e667f3bcdp-1, -0x1.bdd3413b26456p-55},
  {0x1.72d0837efff96p-1, 0x1.0d4ef0f1d915cp-55, 0x1.610b7551d2cdfp-1, -0x1.251b352ff2a37p-56},
  {0x1.7b5df226aafafp-1, -0x1.0f537acdf0ad7p-56, 0x1.57d69348cecap-1, -0x1.75720992bfbb2p-55},
  {0x1.83b0e0bff976ep-1, -0x1.6f420f8ea3475p-56, 0x1.4e6cabbe3e5e9p-1, 0x1.3c293edceb327p-57},
  {0x1.8bc806b151741p-1, -0x1.2c5e12ed1336dp-55, 0x1.44cf325091dd6p-1, 0x1.8076a2cfdc6b3p-57},
  {0x1.93a22499263fbp-1, 0x1.3d419a920df0bp-55, 0x1.3affa292050b9p-1, 0x1.e3e25e3954964p-56},
  {0x1.9b3e047f38741p-1, -0x1.30ee286712474p-55, 0x1.30ff7fce17035p-1, -0x1.efcc626f74a6fp-57},
  {0x1.a29a7a0462782p-1, -0x1.128bb015df175p-56, 0x1.26d054cdd12dfp-1, -0x1.5da743ef3770cp-55},
  {0x1.a9b66290ea1a3p-1, 0x1.9f630e8b6dac8p-60, 0x1.1c73b39ae68c8p-1, 0x1.b25dd267f66p-55},
  {0x1.b090a581502p-1, -0x1.926da300ffccep-55, 0x1.11eb3541b4b23p-1, -0x1.ef23b69abe4f1p-55},
  {0x1.b728345196e3ep-1, -0x1.bc69f324e6d61p-55, 0x1.073879922ffeep-1, -0x1.a5a014347406cp-55},
  {0x1.bd7c0ac6f952ap-1, -0x1.825a732ac700ap-55, 0x1.f8ba4dbf89abap-2, -0x1.2ec1fc1b776b8p-60},
  {0x1.c38b2f180bdb1p-1, -0x1.6e0b1757c8d07p-56, 0x1.e2b5d3806f63bp-2, 0x1.e0d891d3c6841p-58},
  {0x1.c954b213411f5p-1, -0x1.2fb761e946603p-58, 0x1.cc66e9931c45ep-2, 0x1.6850e59c37f8fp-58},
  {0x1.ced7af43cc773p-1, -0x1.e7b6bb5ab58aep-58, 0x1.b5d1009e15ccp-2, 0x1.5b362cb974183p-57},
  {0x1.d4134d14dc93ap-1, -0x1.4ef5295d25af2p-55, 0x1.9ef7943a8ed8ap-2, 0x1.6da81290bdbabp-57},
  {0x1.d906bcf328d46p-1, 0x1.457e610231ac2p-56, 0x1.87de2a6aea963p-2, -0x1.72cedd3d5a61p-57},
  {0x1.ddb13b6ccc23cp-1, 0x1.83c37c6107db3p-55, 0x1.7088530fa459fp-2, -0x1.44b19e0864c5dp-56},
  {0x1.e212104f686e5p-1, -0x1.014c76c126527p-55, 0x1.58f9a75ab1fddp-2, -0x1.efdc0d58cf62p-62},
  {0x1.e6288ec48e112p-1, -0x1.16b56f2847754p-57, 0x1.4135c94176601p-2, 0x1.0c97c4afa2518p-56},
  {0x1.e9f4156c62ddap-1, 0x1.760b1e2e3f81ep-55, 0x1.294062ed59f06p-2, -0x1.5d28da2c4612dp-56},
  {0x1.ed740e7684963p-1, 0x1.e82c791f59cc2p-56, 0x1.111d262b1f677p-2, 0x1.824c20ab7aa9ap-56},
  {0x1.f0a7efb9230d7p-1, 0x1.52c7adc6b4989p-56, 0x1.f19f97b215f1bp-3, -0x1.42deef11da2c4p-57},
  {0x1.f38f3ac64e589p-1, -0x1.d7bafb51f72e6p-56, 0x1.c0b826a7e4f63p-3, -0x1.af1439e521935p-62},
  {0x1.f6297cff75cbp-1, 0x1.562172a361fd3p-56, 0x1.8f8b83c69a60bp-3, -0x1.26d19b9ff8d82p-57},
  {0x1.f8764fa714ba9p-1, 0x1.ab256778ffcb6p-56, 0x1.5e214448b3fc6p-3, 0x1.531ff779ddac6p-57},
  {0x1.fa7557f08a517p-1, -0x1.7a0a8ca13571fp-55, 0x1.2c8106e8e613ap-3, 0x1.13000a89a11ep-58},
  {0x1.fc26470e19fd3p-1, 0x1.1ec8668ecaceep-55, 0x1.f564e56a9730ep-4, 0x1.a2704729ae56dp-59},
  {0x1.fd88da3d12526p-1, -0x1.87df6378811c7p-55, 0x1.917a6bc29b42cp-4, -0x1.e2718d26ed688p-60},
  {0x1.fe9cdad01883ap-1, 0x1.521ecd0c67e35p-57, 0x1.2d52092ce19f6p-4, -0x1.9a088a8bf6b2cp-59},
  {0x1.ff621e3796d7ep-1, -0x1.c57bc2e24aa15p-57, 0x1.91f65f10dd814p-5, -0x1.912bd0d569a9p-61},
  {0x1.ffd886084cd0dp-1, -0x1.1354d4556e4cbp-55, 0x1.92155f7a3667ep-6, -0x1.b1d63091a013p-64},
  {0x1p+0, 0x0p+0, 0x0p+0, 0x0p+0},
  {0x1.ffd886084cd0dp-1, -0x1.1354d4556e4cbp-55, -0x1.92155f7a3667ep-6, 0x1.b1d63091a013p-64},
  {0x1.ff621e3796d7ep-1, -0x1.c57bc2e24aa15p-57, -0x1.91f65f10dd814p-5, 0x1.912bd0d569a9p-61},
  {0x1.fe9cdad01883ap-1, 0x1.521ecd0c67e35p-57, -0x1.2d52092ce19f6p-4, 0x1.9a088a8bf6b2cp-59},
  {0x1.fd88da3d12526p-1, -0x1.87df6378811c7p-55, -0x1.917a6bc29b42cp-4, 0x1.e2718d26ed688p-60},
  {0x1.fc26470e19fd3p-1, 0x1.1ec8668ecaceep-55, -0x1.f564e56a9730ep-4, -0x1.a2704729ae56dp-59},
  {0x1.fa7557f08a517p-1, -0x1.7a0a8ca13571fp-55, -0x1.2c8106e8e613ap-3, -0x1.13000a89a11ep-58},
  {0x1.f8764fa714ba9p-1, 0x1.ab256778ffcb6p-56, -0x1.5e214448b3fc6p-3, -0x1.531ff779ddac6p-57},
  {0x1.f6297cff75cbp-1, 0x1.562172a361fd3p-56, -0x1.8f8b83c69a60bp-3, 0x1.26d19b9ff8d82p-57},
  {0x1.f38f3ac64e589p-1, -0x1.d7bafb51f72e6p-56, -0x1.c0b826a7e4f63p-3, 0x1.af1439e521935p-62},
  {0x1.f0a7efb9230d7p-1, 0x1.52c7adc6b4989p-56, -0x1.f19f97b215f1bp-3, 0x1.42deef11da2c4p-57},
  {0x1.ed740e7684963p-1, 0x1.e82c791f59cc2p-56, -0x1.111d262b1f677p-2, -0x1.824c20ab7aa9ap-56},
  {0x1.e9f4156c62ddap-1, 0x1.760b1e2e3f81ep-55, -0x1.294062ed59f06p-2, 0x1.5d28da2c4612dp-56},
  {0x1.e6288ec48e112p-1, -0x1.16b56f2847754p-57, -0x1.4135c94176601p-2, -0x1.0c97c4afa2518p-56},
  {0x1.e212104f686e5p-1, -0x1.014c76c126527p-55, -0x1.58f9a75ab1fddp-2, 0x1.efdc0d58cf62p-62},
  {0x1.ddb13b6ccc23cp-1, 0x1.83c37c6107db3p-55, -0x1.7088530fa459fp-2, 0x1.44b19e0864c5dp-56},
  {0x1.d906bcf328d46p-1, 0x1.457e610231ac2p-56, -0x1.87de2a6aea963p-2, 0x1.72cedd3d5a61p-57},
  {0x1.d4134d14dc93ap-1, -0x1.4ef5295d25af2p-55, -0x1.9ef7943a8ed8ap-2, -0x1.6da81290bdbabp-57},
  {0x1.ced7af43cc773p-1, -0x1.e7b6bb5ab58aep-58, -0x1.b5d1009e15ccp-2, -0x1.5b362cb974183p-57},
  {0x1.c954b213411f5p-1, -0x1.2fb761e946603p-58, -0x1.cc66e9931c45ep-2, -0x1.6850e59c37f8fp-58},
  {0x1.c38b2f180bdb1p-1, -0x1.6e0b1757c8d07p-56, -0x1.e2b5d3806f63bp-2, -0x1.e0d891d3c6841p-58},
  {0x1.bd7c0ac6f952ap-1, -0x1.825a732ac700ap-55, -0x1.f8ba4dbf89abap-2, 0x1.2ec1fc1b776b8p-60},
  {0x1.b728345196e3ep-1, -0x1.bc69f324e6d61p-55, -0x1.073879922ffeep-1, 0x1.a5a014347406cp-55},
  {0x1.b090a581502p-1, -0x1.926da300ffccep-55, -0x1.11eb3541b4b23p-1, 0x1.ef23b69abe4f1p-55},
  {0x1.a9b66290ea1a3p-1, 0x1.9f630e8b6dac8p-60, -0x1.1c73b39ae68c8p-1, -0x1.b25dd267f66p-55},
  {0x1.a29a7a0462782p-1, -0x1.128bb015df175p-56, -0x1.26d054cdd12dfp-1, 0x1.5da743ef3770cp-55},
  {0x1.9b3e047f38741p-1, -0x1.30ee286712474p-55, -0x1.30ff7fce17035p-1, 0x1.efcc626f74a6fp-57},
  {0x1.93a22499263fbp-1, 0x1.3d419a920df0bp-55, -0x1.3affa292050b9p-1, -0x1.e3e25e3954964p-56},
  {0x1.8bc806b151741p-1, -0x1.2c5e12ed1336dp-55, -0x1.44cf325091dd6p-1, -0x1.8076a2cfdc6b3p-57},
  {0x1.83b0e0bff976ep-1, -0x1.6f420f8ea3475p-56, -0x1.4e6cabbe3e5e9p-1, -0x1.3c293edceb327p-57},
  {0x1.7b5df226aafafp-1, -0x1.0f537acdf0ad7p-56, -0x1.57d69348cecap-1, 0x1.75720992bfbb2p-55},
  {0x1.72d0837efff96p-1, 0x1.0d4ef0f1d915cp-55, -0x1.610b7551d2cdfp-1, 0x1.251b352ff2a37p-56},
  {0x1.6a09e667f3bcdp-1, -0x1.bdd3413b26456p-55, -0x1.6a09e667f3bcdp-1, 0x1.bdd3413b26456p-55},
  {0x1.610b7551d2cdfp-1, -0x1.251b352ff2a37p-56, -0x1.72d0837efff96p-1, -0x1.0d4ef0f1d915cp-55},
  {0x1.57d69348cecap-1, -0x1.75720992bfbb2p-55, -0x1.7b5df226aafafp-1, 0x1.0f537acdf0ad7p-56},
  {0x1.4e6cabbe3e5e9p-1, 0x1.3c293edceb327p-57, -0x1.83b0e0bff976ep-1, 0x1.6f420f8ea3475p-56},
  {0x1.44cf325091dd6p-1, 0x1.8076a2cfdc6b3p-57, -0x1.8bc806b151741p-1, 0x1.2c5e12ed1336dp-55},
  {0x1.3affa292050b9p-1, 0x1.e3e25e3954964p-56, -0x1.93a22499263fbp-1, -0x1.3d419a920df0bp-55},
  {0x1.30ff7fce17035p-1, -0x1.efcc626f74a6fp-57, -0x1.9b3e047f38741p-1, 0x1.30ee286712474p-55},
  {0x1.26d054cdd12dfp-1, -0x1.5da743ef3770cp-55, -0x1.a29a7a0462782p-1, 0x1.128bb015df175p-56},
  {0x1.1c73b39ae68c8p-1, 0x1.b25dd267f66p-55, -0x1.a9b66290ea1a3p-1, -0x1.9f630e8b6dac8p-60},
  {0x1.11eb3541b4b23p-1, -0x1.ef23b69abe4f1p-55, -0x1.b090a581502p-1, 0x1.926da300ffccep-55},
  {0x1.073879922ffeep-1, -0x1.a5a014347406cp-55, -0x1.b728345196e3ep-1, 0x1.bc69f324e6d61p-55},
  {0x1.f8ba4dbf89abap-2, -0x1.2ec1fc1b776b8p-60, -0x1.bd7c0ac6f952ap-1, 0x1.825a732ac700ap-55},
  {0x1.e2b5d3806f63bp-2, 0x1.e0d891d3c6841p-58, -0x1.c38b2f180bdb1p-1, 0x1.6e0b1757c8d07p-56},
  {0x1.cc66e9931c45ep-2, 0x1.6850e59c37f8fp-58, -0x1.c954b213411f5p-1, 0x1.2fb761e946603p-58},
  {0x1.b5d1009e15ccp-2, 0x1.5b362cb974183p-57, -0x1.ced7af43cc773p-1, 0x1.e7b6bb5ab58aep-58},
  {0x1.9ef7943a8ed8ap-2, 0x1.6da81290bdbabp-57, -0x1.d4134d14dc93ap-1, 0x1.4ef5295d25af2p-55},
  {0x1.87de2a6aea963p-2, -0x1.72cedd3d5a61p-57, -0x1.d906bcf328d46p-1, -0x1.457e610231ac2p-56},
  {0x1.7088530fa459fp-2, -0x1.44b19e0864c5dp-56, -0x1.ddb13b6ccc23cp-1, -0x1.83c37c6107db3p-55},
  {0x1.58f9a75ab1fddp-2, -0x1.efdc0d58cf62p-62, -0x1.e212104f686e5p-1, 0x1.014c76c126527p-55},
  {0x1.4135c94176601p-2, 0x1.0c97c4afa2518p-56, -0x1.e6288ec48e112p-1, 0x1.16b56f2847754p-57},
  {0x1.294062ed59f06p-2, -0x1.5d28da2c4612dp-56, -0x1.e9f4156c62ddap-1, -0x1.760b1e2e3f81ep-55},
  {0x1.111d262b1f677p-2, 0x1.824c20ab7aa9ap-56, -0x1.ed740e7684963p-1, -0x1.e82c791f59cc2p-56},
  {0x1.f19f97b215f1bp-3, -0x1.42deef11da2c4p-57, -0x1.f0a7efb9230d7p-1, -0x1.52c7adc6b4989p-56},
  {0x1.c0b826a7e4f63p-3, -0x1.af1439e521935p-62, -0x1.f38f3ac64e589p-1, 0x1.d7bafb51f72e6p-56},
  {0x1.8f8b83c69a60bp-3, -0x1.26d19b9ff8d82p-57, -0x1.f6297cff75cbp-1, -0x1.562172a361fd3p-56},
  {0x1.5e214448b3fc6p-3, 0x1.531ff779ddac6p-57, -0x1.f8764fa714ba9p-1, -0x1.ab256778ffcb6p-56},
  {0x1.2c8106e8e613ap-3, 0x1.13000a89a11ep-58, -0x1.fa7557f08a517p-1, 0x1.7a0a8ca13571fp-55},
  {0x1.f564e56a9730ep-4, 0x1.a2704729ae56dp-59, -0x1.fc26470e19fd3p-1, -0x1.1ec8668ecaceep-55},
  {0x1.917a6bc29b42cp-4, -0x1.e2718d26ed688p-60, -0x1.fd88da3d12526p-1, 0x1.87df6378811c7p-55},
  {0x1.2d52092ce19f6p-4, -0x1.9a088a8bf6b2cp-59, -0x1.fe9cdad01883ap-1, -0x1.521ecd0c67e35p-57},
  {0x1.91f65f10dd814p-5, -0x1.912bd0d569a9p-61, -0x1.ff621e3796d7ep-1, 0x1.c57bc2e24aa15p-57},
  {0x1.92155f7a3667ep-6, -0x1.b1d63091a013p-64, -0x1.ffd886084cd0dp-1, 0x1.1354d4556e4cbp-55},
};

/* for each j, 0 <= j < 128, U2[j] contains sh, sl, ch, cl where
   sh+sl is a double-double approximation of sin(j*pi/2^14) and
   ch+cl is a double-double approximation of cos(j*pi/2^14).
   Generated by U2() from sin.sage */
static const double U2[128][4] = {
  {0x0p+0, 0x0p+0, 0x1p+0, 0x0p+0},
  {0x1.921fb51aeb57cp-13, -0x1.a6e1d4916c435p-67, 0x1.ffffff621619cp-1, -0x1.7507dbbbd8fe6p-55},
  {0x1.921fb49ee4ea6p-12, 0x1.e894d744a453ep-66, 0x1.fffffd8858675p-1, -0x1.79f0e54748eabp-55},
  {0x1.2d97c6dc23a75p-11, -0x1.1c27727253849p-66, 0x1.fffffa72c6e9dp-1, 0x1.ffa0e2a0e5c29p-56},
  {0x1.921fb2aecb36p-11, 0x1.876157e566b4cp-65, 0x1.fffff62161a34p-1, -0x1.136dcb1f9b9c4p-57},
  {0x1.f6a79d8965ebap-11, -0x1.a2a951670a396p-65, 0x1.fffff09428963p-1, 0x1.2f505fba9c0bp-55},
  {0x1.2d97c396f8497p-10, -0x1.45cb4cc3d0fb1p-66, 0x1.ffffe9cb1bc62p-1, 0x1.0fe9fa3478ec9p-56},
  {0x1.5fdbb7af33fbbp-10, 0x1.c0c37a5d5101dp-65, 0x1.ffffe1c63b373p-1, 0x1.a50f155e87e56p-55},
  {0x1.921faaee6472ep-10, -0x1.ee52e284a9df8p-64, 0x1.ffffd88586ee6p-1, 0x1.1af64f173ae5bp-55},
  {0x1.c4639d358815ap-10, 0x1.0ffa0e86b4dccp-69, 0x1.ffffce08fef16p-1, 0x1.92e420a03bf59p-58},
  {0x1.f6a78e659d4b6p-10, 0x1.cb8a9b355ac8p-64, 0x1.ffffc250a346ap-1, 0x1.d099c02280c58p-56},
  {0x1.1475bf2fd13e2p-9, -0x1.4ae6ded9797ccp-63, 0x1.ffffb55c73f56p-1, 0x1.ee2ae1596963ap-55},
  {0x1.2d97b6824b087p-9, -0x1.9dae69dd4f97fp-63, 0x1.ffffa72c7105bp-1, -0x1.564ec4a452371p-55},
  {0x1.46b9ad1abb397p-9, -0x1.ac1b1950e9dc6p-64, 0x1.ffff97c09a803p-1, -0x1.961ec991858f1p-56},
  {0x1.5fdba2e9a1066p-9, 0x1.f1ed8abaaa608p-64, 0x1.ffff8718f06e7p-1, 0x1.793938ee3fc9ap-57},
  {0x1.78fd97df7ba51p-9, -0x1.2f2580c9ba856p-63, 0x1.ffff753572dacp-1, -0x1.5df4d689b3227p-57},
  {0x1.921f8becca4bap-9, 0x1.2ba407bcab5b2p-63, 0x1.ffff621621d02p-1, -0x1.6acfcebc82813p-56},
  {0x1.ab417f020c31p-9, 0x1.7fb6c9ae3781ap-65, 0x1.ffff4dbafd5a6p-1, -0x1.d017de84ebc3fp-55},
  {0x1.c463710fc08c9p-9, -0x1.f621bbef801cfp-64, 0x1.ffff38240586p-1, -0x1.f949383d834ep-61},
  {0x1.dd85620666965p-9, 0x1.5a1a4ecd8345ap-65, 0x1.ffff21513a606p-1, 0x1.ee48b39cc31cfp-56},
  {0x1.f6a751d67d871p-9, -0x1.66bdef183ff59p-63, 0x1.ffff09429bf7ap-1, -0x1.cabd3ded13a63p-55},
  {0x1.07e4a038424c1p-8, 0x1.7e6886ed9bec7p-63, 0x1.fffeeff82a5a7p-1, 0x1.6efc1cf23e273p-55},
  {0x1.147596e27d81fp-8, -0x1.bbfe9edba714ap-62, 0x1.fffed571e5989p-1, 0x1.151c0e15503cfp-55},
  {0x1.21068ce230028p-8, 0x1.6ce44ca79798dp-63, 0x1.fffeb9afcdc25p-1, 0x1.32c5af9456b9fp-57},
  {0x1.2d97822f996bcp-8, 0x1.3e5a15ed6aa3ep-62, 0x1.fffe9cb1e2e8dp-1, -0x1.10f44663fd601p-55},
  {0x1.3a2876c2f95cp-8, 0x1.03f5805a36f7ep-63, 0x1.fffe7e78251dep-1, 0x1.89b8834c8800cp-55},
  {0x1.46b96a948f72p-8, -0x1.3fbff884c87dap-65, 0x1.fffe5f0294744p-1, 0x1.5e809bdc855fcp-55},
  {0x1.534a5d9c9b4dp-8, -0x1.8653d3b3f3bf5p-62, 0x1.fffe3e5130ff5p-1, 0x1.8d94f3d6ca2d6p-57},
  {0x1.5fdb4fd35c8cbp-8, -0x1.4b73d73e9b437p-62, 0x1.fffe1c63fad33p-1, 0x1.429d08eb02c47p-55},
  {0x1.6c6c413112d15p-8, -0x1.c66913985b2a8p-63, 0x1.fffdf93af204ep-1, -0x1.4ade22d9f9483p-56},
  {0x1.78fd31adfdbbap-8, -0x1.03bf7bee2893dp-64, 0x1.fffdd4d616aap-1, -0x1.422710e8c8595p-55},
  {0x1.858e21425cecfp-8, -0x1.5af2befb31c3cp-63, 0x1.fffdaf3568d9p-1, 0x1.99c62dc3a22acp-58},
  {0x1.921f0fe670071p-8, 0x1.ab967fe6b7a9bp-64, 0x1.fffd8858e8a92p-1, 0x1.359c71883bcf7p-55},
  {0x1.9eaffd9276ac8p-8, -0x1.558e4e5ccafc7p-62, 0x1.fffd604096326p-1, -0x1.7d6bdd1c1f436p-60},
  {0x1.ab40ea3eb0803p-8, -0x1.7fe2a8e15f257p-67, 0x1.fffd36ec718d7p-1, -0x1.5c59ffcfba7e4p-56},
  {0x1.b7d1d5e35d25dp-8, 0x1.9cfe67f192cd7p-62, 0x1.fffd0c5c7ad3dp-1, -0x1.19d99906d21e9p-55},
  {0x1.c462c078bc41bp-8, 0x1.9f6db9a6ada3cp-62, 0x1.fffce090b21fcp-1, -0x1.0389138f7e5efp-55},
  {0x1.d0f3a9f70d78cp-8, -0x1.a08e48b732df4p-65, 0x1.fffcb389178c4p-1, 0x1.2883b1fe7d578p-56},
  {0x1.dd84925690709p-8, -0x1.0441c939bd6dap-62, 0x1.fffc8545ab352p-1, 0x1.61157c68f984bp-55},
  {0x1.ea15798f84cf7p-8, -0x1.bd41d0fe0c595p-62, 0x1.fffc55c66d36fp-1, -0x1.ab98d6fcf9cd6p-58},
  {0x1.f6a65f9a2a3c6p-8, -0x1.de8c48783f3aep-62, 0x1.fffc250b5daefp-1, -0x1.13b48657081cdp-55},
  {0x1.019ba237602f9p-7, -0x1.2479c7145927dp-61, 0x1.fffbf3147cbb3p-1, -0x1.6a078c1216e21p-55},
  {0x1.07e41402c3701p-7, -0x1.93e58a2a04d27p-67, 0x1.fffbbfe1ca7a8p-1, -0x1.698b82e41cfaep-56},
  {0x1.0e2c852b5eb46p-7, -0x1.9b8495e4ec41dp-61, 0x1.fffb8b73470c8p-1, -0x1.bc4e4be396442p-55},
  {0x1.1474f5ad51d17p-7, -0x1.a0fb850f9330dp-62, 0x1.fffb55c8f2917p-1, 0x1.6af540576eb51p-55},
  {0x1.1abd6584bc9cbp-7, 0x1.21b927d9670ccp-61, 0x1.fffb1ee2cd2a9p-1, -0x1.3ee83959a5f1dp-56},
  {0x1.2105d4adbeecp-7, 0x1.4e6ed5742dcf7p-61, 0x1.fffae6c0d6f9ap-1, -0x1.102fd2cb5b72p-56},
  {0x1.274e43247895ap-7, -0x1.46ca3539434eap-63, 0x1.fffaad6310214p-1, 0x1.b8629443a759ap-55},
  {0x1.2d96b0e509703p-7, -0x1.1e9131ff52dc9p-63, 0x1.fffa72c978c4fp-1, -0x1.22cb000328f91p-55},
  {0x1.33df1deb9152dp-7, 0x1.f5fb935dbf892p-62, 0x1.fffa36f41108bp-1, 0x1.55d1c829905c6p-57},
  {0x1.3a278a3430152p-7, -0x1.2c26d82855518p-63, 0x1.fff9f9e2d9118p-1, 0x1.102f90f3495ep-57},
  {0x1.406ff5bb058f1p-7, 0x1.6e66ff4bcd327p-61, 0x1.fff9bb95d105p-1, 0x1.7fb680db05f83p-55},
  {0x1.46b8607c31993p-7, 0x1.91d1ac2a2ced6p-61, 0x1.fff97c0cf909bp-1, -0x1.b418b1f88cf02p-57},
  {0x1.4d00ca73d40c8p-7, -0x1.fae7469c8dbd1p-61, 0x1.fff93b485146bp-1, -0x1.43e02406fb373p-55},
  {0x1.5349339e0cc25p-7, -0x1.ec016ea8e8f27p-64, 0x1.fff8f947d9e3fp-1, -0x1.2cefd02ebe75cp-59},
  {0x1.59919bf6fb94bp-7, -0x1.1ca7261e51e18p-62, 0x1.fff8b60b930a3p-1, 0x1.97dbe1fdcab9p-56},
  {0x1.5fda037ac05e1p-7, -0x1.ff59bf4b574eep-61, 0x1.fff871937ce2fp-1, -0x1.38dae49f0be32p-57},
  {0x1.66226a257af95p-7, 0x1.3cd7a66dfafc8p-61, 0x1.fff82bdf97986p-1, -0x1.507b1e2913e3bp-57},
  {0x1.6c6acff34b421p-7, 0x1.b03798bbe7197p-63, 0x1.fff7e4efe3558p-1, 0x1.ebe425d873a16p-57},
  {0x1.72b334e051144p-7, 0x1.caf2df425ef35p-62, 0x1.fff79cc460462p-1, -0x1.6c8d94412b92ap-55},
  {0x1.78fb98e8ac4c8p-7, -0x1.555202816f3a5p-61, 0x1.fff7535d0e96bp-1, -0x1.c5727ff41299bp-56},
  {0x1.7f43fc087cc7dp-7, 0x1.00a2a3a85d253p-61, 0x1.fff708b9ee748p-1, -0x1.7ff8f5cefe85ap-60},
  {0x1.858c5e3be264p-7, -0x1.ee9c7105192d8p-62, 0x1.fff6bcdb000dap-1, -0x1.65348e4f5cd16p-57},
  {0x1.8bd4bf7efcff3p-7, 0x1.66640846e6109p-63, 0x1.fff66fc04390dp-1, 0x1.77f5704f2756cp-55},
  {0x1.921d1fcdec784p-7, 0x1.9878ebe836d9dp-61, 0x1.fff62169b92dbp-1, 0x1.5dda3c81fbd0dp-55},
  {0x1.98657f24d0aeap-7, 0x1.11d4d1d1806d2p-61, 0x1.fff5d1d761149p-1, 0x1.e7b86f31e2875p-63},
  {0x1.9eaddd7fc9825p-7, -0x1.5adb5adfd3669p-62, 0x1.fff581093b768p-1, -0x1.3eeee9513ae4cp-55},
  {0x1.a4f63adaf6d3ep-7, -0x1.f838d8dd726dcp-64, 0x1.fff52eff48855p-1, -0x1.514532e82de3fp-57},
  {0x1.ab3e973278849p-7, 0x1.35b145a18353cp-61, 0x1.fff4dbb98873ap-1, 0x1.857663dc92252p-55},
  {0x1.b186f2826e765p-7, -0x1.6b3b32921e88cp-61, 0x1.fff48737fb74ep-1, -0x1.d462b0c1681e7p-58},
  {0x1.b7cf4cc6f88b8p-7, -0x1.32c6a4623533p-62, 0x1.fff4317aa1bd2p-1, -0x1.6a0039355b637p-55},
  {0x1.be17a5fc36a75p-7, 0x1.c404548edfc51p-63, 0x1.fff3da817b814p-1, -0x1.2b464e14730cp-55},
  {0x1.c45ffe1e48ad9p-7, 0x1.4060e4bd32e79p-63, 0x1.fff3824c88f6fp-1, -0x1.ed820f3fe698p-55},
  {0x1.caa855294e82bp-7, 0x1.eca8fa79834dep-66, 0x1.fff328dbca549p-1, -0x1.6ca1af321dc1cp-55},
  {0x1.d0f0ab19680bdp-7, -0x1.e9c612f8a4102p-64, 0x1.fff2ce2f3fd15p-1, -0x1.62ccf8bdb122fp-56},
  {0x1.d738ffeab52ecp-7, -0x1.c110bc5257c8cp-62, 0x1.fff27246e9a52p-1, -0x1.16d6346f3f55fp-59},
  {0x1.dd81539955d2p-7, -0x1.eaf6d880c7f01p-61, 0x1.fff21522c808bp-1, 0x1.a190e2d2eaf6ep-56},
  {0x1.e3c9a62169dcbp-7, 0x1.6666333e35e32p-61, 0x1.fff1b6c2db358p-1, -0x1.f4664a9ad833cp-56},
  {0x1.ea11f77f1136ep-7, -0x1.7a617506aacb8p-62, 0x1.fff157272365bp-1, 0x1.44ef64386476p-57},
  {0x1.f05a47ae6bc91p-7, -0x1.0935bb936be66p-64, 0x1.fff0f64fa0d45p-1, -0x1.ac9acb9bb8d2dp-56},
  {0x1.f6a296ab997cbp-7, -0x1.f2943d8fe7033p-61, 0x1.fff0943c53bd1p-1, -0x1.47399f361d158p-55},
  {0x1.fceae472ba3bcp-7, -0x1.3a08b1804396p-64, 0x1.fff030ed3c5c7p-1, -0x1.2504b1dc047e8p-55},
  {0x1.0199987ff6f89p-6, 0x1.8daccef2e2556p-60, 0x1.ffefcc625aefbp-1, 0x1.fbb0dd6af1603p-59},
  {0x1.04bdbe27aa444p-6, 0x1.0c87a06c6956dp-60, 0x1.ffef669bafb4ep-1, -0x1.b5a4e5318177bp-58},
  {0x1.07e1e32e86f72p-6, -0x1.f40145dde0463p-60, 0x1.ffeeff993aeacp-1, -0x1.99db334be163bp-58},
  {0x1.0b0607929d07ap-6, 0x1.f989b780d54e7p-61, 0x1.ffee975afcd0ep-1, -0x1.2df82b9320ecdp-55},
  {0x1.0e2a2b51fc6cep-6, -0x1.2ca118aaa212fp-61, 0x1.ffee2de0f5a78p-1, 0x1.9daa217acc1bp-58},
  {0x1.114e4e6ab51e2p-6, 0x1.10d53d5e0ffc5p-62, 0x1.ffedc32b25afcp-1, -0x1.bb89d7b0d4487p-64},
  {0x1.147270dad7133p-6, 0x1.8769e00e018p-63, 0x1.ffed57398d2b7p-1, -0x1.0bb0db6384cf3p-55},
  {0x1.179692a072443p-6, 0x1.d40c977f2fe27p-60, 0x1.ffecea0c2c5d2p-1, -0x1.7a22a10e4cb4fp-55},
  {0x1.1abab3b996a9dp-6, -0x1.e21a2ef2391d3p-60, 0x1.ffec7ba303882p-1, 0x1.b209ec7100e1ap-56},
  {0x1.1dded424543cep-6, 0x1.2ecdb3d03a70bp-61, 0x1.ffec0bfe12f0ap-1, 0x1.8b1483b4090edp-56},
  {0x1.2102f3debaf6fp-6, -0x1.45a1ea37e10eap-62, 0x1.ffeb9b1d5adb7p-1, 0x1.d5b363c2f437p-55},
  {0x1.242712e6dad1cp-6, 0x1.0317d8ce1d1b9p-60, 0x1.ffeb2900db8e4p-1, 0x1.18f6a66301d6p-57},
  {0x1.274b313ac3c7ap-6, 0x1.b29e49acf7e7ep-60, 0x1.ffeab5a8954f6p-1, 0x1.0652cfe4cebfap-55},
  {0x1.2a6f4ed885d35p-6, -0x1.5a8a13d8f9966p-60, 0x1.ffea41148866p-1, 0x1.b77e2735c69f9p-55},
  {0x1.2d936bbe30efdp-6, 0x1.b5f91ee371d64p-61, 0x1.ffe9cb44b51a1p-1, 0x1.5b43366df667p-56},
  {0x1.30b787e9d518dp-6, 0x1.b620c940df42ap-60, 0x1.ffe954391bb43p-1, 0x1.ddd5377beef1ap-56},
  {0x1.33dba359824a6p-6, -0x1.6c4020b5f32a8p-61, 0x1.ffe8dbf1bc7dep-1, -0x1.cc391831cbcdp-55},
  {0x1.36ffbe0b4880ep-6, 0x1.4b88026f4e07bp-61, 0x1.ffe8626e97c13p-1, 0x1.cdb27525fcaf7p-56},
  {0x1.3a23d7fd37b96p-6, -0x1.dd9e849a0a5d5p-60, 0x1.ffe7e7afadc94p-1, -0x1.db1225b06b02cp-55},
  {0x1.3d47f12d5ff12p-6, 0x1.add51eab25dd4p-60, 0x1.ffe76bb4fee1ap-1, -0x1.0c43dae9e89c6p-57},
  {0x1.406c0999d1263p-6, -0x1.cf7cbeea1ddcfp-61, 0x1.ffe6ee7e8b56ep-1, 0x1.7ecac48a7b7bcp-58},
  {0x1.439021409b56cp-6, 0x1.81352c74c7ae2p-62, 0x1.ffe6700c53764p-1, -0x1.4f502302afd49p-55},
  {0x1.46b4381fce81bp-6, 0x1.ff9f89fb65be3p-60, 0x1.ffe5f05e578dbp-1, -0x1.b71f7469ecc13p-56},
  {0x1.49d84e357aa66p-6, -0x1.9143373e36698p-61, 0x1.ffe56f7497ecp-1, -0x1.ddf92a8d40b3p-55},
  {0x1.4cfc637fafc48p-6, -0x1.a3541bb15d5bap-60, 0x1.ffe4ed4f14e0ap-1, 0x1.e48478b2f0ba1p-56},
  {0x1.502077fc7ddc5p-6, 0x1.25236a8347d96p-60, 0x1.ffe469edcebbfp-1, 0x1.8fc7557d6bc53p-55},
  {0x1.53448ba9f4eebp-6, 0x1.e6adde30ec695p-61, 0x1.ffe3e550c5cfp-1, -0x1.58fb3ae72192ep-55},
  {0x1.56689e8624fcep-6, -0x1.b967a2d3ee5b3p-61, 0x1.ffe35f77fa6b8p-1, -0x1.a7661b121429ap-57},
  {0x1.598cb08f1e089p-6, 0x1.56f9c34c216b4p-60, 0x1.ffe2d8636ce41p-1, 0x1.b4243168ab5a4p-57},
  {0x1.5cb0c1c2f0142p-6, 0x1.c2c42ac4d436cp-60, 0x1.ffe250131d8cp-1, 0x1.ed2d09d5789d3p-55},
  {0x1.5fd4d21fab226p-6, -0x1.0c0a91c37851cp-61, 0x1.ffe1c6870cb77p-1, 0x1.89aa14768323ep-55},
  {0x1.62f8e1a35f369p-6, -0x1.44ab93e20133fp-60, 0x1.ffe13bbf3abb3p-1, 0x1.6721eee289b1dp-55},
  {0x1.661cf04c1c548p-6, 0x1.3d6cc7ce6ff57p-60, 0x1.ffe0afbba7ecep-1, 0x1.6ecd981464044p-57},
  {0x1.6940fe17f280bp-6, -0x1.84bf8988cc42p-60, 0x1.ffe0227c54a2dp-1, 0x1.eeed1d70e13e2p-55},
  {0x1.6c650b04f1bfep-6, 0x1.53e382471fa0cp-69, 0x1.ffdf940141344p-1, -0x1.a63fc248c2dacp-55},
  {0x1.6f8917112a179p-6, 0x1.b6de0ad6e9445p-60, 0x1.ffdf044a6df8fp-1, -0x1.798cfb498d8b1p-55},
  {0x1.72ad223aab8ddp-6, -0x1.cf17a0ad13b4ap-60, 0x1.ffde7357db499p-1, 0x1.09337deb3ff5ap-59},
  {0x1.75d12c7f8629p-6, -0x1.07d8473495aaep-60, 0x1.ffdde129897f9p-1, 0x1.44e8786dedb6bp-55},
  {0x1.78f535ddc9f04p-6, 0x1.b194ad9b1aa97p-61, 0x1.ffdd4dbf78f52p-1, 0x1.216679a26323dp-55},
  {0x1.7c193e5386eb4p-6, 0x1.725118c2860e7p-63, 0x1.ffdcb919aa053p-1, -0x1.55c1ce37ae6f6p-56},
  {0x1.7f3d45decd222p-6, 0x1.2fe8cfe80809dp-60, 0x1.ffdc23381d0b6p-1, 0x1.00d15f9603fedp-57},
  {0x1.82614c7dac9dbp-6, 0x1.e95f011ac18c6p-63, 0x1.ffdb8c1ad2643p-1, 0x1.e79471744a0cbp-56},
  {0x1.8585522e35674p-6, -0x1.d5a766c5894c1p-60, 0x1.ffdaf3c1ca6cep-1, -0x1.9c0f0834fa422p-56},
  {0x1.88a956ee7788ap-6, 0x1.19dfbd12ff726p-61, 0x1.ffda5a2d05835p-1, 0x1.6dbc405ff8e25p-55},
  {0x1.8bcd5abc830c7p-6, -0x1.f7d609f14a311p-60, 0x1.ffd9bf5c84066p-1, -0x1.35167a6ce828dp-55},
  {0x1.8ef15d9667fdap-6, -0x1.668dc9ba82a7p-60, 0x1.ffd9235046557p-1, -0x1.c3c6e16616fb2p-56},
};

static double
moderate_exceptions (double x, double y)
{
  static double e[][3] = {
  };
  double ax = __builtin_fabs (x);
  for (unsigned int i = 0; i < sizeof(e)/(3*sizeof(double)); i++) {
    if (ax == e[i][0])
      return (x == ax) ? e[i][1] + e[i][2] : - e[i][1] - e[i][2];
  }
  return y;
}

// accurate path for 0x1.95e4p+1 <= |x| < 2^31
static double
sin_accurate_moderate (double x, double ax)
{
  int sbit = (x > 0) ? 0 : 1;
  static const double invpi = 0x1.45f306dc9c883p+12;
  double k = __builtin_roundeven (invpi * ax);
  static const double pih = -0x1.921fb54442d18p-13,
    pil = -0x1.1a62633145c07p-67, pis = 0x1.f1976b7ed8fbcp-123;
  // |pih+pil+pis+pi/2^14| < 2^-176.620
  double rh = __builtin_fma (k, pih, ax), rl = k * pil,
    rs = k * pis + __builtin_fma (k, pil, -rl);
  rh = fasttwosum (rh, rl, &rl);
  rl += rs;

  uint64_t j = k;
  sbit = sbit ^ ((j >> 14) & 1); // reduction by an odd multiple of pi?
  int i1 = (j >> 7) & 0x7f, i2 = j & 0x7f;
  double s1h, s1l, s2h, s2l;
  // s1h approximates sin(t1)*cos(t2)
  s1h = muldd_acc (U1[i1][0], U1[i1][1], U2[i2][2], U2[i2][3], &s1l);
  // s2h approximates cos(t1)*sin(t2)
  s2h = muldd_acc (U2[i2][0], U2[i2][1], U1[i1][2], U1[i1][3], &s2l);
  double Sh, Sl;
  // Sh approximates sin(t1+t2)
  Sh = fastsum (s1h, s1l, s2h, s2l, &Sl);
  double c1h, c1l, c2h, c2l;
  // c1h+c1l approximates cos(t1)*cos(t2)
  c1h = muldd_acc (U1[i1][2], U1[i1][3], U2[i2][2], U2[i2][3], &c1l);
  // c2h+c2l approximates sin(t1)*sin(t2)
  c2h = muldd_acc (U1[i1][0], U1[i1][1], U2[i2][0], U2[i2][1], &c2l);
  double Ch, Cl;
  // Ch approximates cos(t1+t2)
  Ch = fastsum (c1h, c1l, -c2h, -c2l, &Cl);

  double r2h = rh * rh, r2l = __builtin_fma (rh, rh, -r2h) + 2 * rh * rl;
  /* for |r| <= 2^-13.339, the polynomial
     ps[0]*r+(ps[1]+ps[2])*r^3+ps[3]*r^5+ps[4]*r^7
     approximates sin(r) with relative error < 2^-121.292, and the polynomial
     pc[0]*r^2+pc[1]*r^4+pc[2]*r^6 approximates cos(r)-1 with
     absolute error < 2^-114.946 (cf sinmoderate_acc.sollya).
     The reason why we use the relative error for sin(r) is that for x tiny,
     we have i1=i2=0, thus the approximation is simply sh+sl, and we need
     a small relative error. */
  static double ps[] = {1.0, -0x1.5555555555555p-3, -0x1.555555551de06p-57,
                        0x1.1111111111111p-7, -0x1.a01a006eb9947p-13};
  static double pc[] = {-0x1p-1, 0x1.5555555555555p-5, -0x1.6c16bcc416d44p-10};
  double sh, sl, t;
  // since |ps[4]*r^7/ps[0]*r| < 2^-92, we can compute ps[4]*r^2 as double
  sh = r2h * ps[4];
  // since |ps[3]*r^5/ps[0]*r| < 2^-60, we can still compute in binary64
  // add ps[3]
  sh += ps[3];
  sh = r2h * sh;
  // sh approximates ps[3]*r^2+ps[4]*r^4
  // add ps[1]+ps[2]
  sh = fasttwosum (ps[1], sh, &sl);
  sl += ps[2];
  sh = muldd (r2h, r2l, sh, sl, &sl);
  // sh+sl approximates (ps[1]+ps[2])*r^2+ps[3]*r^4
  // add ps[0]
  sh = fasttwosum (ps[0], sh, &t);
  sl += t;
  // multiply by rh+rl
  sh = muldd (rh, rl, sh, sl, &sl);
  // now sh+sl approximates sin(rh+rl)
  double ch, cl;
  // since |pc[2]*r^6| < 2^-89, we can compute pc[2]*r^2 as double
  ch = pc[2] * r2h;
  // since |pc[1]*r^4| < 2^-57, we can still compute in binary64
  // add pc[1]
  ch += pc[1];
  ch = r2h * ch;
  // ch approximates pc[1]*r^2+pc[2]*r^4
  // add pc[0]
  ch = fasttwosum (pc[0], ch, &cl);
  ch = muldd (r2h, r2l, ch, cl, &cl);
  // add 1
  ch = fasttwosum (1.0, ch, &t);
  cl += t;
  // now ch+cl approximates cos(rh+rl)

  // we now have to compute (Sh+Sl)*(ch+cl) + (Ch+Cl)*(sh+sl)
  Sh = muldd_acc (Sh, Sl, ch, cl, &Sl);
  Ch = muldd_acc (Ch, Cl, sh, sl, &Cl);
  // twosum is faster than a test and fasttwosum by about 1 cycle
  sh = twosum (Sh, Ch, &sl);
  sl += Sl + Cl;
  double ret = (sbit == 0) ? sh + sl : - sh - sl;
  // check worst cases
  b64u64_u z = {.f = sl};
  // if (x == 0x1.920745cc24d5ep+1) printf ("p+ z.u=%lx\n", z.u);
  if (((z.u + 5) & 0x7fffffffffffull) <= 6)
    return moderate_exceptions (x, ret);
  return ret;
}

// fast path for 0x1.95e4p+1 <= |x| < 2^31
// ax = |x| and eps is the error bound for the rounding test
// see proof of correctness in sin.pdf
static double
cr_sin_moderate (double x, double ax)
{
  int sbit = (x > 0) ? 0 : 1;
  static const double invpi = 0x1.45f306dc9c883p+12;
  // |invpi/2^14 - 1/pi| < 2^-55.496
  double k = __builtin_roundeven (invpi * ax);
  static const double pih = -0x1.921fb54442d18p-13,
    pil = -0x1.1a62633145c07p-67;
  // |2^14*(pih + pil) + pi| < 2^-108.041
  double rh = __builtin_fma (k, pih, ax), rl = k * pil; // rh is exact

  double r = rh + rl; // |r| < 2^-13.339 (see sin.pdf)
  double r2 = r * r;
  uint64_t j = k;
  sbit = sbit ^ ((j >> 14) & 1); // reduction by an odd multiple of pi?
  int i1 = (j >> 7) & 0x7f, i2 = j & 0x7f;
  double s1h, s1l, s2h, s2l;
  s1h = muldd (U1[i1][0], U1[i1][1], U2[i2][2], U2[i2][3], &s1l);
  s2h = muldd (U2[i2][0], U2[i2][1], U1[i1][2], U1[i1][3], &s2l);
  double Sh, Sl;
  Sh = fastsum (s1h, s1l, s2h, s2l, &Sl);
  double Ch = U1[i1][2] * U2[i2][2] - U1[i1][0] * U2[i2][0];

  /* for |r| <= 2^-13.339, the polynomial r - 0x1.55555553068fp-3 * r^3
     approximates sin(r) with absolute error < 2^-76.494, and the polynomial
     -0.5 * r^2 + 0x1.55555553bfd3p-5 * r^4 approximates cos(r)-1 with
     absolute error < 2^-92.723 (cf sinmoderate.sollya) */
  double sh = r * (1.0 - 0x1.55555553068fp-3 * r2);
  double ch = r2 * (-0.5 + 0x1.55555553bfd3p-5 * r2);
  double fh = Sh, fl = Sl + Sh*ch + Ch*sh;
  static double Sgn[] = {1.0, -1.0};
  fh = Sgn[sbit] * fh;
  fl = Sgn[sbit] * fl;
  static double eps = 0x1.dep-64;
  double lb = fh + (fl - eps), ub = fh + (fl + eps);
  if (__builtin_expect (lb == ub, 1)) return lb;
  return sin_accurate_moderate (x, ax);
  // return sin_accurate (x);
}

// fast path for |x| >= 2^31
// ax = |x| and eps is the error bound for the rounding test
static double
cr_sin_large (double x, double ax)
{
  double r;
  uint64_t j = reduce_large (&r, ax);
  // now x/(2pi) ~ k + j/2^15 + r with 0 <= r < 2^-15
  int sbit = (x > 0) ? 0 : 1;

  double r2 = r * r;
  sbit = sbit ^ ((j >> 14) & 1); // reduction by an odd multiple of pi?
  int i1 = (j >> 7) & 0x7f, i2 = j & 0x7f;
  double s1h, s1l, s2h, s2l;
  s1h = muldd (U1[i1][0], U1[i1][1], U2[i2][2], U2[i2][3], &s1l);
  s2h = muldd (U2[i2][0], U2[i2][1], U1[i1][2], U1[i1][3], &s2l);
  double Sh, Sl;
  Sh = fastsum (s1h, s1l, s2h, s2l, &Sl);
  double Ch = U1[i1][2] * U2[i2][2] - U1[i1][0] * U2[i2][0];

  double sh = r * (0x1.921fb54442d18p2 - 0x1.4abbcdb6b26d1p5 * r2);
  double ch = r2 * (-0x1.3bd3cc9be45dep4 + 0x1.03c1eee483083p6 * r2);
  double fh = Sh, fl = Sl + Sh*ch + Ch*sh;
  static double Sgn[] = {1.0, -1.0};
  fh = Sgn[sbit] * fh;
  fl = Sgn[sbit] * fl;
  static double eps = 0x1.41p-63;
  double lb = fh + (fl - eps), ub = fh + (fl + eps);
  if (__builtin_expect (lb == ub, 1)) return lb;
  return sin_accurate (x);
}

double
cr_sin (double x)
{
  b64u64_u t = {.f = x};
  double ax = __builtin_fabs (x);
  int e = (t.u >> 52) & 0x7ff;

  // deal with tiny x to avoid underflow
  if (__builtin_expect((t.u<<1) <= 0x7cae26e892247decull, 0)) {
    if (x == 0)
      return x;
    // Taylor expansion of sin(x) is x - x^3/6 around zero
    // for x=-0, fma (x, -0x1p-54, x) returns +0
    /* We have underflow when 0 < |x| < 2^-1022 or when |x| = 2^-1022
       and rounding towards zero. */
    double res = __builtin_fma (x, -0x1p-54, x);
#ifdef CORE_MATH_SUPPORT_ERRNO
    if (ax < 0x1p-1022 || __builtin_fabs (res) < 0x1p-1022)
      errno = ERANGE; // underflow
#endif
    return res;
  }

  if (e < 1054) return cr_sin_moderate (x, ax); // |x| < 2^31

  if (__builtin_expect (e == 0x7ff, 0)) /* NaN, +Inf and -Inf. */
    {
      if ((t.u << 1) == 0x7ffull<<53){ // +/-Inf
#ifdef CORE_MATH_SUPPORT_ERRNO
        errno = EDOM;
#endif
        return x - x; // raises invalid
      }
      return x + x;
    }

  // now |x| >= 2^31
  return cr_sin_large (x, ax);
}
