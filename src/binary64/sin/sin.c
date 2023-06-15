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

#define TRACE -0x1.005023d32fee5p+1

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

/* This table approximates 1/(2pi) downwards with precision 1216:
   1/(2*pi) ~ T[0]/2^64 + T[1]/2^128 + ... + T[i]/2^((i+1)*64) + ...
   Computed with computeT() from sin.sage. */
static const uint64_t T[19] = {
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

/* This table contains 128-bit approximations of sin(i/2^8) for 0 <= i < 256
   (to nearest).
   Each entry is to be interpreted as (hi/2^64+lo/2^128)*2^ex*(-1)*sgn.
   Generated with computeS() from sin.sage. */
static const dint64_t S[256] = {
   {.hi = 0x0, .lo = 0x0, .ex = 128, .sgn=0},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=0},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=0},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=0},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=0},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=0},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=0},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=0},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=0},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=0},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=0},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=0},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=0},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=0},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=0},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=0},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=0},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=0},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=0},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=0},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=0},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=0},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=0},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=0},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=0},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=0},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=0},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=0},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=0},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=0},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=0},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=0},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=0},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=0},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=0},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=0},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=0},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=0},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=0},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=0},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=0},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=0},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=0},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=0},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=0},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=0},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=0},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=0},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=0},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=0},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=0},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=0},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=0},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=0},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=0},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=0},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=0},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=0},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=0},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=0},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=0},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=0},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=0},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=0},
   {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=0},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=0},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=0},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=0},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=0},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=0},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=0},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=0},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=0},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=0},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=0},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=0},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=0},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=0},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=0},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=0},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=0},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=0},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=0},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=0},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=0},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=0},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=0},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=0},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=0},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=0},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=0},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=0},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=0},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=0},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=0},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=0},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=0},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=0},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=0},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=0},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=0},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=0},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=0},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=0},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=0},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=0},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=0},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=0},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=0},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=0},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=0},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=0},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=0},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=0},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=0},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=0},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=0},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=0},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=0},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=0},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=0},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=0},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=0},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=0},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=0},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=0},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=0},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=0},
   {.hi = 0x0, .lo = 0x0, .ex = 128, .sgn=0},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=1},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=1},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=1},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=1},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=1},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=1},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=1},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=1},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=1},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=1},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=1},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=1},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=1},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=1},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=1},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=1},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=1},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=1},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=1},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=1},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=1},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=1},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=1},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=1},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=1},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=1},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=1},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=1},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=1},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=1},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=1},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=1},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=1},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=1},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=1},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=1},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=1},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=1},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=1},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=1},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=1},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=1},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=1},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=1},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=1},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=1},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=1},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=1},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=1},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=1},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=1},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=1},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=1},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=1},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=1},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=1},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=1},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=1},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=1},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=1},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=1},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=1},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=1},
   {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=1},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=1},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=1},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=1},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=1},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=1},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=1},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=1},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=1},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=1},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=1},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=1},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=1},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=1},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=1},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=1},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=1},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=1},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=1},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=1},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=1},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=1},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=1},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=1},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=1},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=1},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=1},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=1},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=1},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=1},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=1},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=1},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=1},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=1},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=1},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=1},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=1},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=1},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=1},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=1},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=1},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=1},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=1},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=1},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=1},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=1},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=1},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=1},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=1},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=1},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=1},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=1},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=1},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=1},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=1},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=1},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=1},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=1},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=1},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=1},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=1},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=1},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=1},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=1},
};

static const dint64_t C[256] = {
   {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=0},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=0},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=0},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=0},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=0},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=0},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=0},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=0},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=0},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=0},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=0},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=0},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=0},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=0},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=0},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=0},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=0},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=0},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=0},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=0},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=0},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=0},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=0},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=0},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=0},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=0},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=0},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=0},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=0},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=0},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=0},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=0},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=0},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=0},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=0},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=0},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=0},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=0},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=0},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=0},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=0},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=0},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=0},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=0},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=0},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=0},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=0},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=0},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=0},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=0},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=0},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=0},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=0},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=0},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=0},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=0},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=0},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=0},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=0},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=0},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=0},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=0},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=0},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=0},
   {.hi = 0x0, .lo = 0x0, .ex = 128, .sgn=0},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=1},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=1},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=1},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=1},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=1},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=1},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=1},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=1},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=1},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=1},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=1},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=1},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=1},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=1},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=1},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=1},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=1},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=1},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=1},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=1},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=1},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=1},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=1},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=1},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=1},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=1},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=1},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=1},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=1},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=1},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=1},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=1},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=1},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=1},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=1},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=1},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=1},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=1},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=1},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=1},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=1},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=1},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=1},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=1},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=1},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=1},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=1},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=1},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=1},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=1},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=1},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=1},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=1},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=1},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=1},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=1},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=1},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=1},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=1},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=1},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=1},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=1},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=1},
   {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=1},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=1},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=1},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=1},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=1},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=1},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=1},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=1},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=1},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=1},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=1},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=1},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=1},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=1},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=1},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=1},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=1},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=1},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=1},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=1},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=1},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=1},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=1},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=1},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=1},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=1},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=1},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=1},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=1},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=1},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=1},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=1},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=1},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=1},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=1},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=1},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=1},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=1},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=1},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=1},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=1},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=1},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=1},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=1},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=1},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=1},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=1},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=1},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=1},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=1},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=1},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=1},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=1},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=1},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=1},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=1},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=1},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=1},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=1},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=1},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=1},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=1},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=1},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=1},
   {.hi = 0x0, .lo = 0x0, .ex = 128, .sgn=0},
   {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=0},
   {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=0},
   {.hi = 0x96a9049670cfae65, .lo = 0xf77574094d3c35c4, .ex = -3, .sgn=0},
   {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=0},
   {.hi = 0xfab272b54b9871a2, .lo = 0x704729ae56d78a37, .ex = -3, .sgn=0},
   {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=0},
   {.hi = 0xaf10a22459fe32a6, .lo = 0x3feef3bb58b1f10d, .ex = -2, .sgn=0},
   {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=0},
   {.hi = 0xe05c1353f27b17e5, .lo = 0xebc61ade6ca83cd, .ex = -2, .sgn=0},
   {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=0},
   {.hi = 0x888e93158fb3bb04, .lo = 0x984156f553344306, .ex = -1, .sgn=0},
   {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=0},
   {.hi = 0xa09ae4a0bb300a19, .lo = 0x2f895f44a303cc0b, .ex = -1, .sgn=0},
   {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=0},
   {.hi = 0xb8442987d22cf576, .lo = 0x9cc3ef36746de3b8, .ex = -1, .sgn=0},
   {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=0},
   {.hi = 0xcf7bca1d476c516d, .lo = 0xa81290bdbaad62e4, .ex = -1, .sgn=0},
   {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=0},
   {.hi = 0xe63374c98e22f0b4, .lo = 0x2872ce1bfc7ad1cd, .ex = -1, .sgn=0},
   {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=0},
   {.hi = 0xfc5d26dfc4d5cfda, .lo = 0x27c07c911290b8d1, .ex = -1, .sgn=0},
   {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=0},
   {.hi = 0x88f59aa0da591421, .lo = 0xb892ca8361d8c84c, .ex = 0, .sgn=0},
   {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=0},
   {.hi = 0x93682a66e896f544, .lo = 0xb17821911e71c16e, .ex = 0, .sgn=0},
   {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=0},
   {.hi = 0x9d7fd1490285c9e3, .lo = 0xe25e39549638ae68, .ex = 0, .sgn=0},
   {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=0},
   {.hi = 0xa73655df1f2f489e, .lo = 0x149f6e75993468a3, .ex = 0, .sgn=0},
   {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=0},
   {.hi = 0xb085baa8e966f6da, .lo = 0xe4cad00d5c94bcd2, .ex = 0, .sgn=0},
   {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=0},
   {.hi = 0xb96841bf7ffcb21a, .lo = 0x9de1e3b22b8bf4db, .ex = 0, .sgn=0},
   {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=0},
   {.hi = 0xc1d8705ffcbb6e90, .lo = 0xbdf0715cb8b20bd7, .ex = 0, .sgn=0},
   {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=0},
   {.hi = 0xc9d1124c931fda7a, .lo = 0x8335241be1693225, .ex = 0, .sgn=0},
   {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=0},
   {.hi = 0xd14d3d02313c0eed, .lo = 0x744fea20e8abef92, .ex = 0, .sgn=0},
   {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=0},
   {.hi = 0xd84852c0a80ffcdb, .lo = 0x24b9fe00663574a4, .ex = 0, .sgn=0},
   {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=0},
   {.hi = 0xdebe05637ca94cfb, .lo = 0x4b19aa71fec3ae6d, .ex = 0, .sgn=0},
   {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=0},
   {.hi = 0xe4aa5909a08fa7b4, .lo = 0x122785ae67f5515d, .ex = 0, .sgn=0},
   {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=0},
   {.hi = 0xea09a68a6e49cd62, .lo = 0x15ad45b4a1b5e823, .ex = 0, .sgn=0},
   {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=0},
   {.hi = 0xeed89db66611e307, .lo = 0x86f8c20fb664b01b, .ex = 0, .sgn=0},
   {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=0},
   {.hi = 0xf314476247088f74, .lo = 0xa5486bdc455d56a2, .ex = 0, .sgn=0},
   {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=0},
   {.hi = 0xf6ba073b424b19e8, .lo = 0x2c791f59cc1ffc23, .ex = 0, .sgn=0},
   {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=0},
   {.hi = 0xf9c79d63272c4628, .lo = 0x4504ae08d19b2980, .ex = 0, .sgn=0},
   {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=0},
   {.hi = 0xfc3b27d38a5d49ab, .lo = 0x256778ffcb5c1769, .ex = 0, .sgn=0},
   {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=0},
   {.hi = 0xfe1323870cfe9a3d, .lo = 0x90cd1d959db674ef, .ex = 0, .sgn=0},
   {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=0},
   {.hi = 0xff4e6d680c41d0a9, .lo = 0xf668633f1ab858a, .ex = 0, .sgn=0},
   {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=0},
   {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=0},
};

/* The following is a degree-13 polynomial with odd coefficients
   approximating sin2pi(x) for 0 <= x < 1/256 with relative error 2^-124.764.
   Generated with sin_accurate.sollya. */
static const dint64_t PS[] = {
  {.hi = 0xc90fdaa22168c234, .lo = 0xc4c6628b80dc1cd0, .ex = 3, .sgn=0}, // 1
  {.hi = 0xa55de7312df295f5, .lo = 0x5dc72f712ae39860, .ex = 6, .sgn=1}, // 3
  {.hi = 0xa335e33bad570e92, .lo = 0x3f3421d4074fb6a9, .ex = 7, .sgn=0}, // 5
  {.hi = 0x9969667315ec2df3, .lo = 0x2c986d9249e41ea2, .ex = 7, .sgn=1}, // 7
  {.hi = 0xa83c1a43f73bfe92, .lo = 0x0, .ex = 6, .sgn=0}, // degree 9
  {.hi = 0xf183a7eef4809d45, .lo = 0x0, .ex = 4, .sgn=1}, // degree 11
  {.hi = 0xf4795452918b54f6, .lo = 0x0, .ex = 2, .sgn=0}, // degree 13
};

/* The following is a degree-14 polynomial with even coefficients
   approximating cos2pi(x) for 0 <= x < 1/256 with relative error 2^-136.297.
   Generated with cos_accurate.sollya. */
static const dint64_t PC[] = {
  {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=0}, // degree 0
  {.hi = 0x9de9e64df22ef2d2, .lo = 0x56e26cd9808c1ab7, .ex = 5, .sgn=1}, // 2
  {.hi = 0x81e0f840dad61d9a, .lo = 0x9980f007d6e9a4f2, .ex = 7, .sgn=0}, // 4
  {.hi = 0xaae9e3f1e5ffcfe2, .lo = 0xa7d6da856a3d7a09, .ex = 7, .sgn=1}, // 6
  {.hi = 0xf0fa83448dd5d7a3, .lo = 0x0, .ex = 6, .sgn=0}, // degree 8
  {.hi = 0xd368f9510253c781, .lo = 0x0, .ex = 5, .sgn=1}, // degree 10
  {.hi = 0xfce9c519909553f6, .lo = 0x0, .ex = 3, .sgn=0}, // degree 12
  {.hi = 0xdb6e0c3401e61fad, .lo = 0x0, .ex = 1, .sgn=1}, // degree 14
};

// put in Y an approximation of sin2pi(X), where X2 approximates X^2
static void
evalPS (dint64_t *Y, dint64_t *X, dint64_t *X2)
{
  mul_dint_21 (Y, X2, PS+6); // degree 13
  add_dint (Y, Y, PS+5);     // degree 11
  mul_dint (Y, Y, X2);
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

// put in Y an approximation of cos2pi(X), where X2 approximates X^2
static void
evalPC (dint64_t *Y, dint64_t *X2)
{
  mul_dint_21 (Y, X2, PC+7); // degree 14
  add_dint (Y, Y, PC+6);     // degree 12
  mul_dint (Y, Y, X2);
  add_dint (Y, Y, PC+5);     // degree 10
  mul_dint (Y, Y, X2);
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
   |Xout - (Xin/(2pi) mod 1)| < 2^-124.34*|Xout| when |Xin| < 2
   |Xout - (Xin/(2pi) mod 1)| < 2^-127           when |Xin| >= 2
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
    /* the ignored part u % 2^64 is at most ulp(X->lo) */
    X->lo = u >> 64;
    u = (u128) X->hi * (u128) T[0];
    X->lo += u;
    X->hi = (u >> 64) + (X->lo < (uint64_t) u);
    /* since X is normalized at input, X->hi >= 2^63, and since T[0] >= 2^61,
       we have X->hi >= 2^(63+61-64) = 2^60, thus the normalize() below
       perform a left shift by at most 3 bits */
    normalize (X);
    /* The error on X->lo is thus bounded by 2^3 = 8, which corresponds to
       a relative error of at most 2^3/2^127 = 2^-124.
       We can get a finer bound as follows:
       (i) if at input X->hi >= 14488038916154245688, then
           X->hi*T0 >= 2^125, thus the left shift is at most by 2 bits,
           and the relative error is bounded by 2^-125.
       (ii) if at input 2^63 <= X->hi < 14488038916154245688, then
           the left shift might be of 3 bits, but after the shift the
           value is lower bounded by 2^63*2^3*T0, thus the relative error
           is bounded by 2^3/(2^63*2^3*T0) < 2^-124.34.
    */
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
  assert (2 <= f && f <= 127);
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
  else /* 65 <= f <= 127 */
  {
    int g = f - 64; /* 1 <= g <= 63 */
    X->hi = (c[3] << g) | (c[2] >> (64 - g));
    X->lo = (c[2] << g) | (c[1] >> (64 - g));
    tiny = (c[1] << g) | (c[0] >> (64 - g));
    /* the ignored part was less than 1 in c[1],
       thus less than 2^g <= 2^63 in tiny */
  }
  /* the approximation error is at most 2 ulps:
     (a) the truncated part in the above shifts, which is less than 1 ulp,
         i.e., less than 2^-128
     (b) the ignored terms hi*T[i+4] + ..., which accumulate to less than
         1 ulp too.
  */
  X->ex = 0;
  /* since X->ex=0, the absolute error of 2 ulps corresponds to 2^-127
     and is not changed after the normalize() call */
  normalize (X);
  /* the worst case (for 2^25 <= x < 2^1024) is X->ex = -61, attained
     for |x| = 0x1.6ac5b262ca1ffp+851 */
  assert (X->ex >= -61);
  if (X->ex < 0) // put the upper -ex bits of tiny into low bits of lo
    X->lo |= tiny >> (64 + X->ex);
}

/* return i and modify X such that Xin = i/2^8 + Xout */
static int
reduce2 (dint64_t *X)
{
  assert (X->ex <= 0);
  int i;
  if (X->ex <= -8)
    i = 0; // X is unchanged
  else
  {
    int sh = 64 - 8 - X->ex;
    i = X->hi >> sh;
    X->hi = X->hi & ((1ul << sh) - 1);
  }
  normalize (X);
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

  double absx = (x > 0) ? x : -x;

  /* now x > 0x1.7137449123ef6p-26 */
  dint64_t X[1];
  dint_fromd (X, absx);

  /* reduce argument */
  reduce (X);

  /* Now X = frac(x/(2pi)) + eps with |eps| < 2^-127, with |X| < 1.
     Write X = i/2^8 + r with r < 2^8. */
  int i = reduce2 (X);
  // if (x0 == TRACE) { printf ("i=%d X2=", i); print_dint (X); }

  /* We use the following identities:
   * sin(x+pi) = -sin(x), thus we can reduce x to [0, pi)
   * sin(x+pi/2) = cos(x), thus we can reduce x to [0, pi/2)
   * sin(pi/2-x) = cos(x), thus we can reduce x to [0, pi/4)
   */

  // approximate sin2pi(x) by sin2pi(i/2^8)*cos2pi(X)+cos2pi(i/2^8)*sin2pi(X)
  dint64_t U[1], V[1], X2[1];
  mul_dint (X2, X, X);       // X2 approximates X^2
  evalPC (U, X2);    // cos2pi(X)
  // if (x0 == TRACE) { printf ("cos2pi(X2)="); print_dint (U); }
  evalPS (V, X, X2); // sin2pi(X)
  // if (x0 == TRACE) { printf ("sin2pi(X2)="); print_dint (V); }
  mul_dint (U, S+i, U); // sin2pi(i/2^8)*cos2pi(X)
  // if (x0 == TRACE) { printf ("sin2pi(i/2^8)*cos2pi(X2)="); print_dint (U); }
  mul_dint (V, C+i, V); // cos2pi(i/2^8)*sin2pi(X)
  // if (x0 == TRACE) { printf ("cos2pi(i/2^8)*sin2pi(X2)="); print_dint (V); }
  add_dint (U, U, V);

  if (x < 0)
    U->sgn = 1 - U->sgn;

  double y = dint_tod (U);

  return y;
}

double
cr_sin (double x)
{
  return sin_accurate (x);
}
