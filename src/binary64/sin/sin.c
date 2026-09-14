/* Correctly-rounded sine function for binary64 value.

Copyright (c) 2022-2026 Paul Zimmermann and Tom Hubrecht and Alexei Sibidanov

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

#include <stdio.h>
#include <assert.h>
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

typedef uint64_t u64;

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

// Return non-zero if a = 0
static inline int
dint_zero_p (const dint64_t *a)
{
  return a->hi == 0;
}

static inline int cmp(int64_t a, int64_t b) { return (a > b) - (a < b); }

static inline int cmpu128 (u128 a, u128 b) { return (a > b) - (a < b); }

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
static const uint64_t _T[20] = {
  0,
  0x28be60db9391054a, // i=1
   0x7f09d5f47d4d3770,
   0x36d8a5664f10e410,
   0x7f9458eaf7aef158,
   0x6dc91b8e909374b8,
   0x1924bba82746487, // i=6
   0x3f877ac72c4a69cf,
   0xba208d7d4baed121,
   0x3a671c09ad17df90,
   0x4e64758e60d4ce7d,
   0x272117e2ef7e4a0e, // i=11
   0xc7fe25fff7816603,
   0xfbcbc462d6829b47,
   0xdb4d9fb3c9f2c26d,
   0xd3d18fd9a797fa8b,
   0x5d49eeb1faf97c5e, // i=16
   0xcf41ce7de294a4ba,
   0x9afed7ec47e35742,
   // 0x1580cc11bf1edaea, // i=19 (only used in reduce_large_acc)
   // 0xfc33ef0826bd0d87, // i=20 (unused)
};

#define U128(l,h) (((u128)h)<<64 | (u128)l)

/* The following is a degree-9 polynomial with odd coefficients
   approximating sin(2*pi*x)/2^7 for 0 <= x < 2^-14 with relative error
   < 2^-126.387. Coefficients of degree 1, 3, 5, 7 are PS[i]/2^128,
   with that of degree 9 is PS[i]/2^64 (fixed point).
   Generated with sinlarge_acc.sollya. */
static const u128 PS[] = {
  // little-endian format
  U128(0x4c4c6628b80dc1cd,0xc90fdaa22168c23),  // degree 1
  U128(0xaee397b871c5480d,0x52aef39896f94afa), // degree 3, implicit - sign
  U128(0x32b2dadb1115ef68,0xa335e33bad570e92), // degree 5
  U128(0x1db4d6855659744c,0x99696673148e38fa), // degree 7, implicit - sign
  U128(0x54125fff015e02e8,0)                   // degree 9
};

/* The following is a degree-6 polynomial with even coefficients
   approximating cos(2*pi*x)/2^7 for 0 <= x < 2^-14 with absolute error
   2^-120.087. Coefficients are PC[i]/2^128 (fixed precision),
   except degree-6 coefficient which is PC[i]/2^64.
   Generated with coslarge_acc.sollya. */
static const u128 PC[] = {
  U128(0xffffffffffffff0f,0x1ffffffffffffff),  // degree 0
  U128(0x95b89954685e30ef,0x277a79937c8bbcb4), // degree 2, implicit - sign
  U128(0xfc971fbbcf438a36,0x81e0f840dad61d03), // degree 4
  U128(0xaae9e3e2d654796c,0),                  // degree 6, implicit - sign
};

static inline u128 mhUU(u128 a, u128 b){
  u64 ah = a>>64, al = a;
  u64 bh = b>>64, bl = b;
  u128 ahbh = (u128)ah*bh;
  u128 ahbl = (u128)ah*bl;
  u128 albh = (u128)al*bh;
  return ahbh += (ahbl>>64)+(albh>>64);
}

static inline void dint_normalize (dint64_t *x)
{
  if (__builtin_expect (x->r == 0, 0)) return;
  uint64_t h = x->r >> 64;
  int sh = (h) ? __builtin_clzll (h) : 64 + __builtin_clzll ((uint64_t) x->r);
  x->r <<= sh;
  x->_ex -= sh;
}

#if 0
// Prints a dint64_t value for debugging purposes
static inline void print_dint(const dint64_t *a) {
  printf("{.hi=0x%"PRIx64", .lo=0x%"PRIx64", .ex=%"PRId64", .sgn=0x%"PRIx64"}\n", a->hi, a->lo, a->ex,
         a->sgn);
}
#endif

/* Put in Y an approximation of sin2pi(X), for 0 <= X < 2^-14,
   where X2 approximates X^2. */
static void
evalPS (dint64_t *Y, dint64_t *X, dint64_t *X2)
{
  u128 u = X->r >> -X->_ex, u2 = X2->r >> -X2->_ex, u2h = u2 >> 64;
  u128 s;
  /* we perform the computation in fixed point, where each variable a is
     interpreted as a/2^128, thus multiplying two variables a and b mean
     taking floor(a*b/2^128) */
  // since signs of coefficients are alternating, and each new coefficient
  // dominates the lower terms, we subtract each time the lower terms from
  // the (absolute value of) the new coefficient
  s = PS[3] - PS[4] * u2h; // PS[3] is degree 7, PS[4] is degree 9
  s = (s>>64) * u2h;
  s = PS[2] - s;
  s = mhUU(s, u2); // multiply by r^2
  s = PS[1] - s;
  s = mhUU(s, u2); // multiply by r^2
  s = PS[0] - s;
  s = mhUU(s, u); // multiply by r
  Y->r = s;
  Y->sgn = X->sgn;
  // since this was for sin(2*pi*r)/2^7, multiply by 2^7
  Y->_ex += 7;
  dint_normalize (Y);
}

/* Put in Y an approximation of cos2pi(X), for 0 <= X < 2^-14,
   where X2 approximates X^2. */
static void
evalPC (dint64_t *Y, dint64_t *X2)
{
  u128 u2 = X2->r >> -X2->_ex, s, u2h = u2 >> 64;
  s = PC[2] - u2h * PC[3];
  s = mhUU(s, u2);
  s = PC[1] - s;
  s = mhUU(s, u2); // multiply by r^2
  s = PC[0] - s;
  Y->r = s;
  // since this was for cos(2*pi*r)/2^7, multiply by 2^7
  Y->_ex += 7;
  dint_normalize (Y);
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
  uint64_t V0, V1;
  if (f == 0) {
    V0 = _T[i];
    V1 = _T[i+1];
  } else {
    V0 = (_T[i] << f) | (_T[i+1] >> (64-f));
    V1 = (_T[i+1] << f) | (_T[i+2] >> (64-f));
  }
  /* Remark: computing directly u with 128-bit arithmetic from _T[i],
     _T[i+1] and _T[i+2] is slower (surely because 128-bit arithmetic is
     emulated.) */
  u128 u = (u128) V1 | (((u128) V0) << 64);
  u = (u128) m * u;
  // round r to nearest, where 0x810000000000000 = 2^59 + 2^52
  static const u128 magic = ((u128) 1 << 112) + 0x810000000000000ull;
  u += magic;
  t.f = (u << 15) >> 75; // next 53 bits of u after the first 15
  *r = t.f * 0x1p-68 - 0x1p-16;
  return u >> 113;
  // since we return 15 bits in i and 53 in h, the accuracy is at most 2^-68
}

static inline double fasttwosum(double x, double y, double *e){
  //  assert (__builtin_fabs (x) >= __builtin_fabs (y));
  double s = x + y, z = s - x;
  *e = y - z;
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
  // the entry U2[1][3] was changed from -0x1.7507dbbbd8fe6p-55 to make
  // x=0x1.997d35866ce04p-14 pass with FMA and rndn
  {0x1.921fb51aeb57cp-13, -0x1.a6e1d4916c435p-67, 0x1.ffffff621619cp-1, -0x1.7507dbbbd8fe7p-55},
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

// deal with 214 exceptional cases in accurate path of moderate branch
// all fail without FMA in at least one rounding mode (either x or -x)
static double __attribute__((noinline))
moderate_exceptions (double x, double y)
{
  static double e[][3] = {
    // for these numbers, sin(x) has >= 45 identical bits after the round bit
    {0x1.d12ed0af1a27fp-26, 0x1.d12ed0af1a27ep-26, 0x1.63676aa7969a8p-134},
    {0x1.250bfe1b082f5p-25, 0x1.250bfe1b082f4p-25, 0x1.951a899fd2c19p-132},
    {0x1.a6a58d55e307cp-25, 0x1.a6a58d55e3079p-25, 0x1.d1d88a022cc12p-130},
    {0x1.f51a62037e956p-25, 0x1.f51a62037e951p-25, -0x1.04f5c2899cd05p-129},
    {0x1.3bacd6561ff5ep-24, 0x1.3bacd6561ff59p-24, 0x1.9ca0914e92066p-129},
    {0x1.4f747439b348bp-24, 0x1.4f747439b3485p-24, 0x1.b29fdd7433125p-129},
    {0x1.9a907c24108e8p-24, 0x1.9a907c24108ddp-24, -0x1.a5478bf8a600fp-129},
    {0x1.b2133770c88d7p-24, 0x1.b2133770c88cap-24, 0x1.5f344ac5b82b3p-129},
    {0x1.c74847a112b67p-24, 0x1.c74847a112b58p-24, 0x1.0417199dde531p-130},
    {0x1.02a3ad2ef6f49p-23, 0x1.02a3ad2ef6f3ep-23, 0x1.362b5aa5d01d8p-128},
    {0x1.3cfc2a006a465p-22, 0x1.3cfc2a006a414p-22, 0x1.1e68bc30335a1p-131},
    {0x1.cbaa95dadb65p-22, 0x1.cbaa95dadb559p-22, -0x1.742593ef4c17bp-128},
    {0x1.fd1109fca2f82p-21, 0x1.fd1109fca2a44p-21, -0x1.379d63c2c675ep-125},
    {0x1.5cd1e10a03f5bp-20, 0x1.5cd1e10a0389cp-20, -0x1.b1860f144f89ap-125},
    {0x1.acd69f89ae57ap-20, 0x1.acd69f89ad8f1p-20, -0x1.b0ddfe3e27547p-126},
    {0x1.e0000000001c2p-20, 0x1.dfffffffff02ep-20, 0x1.dcba692492527p-146},
    {0x1.13841c84df561p-19, 0x1.13841c84de815p-19, 0x1.47a22335e3692p-124},
    {0x1.580646ae65d0ep-19, 0x1.580646ae6432bp-19, 0x1.f1f2e392db463p-125},
    {0x1.e000000000708p-19, 0x1.dffffffffc0b8p-19, 0x1.dcba6924926e6p-139},
    {0x1.7388a06068301p-18, 0x1.7388a06060094p-18, -0x1.32d49d6a8a3c3p-124},
    {0x1.786bbbe176d0ep-18, 0x1.786bbbe16e56ap-18, -0x1.84dcea75432a2p-123},
    {0x1.bd4913a2beabep-18, 0x1.bd4913a2b0a35p-18, 0x1.582219054b3cfp-122},
    {0x1.e000000001c2p-18, 0x1.dffffffff02ep-18, 0x1.dcba692492de2p-132},
    {0x1.eaa3570d59be1p-18, 0x1.eaa3570d46f83p-18, 0x1.c88d1783fb4bbp-123},
    {0x1.6800000002f76p-17, 0x1.67ffffffe54dap-17, 0x1.fd1590076f1cdp-128},
    {0x1.98ec45c2eaa5ep-17, 0x1.98ec45c2bf2c7p-17, 0x1.d77eeaf654169p-122},
    {0x1.e00000000708p-17, 0x1.dfffffffc0b8p-17, 0x1.dcba6924949d1p-125},
    {0x1.2c00000006dddp-16, 0x1.2bffffffc233bp-16, 0x1.1c26f20489067p-122},
    {0x1.680000000bdd8p-16, 0x1.67ffffff95368p-16, 0x1.fd159007734ebp-121},
    {0x1.77cd85e051ebfp-16, 0x1.77cd85dfcaf2bp-16, -0x1.18140d146f5cdp-121},
    {0x1.dbaf171962a04p-16, 0x1.dbaf171850e51p-16, 0x1.afcbd7a64c98cp-121},
    {0x1.58074c0a9021bp-15, 0x1.58074c08f1eep-15, 0x1.065936c29931ep-119},
    {0x1.9a366a206d45cp-15, 0x1.9a366a1daf14bp-15, 0x1.7b4c01529c15p-120},
    {0x1.7bd96ba15c52dp-14, 0x1.7bd96b98a63b4p-14, 0x1.47bbf8c96d894p-120},
    {0x1.f7d7d61b9312fp-14, 0x1.f7d7d6073e9ffp-14, 0x1.55099af61cd9dp-117},
    {0x1.3f69df45a2f3bp-13, 0x1.3f69df30eae3p-13, 0x1.9c7cb344284fap-118},
    {0x1.96350d587e672p-12, 0x1.96350cae09b7fp-12, 0x1.fd72609b5f516p-117},
    {0x1.3bc6ca12143b6p-11, 0x1.3bc6c8d1c549p-11, -0x1.2d3d7b9f94794p-116},
    {0x1.bac75a647203p-11, 0x1.bac756f163452p-11, 0x1.10be21e9a112fp-117},
    {0x1.ce8b994974d0bp-11, 0x1.ce8b955ac69f1p-11, -0x1.38f64e9f1aae3p-115},
    {0x1.933fb67c4d0afp-9, 0x1.933f8ccbc0ea2p-9, 0x1.faaabcb349609p-114},
    {0x1.ab56a7ae04d3bp-9, 0x1.ab5676103cf08p-9, 0x1.ac54455c746b1p-114},
    {0x1.b392f101a024ap-9, 0x1.b392bc7740573p-9, -0x1.fffffffffffffp-63},
    {0x1.b3f2ba40dbc66p-9, 0x1.b3f28593cad2fp-9, 0x1.27c9b58dcaa09p-116},
    {0x1.efe186fe553d9p-9, 0x1.efe13977f9dccp-9, -0x1.f1bfe088e232fp-115},
    {0x1.5dc43f86236ccp-7, 0x1.5dc28c405ade3p-7, -0x1.0965d97639d47p-111},
    {0x1.e17faefac7797p-7, 0x1.e17b3f6bb5e6ep-7, 0x1.159a1d0b5938ep-112},
    {0x1.41db571d96126p-6, 0x1.41d60a76a82edp-6, 0x1.df1d11421219bp-114},
    {0x1.9c412d62c144p-6, 0x1.9c360a8dd681ap-6, -0x1.1ec5e7033170fp-110},
    {0x1.275a3d78c01ecp-5, 0x1.2749dc4d19c6dp-5, -0x1.4974993b17d3dp-110},
    {0x1.4cd45ddee2881p-5, 0x1.4cbced7e1d1aap-5, -0x1.5155493b537eap-110},
    {0x1.69949b3d51fb1p-5, 0x1.69768dc89bbp-5, 0x1.1474147e16d5p-113},
    {0x1.9283586503fep-5, 0x1.9259e3708bd3ap-5, -0x1.c8482a5a8702bp-114},
    {0x1.d7bdcd778049fp-5, 0x1.d77b117f230d6p-5, -0x1.aa47015e8c53ap-114},
    {0x1.0023629fc9899p-4, 0x1.fff150c9d77a2p-5, -0x1.4c880c955dffcp-108},
    {0x1.21857ad584f7fp-4, 0x1.2147c702305d4p-4, -0x1.f672750ebb561p-109},
    {0x1.2e36813a9874p-4, 0x1.2df0542154f1bp-4, 0x1.62b3f158c6d6cp-108},
    {0x1.456ac98461b72p-4, 0x1.45132d642e713p-4, -0x1.70edaec74c0b3p-108},
    {0x1.5231b416ba885p-4, 0x1.51cf5db1b1956p-4, 0x1.e03fb05ea92fp-108},
    {0x1.7f4ea0f3bbc6fp-4, 0x1.7ebf782f9ca35p-4, -0x1.40ea7e2a2e8cbp-107},
    {0x1.a202b3fb84788p-4, 0x1.a1490c8c06ba7p-4, -0x1.59e7bbfe7bafep-110},
    {0x1.c49ac7cde7b4cp-3, 0x1.c0edeb94cef34p-3, 0x1.c1d17c18ed011p-107},
    {0x1.d5064e6fe82c5p-3, 0x1.d0ef799001ba9p-3, 0x1.7a93a488547a6p-109},
    {0x1.dd04b12d498c6p-3, 0x1.d8b784ad454cp-3, 0x1.0131008b93463p-107},
    {0x1.e3095cae52dd7p-3, 0x1.de920f7a4c509p-3, -0x1.98a7842874906p-108},
    {0x1.1c63df63d59dp-2, 0x1.18bf91b163125p-2, -0x1.fp-56},
    {0x1.223c48cd64801p-2, 0x1.1e5d765189691p-2, 0x1.aa13d9a49ca2ap-107},
    {0x1.50954b7bbf87bp-2, 0x1.4a8e1a96e38e3p-2, 0x1.a24f4f00de32bp-108},
    {0x1.e05b0e0a809bcp-2, 0x1.ceee68154d1c9p-2, -0x1.e7e4ad8afd722p-108},
    {0x1.ed25c5eb8c916p-2, 0x1.da4e0e6c717a5p-2, 0x1.5d145a6990971p-109},
    {0x1.fe767739d0f6dp-2, 0x1.e9950730c4696p-2, -0x1.e8b7356bc3fcp-121},
    {0x1.ff92a8ca216cdp-2, 0x1.ea8e8fdf47549p-2, 0x1.104fcd1a5d3d4p-105},
    {0x1.3eb4df7c5532ap-1, 0x1.2a8517f17245dp-1, -0x1.811124fe53558p-104},
    {0x1.41516c909749cp-1, 0x1.2ca340fa7150ap-1, -0x1.e067d5234e30fp-106},
    {0x1.d98c4c612718dp-1, 0x1.98dcd09337793p-1, -0x1.fp-55},
    {0x1.02c2f02bd16d5p+0, 0x1.b1cac622470fep-1, 0x1.45d3df4eaefefp-104},
    {0x1.3aa301f6ebb1ep+0, 0x1.e264357ea0e29p-1, 0x1.38aa76739f9ecp-106},
    {0x1.640ac66708cp+0, 0x1.f7ba32c220693p-1, 0x1.b5f34f48f617dp-107},
    {0x1.bbfa05708792dp+0, 0x1.f92c3e0cf3454p-1, 0x1.7b2d930ede991p-107},
    {0x1.4dbe000d5c1d2p+1, 0x1.04b33eeedd2bap-1, 0x1.41ce88e6ad5c4p-108},
    {0x1.6756745770a51p+1, 0x1.4ff350e412821p-2, 0x1.5d0547bfa347p-110},
    {0x1.6e6198df13b76p+1, 0x1.1a3d49d1edd09p-2, -0x1.f22a3583f4cfdp-107},
    {0x1.91f834a9fd65bp+1, 0x1.3c04cd272f3ap-10, 0x1.e8aeaf00deebp-110},
    {0x1.920745cc24d5ep+1, 0x1.86f77f7fce01bp-11, 0x1.9ea9a411162a9p-112},
    {0x1.aa0b46aa9cc59p+1, -0x1.7c7fcf52c49ffp-3, 0x1.3763216dcf618p-107},
    {0x1.3f9f5c0e54006p+2, -0x1.ebd13d9a1a8a4p-1, 0x1.f402e09fd7859p-105},
    {0x1.8bb935a75e95ap+2, -0x1.98f132485adf4p-4, 0x1.89902bfc4a871p-111},
    {0x1.8ddc167f304fcp+2, -0x1.10b403a8068fcp-4, 0x1.90477149d0335p-110},
    {0x1.9328b6f1a41d5p+2, 0x1.08feb81b316d2p-6, -0x1.364b004fa2531p-114},
    {0x1.e2ed3c4a32bc7p+3, 0x1.280804d036036p-1, -0x1.4a4be32f2ed92p-105},
    {0x1.31ce338ed5035p+4, 0x1.0a80445c9ec1fp-2, -0x1.3db17469eea23p-106},
    {0x1.6054c040f738cp+4, -0x1.e3f48ef9add7p-6, 0x1.3839ced9dbe35p-110},
    {0x1.79bd7d6a366c3p+4, -0x1.ff706171be64dp-1, 0x1.3d3586bdff285p-105},
    {0x1.c98f08086fe46p+4, -0x1.451d394832e3fp-2, -0x1.0d48381eda5p-105},
    {0x1.37545908e1661p+5, 0x1.e04d4acb7e603p-1, 0x1.628c511101268p-104},
    {0x1.3ac50dcd51f43p+5, 0x1.fe828f2acae9bp-1, -0x1.89dce31ffd065p-107},
    {0x1.5b47aec7b6254p+5, -0x1.1547cec05d854p-1, -0x1.7a11836f8d6a9p-104},
    {0x1.8b6ddd78c07c5p+5, -0x1.7c2cb0efde00ep-1, 0x1.23fc50138de96p-104},
    {0x1.c1c27d6b9a3a1p+5, -0x1.4a8ff4dfae289p-2, 0x1.fp-56},
    {0x1.c82f99743c4f8p+5, 0x1.d3ed3e0f6b30cp-2, 0x1.4b9048236b732p-106},
    {0x1.e9cf095b6b0b9p+5, -0x1.ffafe375ca99dp-1, -0x1.3e80ba4507c57p-105},
    {0x1.10efabd43f3f8p+6, -0x1.8af073e535056p-1, -0x1.878347340a4ccp-108},
    {0x1.2fa8afa9cb896p+6, 0x1.f9b327ce15d33p-2, 0x1.a3244ef240646p-106},
    {0x1.366b8efaab0fbp+6, 0x1.9be32f88e58afp-1, 0x1.825ba39c557fp-105},
    {0x1.54fd5658615c7p+6, -0x1.a5a7a38e8ae7dp-2, -0x1.fffffffffffffp-56},
    {0x1.68edefc6c4deep+6, 0x1.8898d69cbcd6fp-1, -0x1.b27b5df50b5a3p-109},
    {0x1.d78583e09cd36p+6, -0x1.feb8e9b51a94cp-1, -0x1.62308d7c1d04ep-104},
    {0x1.14707a375dbb6p+7, -0x1.5498521bb1de3p-7, 0x1.72a227d283221p-111},
    {0x1.a26cc8aa85837p+7, 0x1.e9a67c0e36ffp-1, 0x1.c963f5ff4486cp-104},
    {0x1.dc246138efd53p+7, -0x1.45e6ee297daaap-1, -0x1.47f248201931p-104},
    {0x1.dd90d3fdf1aebp+7, 0x1.659057bc7ad3bp-6, -0x1.250803afeb15bp-115},
    {0x1.1ae62cc66d317p+8, 0x1.3dc0a15265623p-3, -0x1.555807e985bdp-107},
    {0x1.399de575e3ffcp+8, -0x1.0853bf49f48ecp-1, -0x1.077193b84b9dp-105},
    {0x1.a6bd78ce2c7cap+8, 0x1.f634a592fbcc7p-1, -0x1.00af2ebe1a146p-104},
    {0x1.ad0e140380c51p+8, 0x1.f2cbe0b19f8fap-1, 0x1.b462db092fa7dp-105},
    {0x1.bd8248a58eb58p+8, -0x1.1fed829ad52ddp-1, -0x1.4e9cd4d29e6eep-112},
    {0x1.c418e6893ae8ap+8, -0x1.26d9965a85982p-2, -0x1.35c98428aad96p-108},
    {0x1.ce3c94c1a1b94p+8, -0x1.a3e7f1c95fe2fp-2, -0x1.055f76e49c785p-111},
    {0x1.f10d75e8dd823p+8, 0x1.42510abd0daf5p-1, 0x1.bee63ebe5d0cep-106},
    {0x1.fb25557b798c9p+8, -0x1.f3814b654ffdap-1, -0x1.acdd7a9a46b91p-108},
    {0x1.167a64e430d1p+9, -0x1.8f3e7fb33ba9bp-1, 0x1.12d45e6320387p-107},
    {0x1.4853ba521c795p+9, -0x1.f58502aa2f0bep-5, -0x1.e449950db3b1bp-111},
    {0x1.a3f4d96bec30cp+9, -0x1.ca0f88edaf363p-1, 0x1.323ec05794224p-104},
    {0x1.bc9a4937749aep+9, -0x1.12e83cc4a9ee9p-3, 0x1.7ed8b800d0193p-107},
    {0x1.ce3420174a86ap+9, 0x1.67b601146cd34p-1, 0x1.247e6a9e24324p-104},
    {0x1.1984cf8c7f5c4p+11, 0x1.73d186b6ec3c8p-2, -0x1.0977a0839ba9fp-106},
    {0x1.216584aec736bp+11, 0x1.7159b0b10c469p-3, -0x1.918b9015e8c28p-107},
    {0x1.3643cf7c19355p+11, 0x1.081f62648dd1bp-2, -0x1.813b6472badd8p-107},
    {0x1.66dc0df9d9395p+11, -0x1.06ca0d6121ec5p-1, 0x1.f884a06a62a42p-105},
    {0x1.77e27dc713297p+11, -0x1.15e0b9725f536p-1, 0x1.b543afc21de79p-110},
    {0x1.ecaecf6b1844bp+11, 0x1.e390d3aa60845p-1, 0x1.f0869f433e539p-107},
    {0x1.49a9bedab6603p+12, 0x1.fe8bcd7736f0dp-4, 0x1.ffffffffffffep-58},
    {0x1.76f32cd757e0bp+12, -0x1.e4d2eb11327aap-1, 0x1.9fba18fb6fe43p-105},
    {0x1.a0b5478918088p+12, 0x1.8774050b578e5p-1, 0x1.1c01a1e42782ep-105},
    {0x1.a245d28deeedbp+12, 0x1.64fdd683ccdd7p-1, -0x1.07569b7eb8a37p-104},
    {0x1.ab3866badab77p+12, -0x1.18d117d156af5p-1, 0x1.fffffffffffffp-55},
    {0x1.b811e18ffb9cp+12, -0x1.7247f80d3110fp-1, 0x1.410c7a61da58ep-111},
    {0x1.b98bf86ee22cp+13, -0x1.f76b210273e48p-1, 0x1.a0c8f9743e6cep-108},
    {0x1.7c1e1cf387129p+14, -0x1.a4f5418717703p-1, 0x1.6896ce4ef676ap-105},
    {0x1.f63eb9dc79f2ep+14, -0x1.c705286a68586p-1, -0x1.26df58ca25281p-105},
    {0x1.f7f510d27168ap+14, 0x1.fd27f03a7a4b4p-1, 0x1.525169ee6e557p-105},
    {0x1.05bcab0c7116ep+15, 0x1.85554df3ec6b4p-2, -0x1.4e7fd187ed6e3p-106},
    {0x1.184b1c13010b1p+15, 0x1.12e143b7ae92ep-1, -0x1.374b39c6cf3e8p-104},
    {0x1.1dcc56d882b0dp+15, 0x1.fd1fda7302095p-1, -0x1.3764f7203fa61p-107},
    {0x1.25a073e80e68cp+15, -0x1.f3f4b5032b111p-1, -0x1.22af4f2e39003p-105},
    {0x1.420a5f7b04247p+15, -0x1.5e6eb3fd00b8fp-2, 0x1.bb222cbf0ae0cp-112},
    {0x1.4fe296a7a1a3ap+15, -0x1.209524923381bp-1, -0x1.ed97a2b0913ddp-106},
    {0x1.5cfbf1eded8fbp+15, 0x1.4f2257eb476b3p-2, 0x1.71f3668e19a38p-107},
    {0x1.889fe79bbd1a3p+15, 0x1.ae975b3b2b213p-4, -0x1.d772b2c798e38p-108},
    {0x1.cadabd4f3f112p+15, -0x1.eca9b1d347777p-1, 0x1.06bbbd1875b4bp-106},
    {0x1.2c2119f5b8437p+16, 0x1.79dd09e87ab32p-1, -0x1.5f96d57dc7b0ep-110},
    {0x1.36f8050dfba07p+16, 0x1.fac85db1ec4fp-5, -0x1.a8b11007d9008p-110},
    {0x1.a66e19f35cf61p+16, 0x1.9e3ca40d4fe86p-1, -0x1.92f4669431b5bp-105},
    {0x1.1fbc42d445c64p+17, -0x1.f075f0488b1dfp-1, -0x1.e5d4b0db372e9p-105},
    {0x1.832f57afef571p+17, -0x1.bb5c342ad4866p-1, 0x1.2097af06e20a1p-107},
    {0x1.d93d1bee0f241p+17, -0x1.0414ad5e2e479p-2, -0x1.adc4741c1b685p-107},
    {0x1.01ee93b055ec6p+18, 0x1.735fb614cb857p-1, -0x1.95920b81e76dcp-105},
    {0x1.17f9b4fd7de4dp+18, -0x1.2f76a07410eb4p-1, -0x1.f567734984f1fp-105},
    {0x1.21bdc29ab528dp+18, 0x1.c62690b7d5d85p-4, -0x1.67d1e62a3654p-107},
    {0x1.66aac48b525f4p+18, -0x1.90a229f044c2bp-1, 0x1.7be73ec2309f5p-104},
    {0x1.77f4b6490e483p+18, 0x1.f276bfcc06527p-1, -0x1.7b41afb4a9cd5p-104},
    {0x1.b821d8a71218bp+18, 0x1.3123a70f01f36p-1, 0x1.39a091704bce9p-105},
    {0x1.019d3293a13b8p+19, 0x1.6cf24f7cca579p-1, -0x1.853825e416bb2p-106},
    {0x1.085598417407dp+19, -0x1.37495f2965638p-1, -0x1.47abc76f99d9fp-104},
    {0x1.7774a26a78f4dp+19, 0x1.a4ed4fbe264cap-11, 0x1.36c9178e4d07ep-110},
    {0x1.fbcef53002aefp+19, -0x1.d7693b7007dd5p-1, -0x1.adfe6fa304199p-107},
    {0x1.be5309b4f110cp+20, 0x1.fffe2d9a39c1fp-1, 0x1.b384b22af99bap-106},
    {0x1.d93a1823f3881p+20, -0x1.a5c8b9fd8b9d3p-6, 0x1.8bb1b7e5d0d34p-111},
    {0x1.fa793e7a91851p+20, 0x1.f372655c7afb8p-3, 0x1.7d30b2f58c0cep-106},
    {0x1.3a6afe77b2d01p+21, 0x1.fd51a93d4c345p-1, 0x1.fffffffffffffp-55},
    {0x1.b848d508ee21ap+21, -0x1.ff148e505203ep-1, 0x1.180d6335d4bf7p-104},
    {0x1.bda1aa19cccb5p+21, 0x1.96a390941a1aep-1, -0x1.7e01c5ee463e6p-105},
    {0x1.c4fe5cfcf07cp+21, -0x1.b382be3199d86p-1, -0x1.b872771c91139p-105},
    {0x1.e1afb7cedcff1p+21, 0x1.ddbf19deb6076p-2, 0x1.6299d0d5fd56p-106},
    {0x1.1b25a9807cee2p+22, -0x1.42b5a1adc680dp-1, -0x1.fffffffffffffp-55},
    {0x1.22917a498c836p+22, 0x1.ffd7e0aaecf89p-1, 0x1.2e75d582f65ap-108},
    {0x1.6131e69b9fb8p+22, 0x1.c816ebcd50d2ep-1, -0x1.4f15c864d4d75p-104},
    {0x1.5e228c8baa93cp+23, 0x1.e32c546ab2541p-3, 0x1.10397ab87c74bp-110},
    {0x1.c27086872677fp+23, 0x1.06e75724debd3p-3, 0x1.4064fa785a243p-107},
    {0x1.1cc0782f531cdp+24, 0x1.83a0b2bb097b6p-1, -0x1.b0b63678f2e17p-106},
    {0x1.1e3b5f4f04c5fp+24, 0x1.ffe321a7aba3dp-1, -0x1.e135b1b26fdcap-106},
    {0x1.7ae8db3183dp+24, 0x1.a4f52902fcf34p-2, -0x1.4af1b599485e2p-106},
    {0x1.e130f2f211197p+24, 0x1.a106093208996p-1, 0x1.43714f6352b66p-104},
    {0x1.1ee987b371a4p+25, 0x1.af75397f7fa35p-1, -0x1.5cb828223ed66p-105},
    {0x1.4585b9930ccffp+25, -0x1.d86e830acdf3fp-1, -0x1.b382fda9b6ff6p-106},
    {0x1.5a853150b67b8p+25, 0x1.befefb7d4d3c8p-1, 0x1.14b77cdce797ap-105},
    {0x1.991a57be5cadcp+25, 0x1.14e67421a4636p-1, -0x1.24449af8d9e21p-105},
    {0x1.b48f087b067f1p+25, 0x1.dcc907c94dbb5p-1, -0x1.d762e51a18f61p-108},
    {0x1.b71fda220e13cp+25, -0x1.fd09705b8096fp-1, -0x1.11226711abaf7p-105},
    {0x1.fca67be1d4c5dp+25, -0x1.fd6332010ef33p-1, -0x1.fp-55},
    {0x1.20b6ac558eb5fp+26, -0x1.0184b0a507793p-2, -0x1.43d2e7535ce84p-106},
    {0x1.388526bd17276p+26, 0x1.578546d46d091p-3, -0x1.df3dc6b46e527p-107},
    {0x1.acf716896cba3p+26, -0x1.c8700372d6f7bp-1, 0x1.ffffffffffffep-55},
    {0x1.ad40ecd7a777bp+26, 0x1.1c4942eacc3f7p-1, 0x1.d2ab120b22ce4p-106},
    {0x1.ddc1231000bfp+26, 0x1.c9caeffdca98dp-1, 0x1.8dce5698e77eep-106},
    {0x1.ed8615b697b8cp+26, -0x1.bd4ffc9e2e402p-3, 0x1.0db138687b0b3p-106},
    {0x1.009003ab60bf1p+27, 0x1.fe235f20ea46p-1, -0x1.040a8a8ca17bep-104},
    {0x1.342892a22c1aap+27, -0x1.bf51d5bffa216p-2, -0x1.389ef8732b4f4p-107},
    {0x1.7d04a65e32e63p+27, -0x1.b762ee3b4c2a9p-1, -0x1.6bc16745e68p-104},
    {0x1.b597e803ac69cp+27, 0x1.73fc31300d212p-3, -0x1.3a2c50850bc59p-106},
    {0x1.54061da9ce9a8p+28, 0x1.a8409d2b8f93bp-1, 0x1.fffffffffffffp-55},
    {0x1.62fcf3016ac5p+28, 0x1.ffb489e611b6p-1, -0x1.8310b7215a9f8p-107},
    {0x1.b746170d64856p+28, 0x1.8d0156f261a66p-1, -0x1.fp-55},
    {0x1.2d76bd73bdb4bp+29, 0x1.f93b20d5eb003p-1, -0x1.4dc8679d91eb9p-105},
    {0x1.3f02b26df5f02p+29, -0x1.841be3d2a23f5p-3, 0x1.ae2951c6bdd9bp-110},
    {0x1.445bd4e63db0bp+29, -0x1.e7eee6dc7e973p-2, 0x1.3ab14c833b15p-105},
    {0x1.a3185a4ed864fp+29, -0x1.69991c5c86111p-5, 0x1.50a9dc960a30fp-111},
    {0x1.f90431e0f8215p+29, 0x1.838c0b4fac125p-1, -0x1.31db8329db29cp-108},
    {0x1.098bad0eee782p+30, 0x1.fbac4d8c08565p-1, -0x1.80130739ed7dbp-105},
    {0x1.0b8ecd510ab3bp+30, 0x1.162e36080f948p-3, 0x1.376adda23bcb5p-107},
    {0x1.3711314ca6203p+30, -0x1.1b1906690f41fp-2, 0x1.13b877cf7ec02p-106},
    {0x1.6a4eb91f233bdp+30, -0x1.8007612316001p-1, -0x1.30b968b2da63ep-105},
    {0x1.742fc441fa8b3p+30, -0x1.f7d9e491a7afbp-1, -0x1.56d6397a8a75dp-106},
    {0x1.f2d1b8ac49c5bp+30, 0x1.5b889d955c0d4p-1, -0x1.8a21d6f4fda52p-105},
  };
  double ax = __builtin_fabs (x);
  unsigned int a = 0, b = sizeof(e)/(3*sizeof(double)), c;
  // invariant: e[i][a] <= ax < e[i][b]
  while (a + 1 < b) {
    c = (a + b) / 2;
    if (e[c][0] <= ax)
      a = c;
    else
      b = c;
  }
  if (ax == e[a][0])
    return (x == ax) ? e[a][1] + e[a][2] : - e[a][1] - e[a][2];
  return y;
}

static const double pih = -0x1.921fb54442d18p-13;
static const double pil = -0x1.1a62633145c07p-67;
static const double pis = 0x1.f1976b7ed8fbcp-123;

// accurate path for 0x1.95e4p+1 <= |x| < 2^31,
// reusing values computed in the fast path
static double __attribute__((noinline))
sin_accurate_moderate (double *ctx, int sbit, int i1, int i2)
{
  double x = ctx[0], k = ctx[1], Sh = ctx[2], Sl = ctx[3], rh = ctx[4],
    rl = ctx[5];
  double rs = k * pis + __builtin_fma (k, pil, -rl);
  rh = fasttwosum (rh, rl, &rl);
  rl += rs;
  rh = fasttwosum (rh, rl, &rl);

  // Sh approximates sin(t1+t2)
  double c1h, c1l, c2h, c2l;
  // c1h+c1l approximates cos(t1)*cos(t2)
  c1h = muldd (U1[i1][2], U1[i1][3], U2[i2][2], U2[i2][3], &c1l);
  // c2h+c2l approximates sin(t1)*sin(t2)
  c2h = muldd (U1[i1][0], U1[i1][1], U2[i2][0], U2[i2][1], &c2l);
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
  static const double ps[] = {1.0, -0x1.5555555555555p-3,
         -0x1.555555551de06p-57, 0x1.1111111111111p-7, -0x1.a01a006eb9947p-13};
  static const double pc[] = {-0x1p-1, 0x1.5555555555555p-5,
                              -0x1.6c16bcc416d44p-10};
  double sh, sl, t;
  // since |ps[4]*r^7/ps[0]*r| < 2^-92, we can compute ps[4]*r^2 as double
  // since |ps[3]*r^5/ps[0]*r| < 2^-60, we can still compute in binary64
  // add ps[3]
  sh = ps[3] + r2h * ps[4];
  sh = r2h * sh;
  // sh approximates ps[3]*r^2+ps[4]*r^4
  // add ps[1]+ps[2]
  sh = fasttwosum (ps[1], sh, &sl);
  sl += ps[2];
  sh = muldd (r2h, r2l, sh, sl, &sl);
  // sh+sl approximates (ps[1]+ps[2])*r^2+ps[3]*r^4+ps[4]*r^6
  // add ps[0]
  sh = fasttwosum (ps[0], sh, &t);
  sl += t;
  // multiply by rh+rl
  sh = muldd (rh, rl, sh, sl, &sl);
  // now sh+sl approximates sin(rh+rl)
  double ch, cl;
  // since |pc[2]*r^6| < 2^-89, we can compute pc[2]*r^2 as double
  // since |pc[1]*r^4| < 2^-57, we can still compute in binary64
  // add pc[1]
  ch = pc[1] + pc[2] * r2h;
  ch = r2h * ch;
  // ch approximates pc[1]*r^2+pc[2]*r^4
  // add pc[0]
  ch = fasttwosum (pc[0], ch, &cl);
  ch = muldd (r2h, r2l, ch, cl, &cl);
  // now ch+cl approximates cos(rh+rl)-1

  /* We now have to compute (Sh+Sl)*(1+ch+cl) + (Ch+Cl)*(sh+sl),
     where Sh+Sl approximates sin(z) where z = pi*(i1/2^7+i2/2^14),
     Ch+Cl approximates cos(z),
     sh+sl approximates sin(r) and ch+cl approximates cos(r)-1.
     Since |r| <= 2^-13.339, |sh+sl| <= 2^-13.339 and |1-cos(r)| < 2^-27.678.
  */
  ch = muldd (Sh, Sl, ch, cl, &cl);
  sh = muldd (Ch, Cl, sh, sl, &sl);
  // we can use fasttwosum since either sh=0 or |sh| >= |ch|
  sh = fasttwosum (sh, ch, &t);
  sl += cl + t;
  // now |sh+sl| < 2^-13.339 + 2^-27.678 < 2^-13.338
  /* Now add Sh+Sl: since Sh+Sl is either zero or |Sh+Sl| >= sin(pi/2^14),
     we can use fasttwosum. */
  sh = fasttwosum (Sh, sh, &t);
  // now both Sl and t are of the order of ulp(sh), while sl is smaller
  /* Add sl and Sl first: since they are of the same order of magnitude,
     it might be that the sum is exact. And indeed it gives less failures. */
  sl = (sl + Sl) + t;
  double ret = (sbit == 0) ? sh + sl : - sh - sl;
  // check worst cases
  b64u64_u z = {.f = sl};
  if (__builtin_expect (((z.u + 2) & 0x7fffffffffffull) <= 3, 0))
    return moderate_exceptions (x, ret);
  return ret;
}

// fast path for 0x1.95e4p+1 <= |x| < 2^31
// ax = |x| and eps is the error bound for the rounding test
// see proof of correctness in sin.pdf
static inline double
cr_sin_moderate (double x, int sbit)
{
  double ax = __builtin_fabs(x);
  static const double invpi = 0x1.45f306dc9c883p+12;
  // |invpi/2^14 - 1/pi| < 2^-55.496
  double k = __builtin_roundeven (invpi * ax);
  // |2^14*(pih + pil) + pi| < 2^-108.041
  double rh = __builtin_fma (k, pih, ax), rl = k * pil; // rh is exact

  double r = rh + rl; // |r| < 2^-13.339 (see sin.pdf)
  double r2 = r * r;
  int64_t j = k;
  sbit ^= (j >> 14) & 1; // reduction by an odd multiple of pi?
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
  double sh =  r * ( 1.0 - 0x1.55555553068fp-3 * r2);
  double ch = r2 * (-0.5 + 0x1.55555553bfd3p-5 * r2);
  double fh = Sh, fl = Sl + Sh*ch + Ch*sh;
  static const double Sgn[] = {1.0, -1.0};
  const double eps = 0x1.dep-64, eps2 = 0x1.dep-63;
  fh = Sgn[sbit] * fh;
  fl = Sgn[sbit] * fl - eps;
  double lb = fh + fl, ub = fh + (fl + eps2);
  if (__builtin_expect (ub == lb, 1)) return lb;
  double ctx[] = {x, k, Sh, Sl, rh, rl};
  return sin_accurate_moderate (ctx, sbit, i1, i2);
}

// argument reduction for |x| >= 2^31 (accurate path)
// return k and r such that
// x/(2pi) mod 1 = k/2^13 + r + eps with |r| <= 2^-14
// and 0 <= eps < 2^-128 + 2^-139 < 2^-127.999
static uint64_t
reduce_large_acc (dint64_t *r, double x)
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
  // shift the table entries
  uint64_t V0, V1, V2;
  if (f == 0) {
    V0 = _T[i];
    V1 = _T[i+1];
    V2 = _T[i+2];
  } else {
    V0 = (_T[i] << f) | (_T[i+1] >> (64-f));
    V1 = (_T[i+1] << f) | (_T[i+2] >> (64-f));
    V2 = (_T[i+2] << f) | (_T[i+3] >> (64-f));
  }
  /* m*V0 contributes to 64 bits to the fractional part of x/(2pi),
     from weight 2^-1 to 2^-64,
     m*V1 contributes to 53+64 bits from weight 2^-12 to 2^-128,
     m*V2 contributes to 53+64 bits from weight 2^-76 to 2^-192,
     the ignored part m*(V3+V4+...) contributes to less than 2^-139 */
  u128 u = (u128) V1 | (((u128) V0) << 64);
  u = (u128) m * u;
  u128 v = (u128) m * (u128) V2;
  u += v >> 64; // add contribution of m*V2 from 2^-76 to 2^-128
  // the ignored part of v contributes to less than 2^-128
#define SHIFT 13
  uint64_t k = u >> (128-SHIFT);
  u = u << SHIFT; // ignore leading SHIFT bits (returned in k)
  // round k to nearest to have |r| < 2^-14
  int neg = u >> 127;
  if (neg) { // add 1 to k and subtract 1/2^13 to r
    k = (k+1) & ((1ull<<SHIFT)-1);
    u = -u;
  }
  // now store u in r
  if (u == 0) { cp_dint (r, &ZERO); return k; }
  // compute the number of leading zeros in u
  int sh = (u>>64) ? __builtin_clzll (u>>64) : 64 + __builtin_clzll((uint64_t) u);
  r->r = u << sh; // make significand normalized
  r->_ex = -SHIFT - sh;
  r->_sgn = neg;
  return k;
}

/* Table containing 128-bit approximations of sin(pi*i/2^6) for 0 <= i < 64
   (to nearest).
   Each entry is to be interpreted as (hi/2^64+lo/2^128)*2^ex*(-1)^sgn.
   Generated with computeS1() from sin.sage. */
static const dint64_t S1[64] = {
  {.hi = 0x0, .lo = 0x0, .ex = 0, .sgn=0},
  {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=0},
  {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=0},
  {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=0},
  {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=0},
  {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=0},
  {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=0},
  {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=0},
  {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=0},
  {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=0},
  {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=0},
  {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=0},
  {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=0},
  {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=0},
  {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=0},
  {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=0},
  {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=0},
  {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=0},
  {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=0},
  {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=0},
  {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=0},
  {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=0},
  {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=0},
  {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=0},
  {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=0},
  {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=0},
  {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=0},
  {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=0},
  {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=0},
  {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=0},
  {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=0},
  {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=0},
  {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=0},
  {.hi = 0xffb10f1bcb6bef1d, .lo = 0x421e8edaaf59453e, .ex = 0, .sgn=0},
  {.hi = 0xfec46d1e89292cf0, .lo = 0x41390efdc726e9ef, .ex = 0, .sgn=0},
  {.hi = 0xfd3aabf84528b50b, .lo = 0xeae6bd951c1dabbe, .ex = 0, .sgn=0},
  {.hi = 0xfb14be7fbae58156, .lo = 0x2172a361fd2a722f, .ex = 0, .sgn=0},
  {.hi = 0xf853f7dc9186b952, .lo = 0xc7adc6b4988891bb, .ex = 0, .sgn=0},
  {.hi = 0xf4fa0ab6316ed2ec, .lo = 0x163c5c7f03b718c5, .ex = 0, .sgn=0},
  {.hi = 0xf1090827b43725fd, .lo = 0x67127db35b287316, .ex = 0, .sgn=0},
  {.hi = 0xec835e79946a3145, .lo = 0x7e610231ac1d6181, .ex = 0, .sgn=0},
  {.hi = 0xe76bd7a1e63b9786, .lo = 0x125129529d48a92f, .ex = 0, .sgn=0},
  {.hi = 0xe1c5978c05ed8691, .lo = 0xf4e8a8372f8c5810, .ex = 0, .sgn=0},
  {.hi = 0xdb941a28cb71ec87, .lo = 0x2c19b63253da43fc, .ex = 0, .sgn=0},
  {.hi = 0xd4db3148750d1819, .lo = 0xf630e8b6dac83e69, .ex = 0, .sgn=0},
  {.hi = 0xcd9f023f9c3a059e, .lo = 0x23af31db7179a4aa, .ex = 0, .sgn=0},
  {.hi = 0xc5e40358a8ba05a7, .lo = 0x43da25d99267326b, .ex = 0, .sgn=0},
  {.hi = 0xbdaef913557d76f0, .lo = 0xac85320f528d6d5d, .ex = 0, .sgn=0},
  {.hi = 0xb504f333f9de6484, .lo = 0x597d89b3754abe9f, .ex = 0, .sgn=0},
  {.hi = 0xabeb49a46764fd15, .lo = 0x1becda8089c1a94c, .ex = 0, .sgn=0},
  {.hi = 0xa267992848eeb0c0, .lo = 0x3b5167ee359a234e, .ex = 0, .sgn=0},
  {.hi = 0x987fbfe70b81a708, .lo = 0x19cec845ac87a5c6, .ex = 0, .sgn=0},
  {.hi = 0x8e39d9cd73464364, .lo = 0xbba4cfecbff54867, .ex = 0, .sgn=0},
  {.hi = 0x839c3cc917ff6cb4, .lo = 0xbfd79717f2880abf, .ex = 0, .sgn=0},
  {.hi = 0xf15ae9c037b1d8f0, .lo = 0x6c48e9e3420b0f1e, .ex = -1, .sgn=0},
  {.hi = 0xdae8804f0ae6015b, .lo = 0x362cb974182e3030, .ex = -1, .sgn=0},
  {.hi = 0xc3ef1535754b168d, .lo = 0x3122c2a59efddc37, .ex = -1, .sgn=0},
  {.hi = 0xac7cd3ad58fee7f0, .lo = 0x811f953984eff83e, .ex = -1, .sgn=0},
  {.hi = 0x94a03176acf82d45, .lo = 0xae4ba773da6bf754, .ex = -1, .sgn=0},
  {.hi = 0xf8cfcbd90af8d57a, .lo = 0x4221dc4ba772598d, .ex = -2, .sgn=0},
  {.hi = 0xc7c5c1e34d3055b2, .lo = 0x5cc8c00e4fccd850, .ex = -2, .sgn=0},
  {.hi = 0x964083747309d113, .lo = 0xa89a11e07c1fe, .ex = -2, .sgn=0},
  {.hi = 0xc8bd35e14da15f0e, .lo = 0xc7396c894bbf7389, .ex = -3, .sgn=0},
  {.hi = 0xc8fb2f886ec09f37, .lo = 0x6a17954b2b7c5171, .ex = -4, .sgn=0},
};

/* Table containing 128-bit approximations of sin(pi*i/2^12) for 0 <= i < 64
   (to nearest).
   Each entry is to be interpreted as (hi/2^64+lo/2^128)*2^ex*(-1)^sgn.
   Generated with computeS2() from sin.sage. */
static const dint64_t S2[64] = {
  {.hi = 0x0, .lo = 0x0, .ex = 0, .sgn=0},
  {.hi = 0xc90fd957659b030e, .lo = 0xc2afcacd698bdf76, .ex = -10, .sgn=0},
  {.hi = 0xc90fd57732396c23, .lo = 0x5a3af6ac41069419, .ex = -9, .sgn=0},
  {.hi = 0x96cbdb41258434c4, .lo = 0xa32c4560d01e7004, .ex = -8, .sgn=0},
  {.hi = 0xc90fc5f66525d257, .lo = 0x480f7956b6470765, .ex = -8, .sgn=0},
  {.hi = 0xfb53a8eb3ec38532, .lo = 0x8421cf8014d69789, .ex = -8, .sgn=0},
  {.hi = 0x96cbc117ccb5e27c, .lo = 0xb42bdad547b10410, .ex = -7, .sgn=0},
  {.hi = 0xafeda7e9ae465569, .lo = 0x185182c979282464, .ex = -7, .sgn=0},
  {.hi = 0xc90f87f3380388d5, .lo = 0xcb3ff35bd4d81baa, .ex = -7, .sgn=0},
  {.hi = 0xe231603c5e20db3e, .lo = 0xdb734d5b478687e2, .ex = -7, .sgn=0},
  {.hi = 0xfb532fcd151e2c42, .lo = 0xe76f0f818a303da0, .ex = -7, .sgn=0},
  {.hi = 0x8a3a7ad6a8e8b65f, .lo = 0x47af06ccf2bed25, .ex = -6, .sgn=0},
  {.hi = 0x96cb587284b81770, .lo = 0xb767005691b9d9d1, .ex = -6, .sgn=0},
  {.hi = 0xa35c303e18cc9b23, .lo = 0xa3585459daba5337, .ex = -6, .sgn=0},
  {.hi = 0xafed01bd602f0401, .lo = 0x4c816951623c2ac3, .ex = -6, .sgn=0},
  {.hi = 0xbc7dcc7456263d55, .lo = 0x5bfafd218b6a5bec, .ex = -6, .sgn=0},
  {.hi = 0xc90e8fe6f63c2330, .lo = 0xf1d7d06db39ea9fc, .ex = -6, .sgn=0},
  {.hi = 0xd59f4b993c424a6b, .lo = 0x628b4306a78ea2f5, .ex = -6, .sgn=0},
  {.hi = 0xe22fff0f2456c8a0, .lo = 0x30725e9973c8a45e, .ex = -6, .sgn=0},
  {.hi = 0xeec0a9ccaae8fc2a, .lo = 0x124efe701fddd951, .ex = -6, .sgn=0},
  {.hi = 0xfb514b55ccbe541a, .lo = 0xd784e031f9af76d6, .ex = -6, .sgn=0},
  {.hi = 0x83f0f197437b8c17, .lo = 0xfd74443f73a0b7b3, .ex = -5, .sgn=0},
  {.hi = 0x8a39386d6b899861, .lo = 0xda7803806001176f, .ex = -5, .sgn=0},
  {.hi = 0x908179ef5d7b775d, .lo = 0x2f0ae40f78b13b4f, .ex = -5, .sgn=0},
  {.hi = 0x96c9b5df1877e9b5, .lo = 0xf91ee371d6467dca, .ex = -5, .sgn=0},
  {.hi = 0x9d11ebfe9bdcac44, .lo = 0xc2f6cbeb4551eb34, .ex = -5, .sgn=0},
  {.hi = 0xa35a1c0fe740dbff, .lo = 0x3f13f6cb7c548aa5, .ex = -5, .sgn=0},
  {.hi = 0xa9a245d4fa7759e6, .lo = 0xadde30ec694a6399, .ex = -5, .sgn=0},
  {.hi = 0xafea690fd5912ef3, .lo = 0xf56e3c87ae3c56df, .ex = -5, .sgn=0},
  {.hi = 0xb632858278dff001, .lo = 0x53e382471fa0ba92, .ex = -5, .sgn=0},
  {.hi = 0xbc7a9aeee4f821b1, .lo = 0x94ad9b1aa972fe03, .ex = -5, .sgn=0},
  {.hi = 0xc2c2a9171ab39c54, .lo = 0xb13274ed67e9ffc7, .ex = -5, .sgn=0},
  {.hi = 0xc90aafbd1b33efc9, .lo = 0xc539edcbfda0cf2c, .ex = -5, .sgn=0},
  {.hi = 0xcf52aea2e7e4c75e, .lo = 0x3f87db6f41e2ead3, .ex = -5, .sgn=0},
  {.hi = 0xd59aa58a827e4daa, .lo = 0x370d90684597567c, .ex = -5, .sgn=0},
  {.hi = 0xdbe29435ed079069, .lo = 0xcd1c0c5d62efc5af, .ex = -5, .sgn=0},
  {.hi = 0xe22a7a6729d8e453, .lo = 0x850021e392744a4f, .ex = -5, .sgn=0},
  {.hi = 0xe87257e03b9e48eb, .lo = 0x7971fa8396212dd9, .ex = -5, .sgn=0},
  {.hi = 0xeeba2c632559cc53, .lo = 0x58418067afeb9868, .ex = -5, .sgn=0},
  {.hi = 0xf501f7b1ea65ef17, .lo = 0xca955048af12e9d, .ex = -5, .sgn=0},
  {.hi = 0xfb49b98e8e7807f6, .lo = 0xb21ccebc9caac3, .ex = -5, .sgn=0},
  {.hi = 0x80c8b8dd8ad153d4, .lo = 0x6f0804dae5f13b9b, .ex = -4, .sgn=0},
  {.hi = 0x83ec8ffcc22bfe51, .lo = 0xdb725856ff8c2846, .ex = -4, .sgn=0},
  {.hi = 0x87106205efb61b6a, .lo = 0x3f67ad05a5be69e5, .ex = -4, .sgn=0},
  {.hi = 0x8a342eda160bf5ae, .lo = 0xde5b1068d174be9c, .ex = -4, .sgn=0},
  {.hi = 0x8d57f65a37fd3c26, .lo = 0xaed92d34c17df538, .ex = -4, .sgn=0},
  {.hi = 0x907bb867588e3427, .lo = 0xc8c5732d89e910b0, .ex = -4, .sgn=0},
  {.hi = 0x939f74e27af8eb2e, .lo = 0xcc93966728e412d3, .ex = -4, .sgn=0},
  {.hi = 0x96c32baca2ae68b4, .lo = 0x37b2dd49d5fca3c0, .ex = -4, .sgn=0},
  {.hi = 0x99e6dca6d357dfff, .lo = 0x9a60c93317892c3b, .ex = -4, .sgn=0},
  {.hi = 0x9d0a87b210d7e1f8, .lo = 0xa318ba775fd7b039, .ex = -4, .sgn=0},
  {.hi = 0xa02e2caf5f4b8ef5, .lo = 0xf3d645e787b94f1d, .ex = -4, .sgn=0},
  {.hi = 0xa351cb7fc30bc889, .lo = 0xb56007d16d4ad5a3, .ex = -4, .sgn=0},
  {.hi = 0xa675640440ae634b, .lo = 0xdcd0d6bb4b9cd3e8, .ex = -4, .sgn=0},
  {.hi = 0xa998f61ddd0758a2, .lo = 0x17954ed6093c44c0, .ex = -4, .sgn=0},
  {.hi = 0xacbc81ad9d29f885, .lo = 0x5213c653bfb79b78, .ex = -4, .sgn=0},
  {.hi = 0xafe00694866a1b44, .lo = 0xcd34d2751c2e1da7, .ex = -4, .sgn=0},
  {.hi = 0xb30384b39e5d5346, .lo = 0xb7029d39efd76818, .ex = -4, .sgn=0},
  {.hi = 0xb626fbebeadc1ec6, .lo = 0x3a95642f565102f2, .ex = -4, .sgn=0},
  {.hi = 0xb94a6c1e7203198e, .lo = 0xfb8391d83d6da17d, .ex = -4, .sgn=0},
  {.hi = 0xbc6dd52c3a342eb5, .lo = 0xf10bfca3d6464012, .ex = -4, .sgn=0},
  {.hi = 0xbf9136f64a17ca4f, .lo = 0x9530f050886d7566, .ex = -4, .sgn=0},
  {.hi = 0xc2b4915da89e0b23, .lo = 0x5bfac0f965612e23, .ex = -4, .sgn=0},
  {.hi = 0xc5d7e4435cfff45c, .lo = 0x6718c1dfd2aa611c, .ex = -4, .sgn=0},
};

/* Table C1 is not needed, since cos(x) = sin(pi/2+x) for 0 <= x < pi/2,
   and cos(x) = -sin(x-pi/2) for pi/2 <= x < pi, thus
   C1[i] = S1[32+i] for 0 <= i < 32, and C1[i] = -S1[i-32] for 32 <= i < 64. */

/* Table containing 128-bit approximations of cos(pi*i/2^12) for 0 <= i < 64
   (to nearest).
   Each entry is to be interpreted as (hi/2^64+lo/2^128)*2^ex*(-1)^sgn.
   Generated with computeC2() from sin.sage. */
static const dint64_t C2[64] = {
  {.hi = 0x8000000000000000, .lo = 0x0, .ex = 1, .sgn=0},
  {.hi = 0xfffffb10b0d19f76, .lo = 0x491a703231e0a12e, .ex = 0, .sgn=0},
  {.hi = 0xffffec42c3773235, .lo = 0xec9e2e75cb525f2c, .ex = 0, .sgn=0},
  {.hi = 0xffffd3963882d553, .lo = 0x6276b75b91de105e, .ex = 0, .sgn=0},
  {.hi = 0xffffb10b10e80e95, .lo = 0x3031437d7eccb9df, .ex = 0, .sgn=0},
  {.hi = 0xffff84a14dfbcc6a, .lo = 0x858425d8b397dee6, .ex = 0, .sgn=0},
  {.hi = 0xffff4e58f17465de, .lo = 0x177338053fd93920, .ex = 0, .sgn=0},
  {.hi = 0xffff0e31fd699a85, .lo = 0x3a11d60588d8b96e, .ex = 0, .sgn=0},
  {.hi = 0xfffec42c7454926b, .lo = 0x38e310779edfec68, .ex = 0, .sgn=0},
  {.hi = 0xfffe7048590fddf8, .lo = 0xedd8e1034213f22b, .ex = 0, .sgn=0},
  {.hi = 0xfffe1285aed775d8, .lo = 0x96f351efc65556cf, .ex = 0, .sgn=0},
  {.hi = 0xfffdaae47948bad5, .lo = 0xea80aedd6a19710f, .ex = 0, .sgn=0},
  {.hi = 0xfffd3964bc6275ba, .lo = 0x69fff9ae0dedb047, .ex = 0, .sgn=0},
  {.hi = 0xfffcbe067c84d725, .lo = 0xf3a703b987eca44a, .ex = 0, .sgn=0},
  {.hi = 0xfffc38c9be717763, .lo = 0x928db07a0e70ba36, .ex = 0, .sgn=0},
  {.hi = 0xfffba9ae874b563a, .lo = 0x8d800bed6653dcba, .ex = 0, .sgn=0},
  {.hi = 0xfffb10b4dc96dabb, .lo = 0xb47903f7a19f8ee2, .ex = 0, .sgn=0},
  {.hi = 0xfffa6ddcc439d30a, .lo = 0xecc7b9244a48eb19, .ex = 0, .sgn=0},
  {.hi = 0xfff9c126447b7424, .lo = 0xfbe18032d0016082, .ex = 0, .sgn=0},
  {.hi = 0xfff90a91640459a1, .lo = 0x90e2d2eaf6da4d1e, .ex = 0, .sgn=0},
  {.hi = 0xfff84a1e29de8571, .lo = 0x8cc193c5d508e13f, .ex = 0, .sgn=0},
  {.hi = 0xfff77fcc9d755f99, .lo = 0x89332d07a713477e, .ex = 0, .sgn=0},
  {.hi = 0xfff6ab9cc695b5e8, .lo = 0x9e4938f661aa140c, .ex = 0, .sgn=0},
  {.hi = 0xfff5cd8ead6dbbab, .lo = 0x66c785e86dfbb75f, .ex = 0, .sgn=0},
  {.hi = 0xfff4e5a25a8d095b, .lo = 0x43366df666fd54ff, .ex = 0, .sgn=0},
  {.hi = 0xfff3f3d7d6e49c49, .lo = 0xdbb49f29fa872a83, .ex = 0, .sgn=0},
  {.hi = 0xfff2f82f2bc6d648, .lo = 0xe08b96133ecce0bd, .ex = 0, .sgn=0},
  {.hi = 0xfff1f2a862e77d4e, .lo = 0x98a31bcda3def20, .ex = 0, .sgn=0},
  {.hi = 0xfff0e343865bbb13, .lo = 0x5428ed0647c9e5d1, .ex = 0, .sgn=0},
  {.hi = 0xffefca00a09a1cb3, .lo = 0x807b6e7a4a723dae, .ex = 0, .sgn=0},
  {.hi = 0xffeea6dfbc7a9242, .lo = 0xccf344c647917821, .ex = 0, .sgn=0},
  {.hi = 0xffed79e0e5366e63, .lo = 0xf0f7cb05bde024f1, .ex = 0, .sgn=0},
  {.hi = 0xffec4304266865d9, .lo = 0x5657552366961732, .ex = 0, .sgn=0},
  {.hi = 0xffeb02498c0c8f12, .lo = 0x9195e99fbca2f10a, .ex = 0, .sgn=0},
  {.hi = 0xffe9b7b1228061b6, .lo = 0x191df31aaa6f7f45, .ex = 0, .sgn=0},
  {.hi = 0xffe8633af682b627, .lo = 0x3b57790bf77b6b2d, .ex = 0, .sgn=0},
  {.hi = 0xffe704e71533c508, .lo = 0x53aa9423bb0adc21, .ex = 0, .sgn=0},
  {.hi = 0xffe59cb58c1526b9, .lo = 0x3e71f7d99688082a, .ex = 0, .sgn=0},
  {.hi = 0xffe42aa66909d2d2, .lo = 0xbe28fbec7cfb8a6, .ex = 0, .sgn=0},
  {.hi = 0xffe2aeb9ba561f99, .lo = 0xf1ed54343fe7be24, .ex = 0, .sgn=0},
  {.hi = 0xffe128ef8e9fc17a, .lo = 0x7d209f32d42d864e, .ex = 0, .sgn=0},
  {.hi = 0xffdf9947f4edca6f, .lo = 0x8e6ee05573d420, .ex = 0, .sgn=0},
  {.hi = 0xffddffc2fca8a970, .lo = 0x44bd28b8d85b530a, .ex = 0, .sgn=0},
  {.hi = 0xffdc5c60b59a29dc, .lo = 0x75a8951fc304b914, .ex = 0, .sgn=0},
  {.hi = 0xffdaaf212fed72db, .lo = 0x4fd8f038449ec436, .ex = 0, .sgn=0},
  {.hi = 0xffd8f8047c2f06be, .lo = 0x8c9611f0b1d8f023, .ex = 0, .sgn=0},
  {.hi = 0xffd7370aab4cc25e, .lo = 0x8d3cd437dc7fa9d2, .ex = 0, .sgn=0},
  {.hi = 0xffd56c33ce95dc73, .lo = 0x45bd035edb0a65f5, .ex = 0, .sgn=0},
  {.hi = 0xffd3977ff7bae4e9, .lo = 0x664649b4d541b9c5, .ex = 0, .sgn=0},
  {.hi = 0xffd1b8ef38cdc433, .lo = 0xc42aac754bedcfde, .ex = 0, .sgn=0},
  {.hi = 0xffcfd081a441ba99, .lo = 0x1fd552bf146b4a6, .ex = 0, .sgn=0},
  {.hi = 0xffcdde374ceb5f7d, .lo = 0x76f487bb853a6989, .ex = 0, .sgn=0},
  {.hi = 0xffcbe2104600a0a9, .lo = 0x5595ca3f421ae09c, .ex = 0, .sgn=0},
  {.hi = 0xffc9dc0ca318c18b, .lo = 0x11b369083a7a62e5, .ex = 0, .sgn=0},
  {.hi = 0xffc7cc2c782c5a76, .lo = 0x5c2a6019679e41f, .ex = 0, .sgn=0},
  {.hi = 0xffc5b26fd99557dd, .lo = 0x579207cfe424dcb7, .ex = 0, .sgn=0},
  {.hi = 0xffc38ed6dc0ef98b, .lo = 0x1c676208aa3be545, .ex = 0, .sgn=0},
  {.hi = 0xffc1616194b5d1d3, .lo = 0xbc8d54e81d94f831, .ex = 0, .sgn=0},
  {.hi = 0xffbf2a101907c4c5, .lo = 0x965827f33d906c7c, .ex = 0, .sgn=0},
  {.hi = 0xffbce8e27ee40754, .lo = 0xe0aa07fcb29eef39, .ex = 0, .sgn=0},
  {.hi = 0xffba9dd8dc8b1e83, .lo = 0xccfed60a91097c48, .ex = 0, .sgn=0},
  {.hi = 0xffb848f3489ede86, .lo = 0xe907d9a298ab1feb, .ex = 0, .sgn=0},
  {.hi = 0xffb5ea31da2269e5, .lo = 0xbfdfce09aea7ac02, .ex = 0, .sgn=0},
  {.hi = 0xffb38194a87a3097, .lo = 0xbadfe70a1ef51116, .ex = 0, .sgn=0},
};


// accurate path for |x| >= 2^31
static double __attribute__((cold,noinline))
sin_large_accurate (double x)
{
  dint64_t r[1];
  uint64_t k = reduce_large_acc (r, x);
  /* x/(2*pi) mod 1 = k/2^13 + r + eps with |r| <= 2^-14 and 0 <= eps < 2^-127.999
     then sin(x) ~ sin(pi*k/2^12 + 2*pi*r)
                 ~ sin(pi*k/2^12)*cos(2*pi*r) + cos(pi*k/2^12)*sin(2*pi*r)
     Write k = 2^12*s + 2^6*i1 + i2 and t1=pi*i1/2^6, t2 = pi*i2/2^12
     then sin(pi*k/2^12) = (-1)^s*[sin(t1)*cos(t2)+cos(t1)*sin(t2)]
     and  cos(pi*k/2^12) = (-1)^s*[cos(t1)*cos(t2)-sin(t1)*sin(t2)]
  */
  int sbit = (x > 0) ? 0 : 1;
  sbit = sbit ^ (k >> (SHIFT-1));
  int i1 = (k >> 6) & 0x3f, i2 = k & 0x3f;
  dint64_t s1[1], s2[1], c1[1], c2[2];
  // approximate sin(t1)*cos(t2) in s1, and cos(t1)*sin(t2) in s2
  mul_dint (s1, S1+i1, C2+i2);
  mul_dint (s2, S1+((i1+32)&0x3f), S2+i2);
  s2->sgn ^= i1>=32;
  // add s1+s2
  add_dint (s1, s1, s2);
  // s1 approximates sin(t1+t2)
  // approximate cos(t1)*cos(t2) in c1, and sin(t1)*sin(t2) in c2
  mul_dint (c1, S1+((i1+32)&0x3f), C2+i2);
  c1->sgn ^= i1>=32;
  mul_dint (c2, S1+i1, S2+i2);
  // compute c1-c2
  c2->sgn ^= 1;
  add_dint (c1, c1, c2);
  // c1 approximates cos(t1+t2)

  dint64_t Sr[1], Cr[1], r2[1];
  mul_dint (r2, r, r); // r2 approximates r^2
  evalPS (Sr, r, r2); // Sr approximates sin(2*pi*r)
  evalPC (Cr, r2);    // Cr approximates cos(2*pi*r)

  // now combine: sin(x) ~ s1*C + c1*S
  mul_dint (s1, s1, Cr);
  mul_dint (c1, c1, Sr);
  add_dint (s1, s1, c1);
  s1->sgn ^= sbit;
  return dint_tod (s1);
}

// fast path for |x| >= 2^31
// ax = |x| and eps is the error bound for the rounding test
static double __attribute__((noinline))
cr_sin_large (double x)
{
  double ax = __builtin_fabs(x);
  double r;
  uint64_t j = reduce_large (&r, ax);
  // now x/(2pi) ~ k + j/2^15 + r with 0 <= r < 2^-15
  int sbit = (x > 0) ? 0 : 1;

  double r2 = r * r;
  sbit = sbit ^ (j >> 14); // reduction by an odd multiple of pi?
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
  static const double Sgn[] = {1.0, -1.0};
  fh = Sgn[sbit] * fh;
  fl = Sgn[sbit] * fl;
  static const double eps = 0x1.41p-63;
  double lb = fh + (fl - eps), ub = fh + (fl + eps);
  // if (__builtin_expect (lb == ub, 1)) return lb;
  return sin_large_accurate (x);
}

double
cr_sin (double x)
{
  b64u64_u t = {.f = x};
  // deal with tiny x to avoid underflow
  uint64_t au = t.u<<1;
  if (__builtin_expect(au <= 0x7cae26e892247decull, 0)) {
    // |x| <= 0x1.7137449123ef6p-26
    if (au == 0) return x;
    // Taylor expansion of sin(x) is x - x^3/6 around zero
    // for x=-0, fma (x, -0x1p-54, x) returns +0
    /* We have underflow when 0 < |x| < 2^-1022 or when |x| = 2^-1022
       and rounding towards zero. */
    double res = __builtin_fma (x, -0x1p-54, x);
#ifdef CORE_MATH_SUPPORT_ERRNO
    if ((t.u<<1)<(1ull<<53) || __builtin_fabs (res) < 0x1p-1022)
      errno = ERANGE; // underflow
#endif
    return res;
  }
  int e = (t.u>>52)&0x7ff;
  if (__builtin_expect(e < 1054, 1)) return cr_sin_moderate(x, t.u>>63); // |x| < 2^31
  if (__builtin_expect (e == 0x7ff, 0)) /* NaN, +Inf and -Inf. */
    {
#ifdef CORE_MATH_SUPPORT_ERRNO
      if ((t.u<<1) == 0x7ffull<<53) // +/-Inf
        errno = EDOM;
#endif
      return x - x; // raises invalid
    }
  // now |x| >= 2^31
  return cr_sin_large (x);
}
