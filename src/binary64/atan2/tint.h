/* Fast 192-bit arithmetic routines.

Copyright (c) 2023 Paul Zimmermann

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

#include <assert.h>

typedef unsigned __int128 u128;

// the following represent (-1)^sgn*(h/2^64+m/2^128+l/2^192)*2^ex
// we have either h=m=l=0 to represent +0 or -0
// or the most significant bit of h is 1
typedef union {
  struct {
    uint64_t m, h, l; // put m before h on little-endian processor
    int64_t ex;
    uint64_t sgn;
  };
  struct {
    u128 _h;
    uint64_t _l;
    int64_t _ex;
    uint64_t _sgn;
  };
} tint_t;

// ZERO is a tint_t representation of 0
static const tint_t ZERO = {.h = 0, .m = 0, .l = 0, .ex = -1076, .sgn = 0};

// ONE is a tint_t representation of 1
static const tint_t ONE = {
  .h = 0x8000000000000000, .m = 0, .l = 0, .ex = 1, .sgn = 0};

// PI2 is a tint_t representation of pi/2
static const tint_t PI2 = {
  .h = 0xc90fdaa22168c234, .m = 0xc4c6628b80dc1cd1, .l = 0x29024e088a67cc74,
  .ex = 1, .sgn = 0};

// Print a tint_t value for debugging purposes
static inline void print_tint (const tint_t *a) {
  printf("{.h=0x%lx, .m=0x%lx, .l=0x%lx, .ex=%ld, .sgn=0x%lx}\n",
         a->h, a->m, a->l, a->ex, a->sgn);
}
// Copy a tint_t value
static inline void cp_tint(tint_t *r, const tint_t *a) {
  r->_h = a->_h;
  r->_l = a->_l;
  r->_ex = a->_ex;
  r->_sgn = a->_sgn;
}

static inline int
is_normalized (const tint_t *a)
{
  if (a->h == 0 && a->m == 0 && a->l == 0)
    return 1;
  return a->h >> 63;
}

// Multiply two tint_t numbers, with error < 10 ulps
// Overlap between r and a or b is allowed
static inline void
mul_tint (tint_t *r, const tint_t *a, const tint_t *b)
{
  // printf ("enter mul_tint, a="); print_tint (a); printf ("b="); print_tint (b);
  assert (is_normalized (a));
  assert (is_normalized (b));
  u128 ah = a->h, am = a->m, al = a->l;
  u128 bh = b->h, bm = b->m, bl = b->l;
  // printf ("enter mul_tint: ah=%lx bh=%lx\n", (uint64_t) ah, (uint64_t) bh);
  u128 rh = ah * bh, rm1 = ah * bm, rm2 = am * bh;
  u128 rl1 = ah * bl, rl2 = am * bm, rl3 = al * bh;
  uint64_t h, l, cm;
  r->h = rh >> 64;
  r->m = rh; // cast to low 64 bits 
  // accumulate rm1
  r->l = rm1; // cast to low 64 bits
  h = rm1 >> 64;
  r->m += h;
  rh += r->m < h;
  // accumulate rm2
  l = rm2; // cast to low 64 bits
  r->l += l;
  cm = r->l < l;
  h = rm2 >> 64;
  r->m += h;
  rh += r->m < h;
  // accumulate rl1+rl2+rl3
  rl1 = (rl1 >> 64) + (rl2 >> 64) + (rl3 >> 64);
  l = rl1; // cast to low 64 bits
  r->l += l;
  cm += r->l < l;
  // accumulate cm
  r->m += cm;
  r->h += r->m < cm;

  r->ex = a->ex + b->ex;
  r->sgn = a->sgn ^ b->sgn;

  /* Note: if one of the operands was zero, then r->h = r->m = r->l = 0,
     and the normalization keeps r=0. */
  if (!(r->h >> 63)) // normalize
  {
    r->h = (r->h << 1) | (r->m >> 63);
    r->m = (r->m << 1) | (r->l >> 63);
    r->l = r->l << 1;
    r->ex --;
  }

  // printf ("exit mul_tint, r="); print_tint (r);
  assert (is_normalized (r));

  /* We ignored the following terms, denoting B=2^64, related to r->l:
     am*bl + al*bm <= 2*(B-1)^2/B^2 = 2 - 4/B + 2/B^2
     al*bl <= (B-1)^2/B^3 = 1/B - 2/B^2 + 1/B^3
     And we truncated rl1+rl2+rl3:
     low(rl1) + low(rl2) + low(rl3) <= 3*(B-1)/B = 3-3/B
     This sums up to: 5 - 6/B + 1/B^3 < 5
     thus the rounding error is bounded by 5 ulps,
     and after normalization by 10 ulps. */
}

// Return non-zero if a = 0
static inline int
tint_zero_p (const tint_t *a)
{
  return a->h == 0;
}

static inline int cmp(int64_t a, int64_t b) { return (a > b) - (a < b); }
static inline char cmpu64(uint64_t a, uint64_t b) { return (a > b) - (a < b); }
static inline char cmpu128(u128 a, u128 b) { return (a > b) - (a < b); }

// Compare the absolute values of a and b
// Return -1 if |a| < |b|
// Return  0 if |a| = |b|
// Return +1 if |a| > |b|
static inline int
cmp_tint_abs (const tint_t *a, const tint_t *b) {
  if (tint_zero_p (a))
    return tint_zero_p (b) ? 0 : -1;
  if (tint_zero_p (b))
    return +1;
  int c = cmp (a->ex, b->ex);
  if (c)
    return c;
  // now a->ex = b->ex
  c = cmpu128 (a->_h, b->_h);
  if (c)
    return c;
  return cmpu64 (a->_l, b->_l);
}

// shift right by k bits (only deal with the significand)
static inline void
rshift (tint_t *a, const tint_t *b, int k)
{
  if (k == 0)
  {
    a->_h = b->_h;
    a->_l = b->_l;
  }
  else if (k < 64)
  {
    a->_h = b->_h >> k;
    a->_l = (b->_h << (64 - k)) | (b->_l >> k);
  }
  else if (k == 64)
  {
    a->_h = b->_h >> k;
    a->_l = b->_h;
  }
  else if (k < 128)
  {
    a->_h = b->_h >> k;
    a->_l = b->_h >> (k - 64);
  }
  else if (k < 192)
  {
    a->_h = 0;
    a->_l = b->_h >> (k - 128);
  }
  else
    a->_h = a->_l = 0;
  // printf ("exit rshift a="); print_tint (a);
}

// shift left by k bits (only deal with the significand)
static inline void
lshift (tint_t *a, const tint_t *b, int k)
{
  if (k == 0)
  {
    a->_h = b->_h;
    a->_l = b->_l;
  }
  else if (k < 64)
  {
    a->_h = (b->_h << k) | (b->_l >> (64 - k));
    a->_l = b->_l << k;
  }
  else if (k == 64)
  {
    a->_h = (b->_h << k) | (u128) b->_l;
    a->_l = 0;
  }
  else if (k < 128)
  {
    a->_h = b->_h << k | ((u128) b->_l << (k - 64));
    a->_l = 0;
  }
  else if (k < 192)
  {
    a->_h = (u128) b->_l << (k - 64);
    a->_l = 0;
  }
  else
    a->_h = a->_l = 0;
}

// Add two tint_t values
static inline void
add_tint (tint_t *r, const tint_t *a, const tint_t *b)
{
  // printf ("enter add_tint, a="); print_tint (a); printf ("b="); print_tint (b);
  assert (is_normalized (a));
  assert (is_normalized (b));
  switch (cmp_tint_abs (a, b))
  {
  case 0: // |a| = |b|
    if (a->sgn ^ b->sgn) {
      cp_tint (r, &ZERO);
      goto end;
    }
    cp_tint (r, a);
    r->ex++;
    goto end;

  case -1: // |a| < |b|
    {
      // swap operands
      const tint_t *tmp = a; a = b; b = tmp;
      break; // fall through the case |a| > |b|
    }
  }

  // From now on, |a| > |b| thus a->ex >= b->ex
  tint_t t[1];
  rshift (t, b, a->ex - b->ex);

  if (a->sgn ^ b->sgn) { // opposite signs, it's a subtraction
    t->_l = a->_l - t->_l;
    t->_h = a->_h - t->_h - (t->_l > a->_l);
    uint64_t th = t->_h >> 64;
    uint64_t ex =
      th ? __builtin_clzl (th)
      : (t->_h ? 64 + __builtin_clzl (t->_h) : 128 + __builtin_clzl (t->_l));
    r->ex = a->ex - ex;
    lshift (r, t, ex);
  }
  else { // same signs, it's an addition
    r->_l = a->_l + t->_l;
    uint64_t cl = t->_l < a->_l;
    r->_h = a->_h + t->_h;
    uint64_t ch = r->_h < a->_h;
    r->_h += cl;
    ch += r->_h < cl;
    if (ch) { // can be at most 1
      r->ex = a->ex + 1;
      r->_l = (r->_h << 127) | (r->_l >> 1);
      r->_h = ((u128) ch << 127) | (r->_h >> 1);
    }
    else
      r->ex = a->ex;
  }
  r->sgn = a->sgn;
 end:
  // printf ("exit add_tint, r="); print_tint (r);
  assert (is_normalized (r));
}

// a <- x, assuming x is not NaN, Inf or 0
static inline void tint_fromd (tint_t *a, double x)
{
  d64u64 u = {.f = x};
  a->sgn = u.u >> 63;
  uint64_t ax = u.u & 0x7ffffffffffffffful;
  int64_t e = ax >> 52;
  if (__builtin_expect (e, 1)) { // normal
    // 1 has e=0x3ff
    a->ex = e - 0x3fe;
    a->h = (1ul << 63) | (ax << 11);
  }
  else { // subnormal
    // 2^-1074 has ax=1
    e = __builtin_clzl (ax);
    a->ex = -0x3f2 - e;
    a->h = ax << e;
  }
  a->m = a->l = 0;
}

static inline double tint_tod (tint_t *a)
{
  if (a->ex >= 1025) // overflow: |a| >= 2^1024
    return a->sgn ? -0x1p1023 - 0x1p1023 : 0x1p1023 + 0x1p1023;
  if (a->ex <= -1074) // underflow: |a| < 2^-1074
  {
    if (a->ex < -1074) // |a| < 2^-1075
      return (a->sgn ? -0x1p-1074 : 0x1p-1074) * 0.5;
    // 2^-1075 <= |a| < 2^-1074
    int mid = a->h == (1ul << 63) && a->m == 0 && a->l == 0;
    // if mid, |a| = 2^-1075
    return (a->sgn ? -0x1p-1074 : 0x1p-1074) * (mid ? 0.5 : 0.75);
  }
#define MASK53 0x1ffffffffffffful
  double r[4];
  r[3] = a->h >> 11; // 53 bits from a->h
  r[2] = ((a->h << 42) & MASK53) | (a->m >> 22); // a->h:11 bits, a->m:42 bits
  r[1] = ((a->m << 31) & MASK53) | (a->l >> 33); // a->m:22 bits, a->l:31 bits
  r[0] = (a->l << 20) & MASK53; // 33 bits from a->l
  static const double S[2] = {1.0, -1.0};
  double s = S[a->sgn];
  r[1] = __builtin_fma (s * r[0], 0x1p-53, s * r[1]);
  r[2] = __builtin_fma (r[1], 0x1p-53, s * r[2]);
  r[3] = __builtin_fma (r[2], 0x1p-53, s * r[3]);
  r[3] *= 0x1p-53;
  return r[3] * __builtin_ldexp (1.0, a->ex);
}

// put in r a 106-bit approximation of 1/a, assuming a is not zero
// assume A = tint_fromd (a)
static inline void inv_tint (tint_t *r, const tint_t *A, double a)
{
  // printf ("enter inv_tint a=%la A=", a); print_tint (A);
  tint_t q[1];
  tint_fromd (r, 1.0 / a); // accurate to about 53 bits
  // we use Newton's iteration: r -> r + r*(1-a*r)
  mul_tint (q, A, r);      // a*r
  q->sgn = 1 - q->sgn;     // -a*r
  add_tint (q, &ONE, q);   // 1-a*r
  mul_tint (q, r, q);      // r*(1-a*r)
  add_tint (r, r, q);
  // printf ("#### exit inv_tint r="); print_tint (r);
}

// put in r an approximation of b/a, assuming a is not zero
static inline void div_tint (tint_t *r, double b, double a)
{
  tint_t A[1], B[1], Y[1], Z[1];
  tint_fromd (A, a);
  tint_fromd (B, b);
  inv_tint (Y, A, a); // Y = 1/a to about 106 bits
  mul_tint (r, Y, B); // r = b/a to about 106 bits
  // we use Karp-Markstein's trick: r' = r + y*(b-a*r)
  mul_tint (Z, A, r);  // a*r
  Z->sgn = 1 - Z->sgn; // -a*r
  add_tint (Z, B, Z);  // b-a*r
  mul_tint (Z, Y, Z);  // y*(b-a*r)
  add_tint (r, r, Z);
}


