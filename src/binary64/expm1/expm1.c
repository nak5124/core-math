/* Correctly-rounded expm1 function for binary64 value.

Copyright (c) 2022-2023 Paul Zimmermann, Tom Hubrecht and Claude-Pierre Jeannerod

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
#include <stdint.h>
#include <math.h>

/****************** code copied from pow.[ch] ********************************/

typedef union {
  double f;
  uint64_t u;
} f64_u;

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma (a, b, -*hi);
}

// Multiply a double with a double double : a * (bh + bl)
static inline void s_mul (double *hi, double *lo, double a, double bh,
                          double bl) {
  a_mul (hi, lo, a, bh); /* exact */
  *lo = __builtin_fma (a, bl, *lo);
}

// Returns (ah + al) * (bh + bl) - (al * bl)
static inline void d_mul(double *hi, double *lo, double ah, double al,
                         double bh, double bl) {
  a_mul (hi, lo, ah, bh);
  *lo = __builtin_fma (ah, bl, *lo);
  *lo = __builtin_fma (al, bh, *lo);
}

// Add a + b, assuming |a| >= |b|
static inline void
fast_two_sum(double *hi, double *lo, double a, double b)
{
  double e;

  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
}

// Add a + (bh + bl), assuming |a| >= |bh|
static inline void fast_sum(double *hi, double *lo, double a, double bh,
                            double bl) {
  fast_two_sum(hi, lo, a, bh);
  /* |(a+bh)-(hi+lo)| <= 2^-105 |hi| and |lo| < ulp(hi) */
  *lo += bl;
  /* |(a+bh+bl)-(hi+lo)| <= 2^-105 |hi| + ulp(lo),
     where |lo| <= ulp(hi) + |bl|. */
}

/* For 0 <= i < 64, T1[i] = (h,l) such that h+l is the best double-double
   approximation of 2^(i/64). The approximation error is bounded as follows:
   |h + l - 2^(i/64)| < 2^-107. */
static const double T1[][2] = {
    {              0x1p+0,                 0x0p+0},
    {0x1.02c9a3e778061p+0, -0x1.19083535b085dp-56},
    {0x1.059b0d3158574p+0,  0x1.d73e2a475b465p-55},
    {0x1.0874518759bc8p+0,  0x1.186be4bb284ffp-57},
    {0x1.0b5586cf9890fp+0,  0x1.8a62e4adc610bp-54},
    {0x1.0e3ec32d3d1a2p+0,  0x1.03a1727c57b53p-59},
    {0x1.11301d0125b51p+0, -0x1.6c51039449b3ap-54},
    { 0x1.1429aaea92dep+0, -0x1.32fbf9af1369ep-54},
    {0x1.172b83c7d517bp+0, -0x1.19041b9d78a76p-55},
    {0x1.1a35beb6fcb75p+0,  0x1.e5b4c7b4968e4p-55},
    {0x1.1d4873168b9aap+0,  0x1.e016e00a2643cp-54},
    {0x1.2063b88628cd6p+0,  0x1.dc775814a8495p-55},
    {0x1.2387a6e756238p+0,  0x1.9b07eb6c70573p-54},
    {0x1.26b4565e27cddp+0,  0x1.2bd339940e9d9p-55},
    {0x1.29e9df51fdee1p+0,  0x1.612e8afad1255p-55},
    {0x1.2d285a6e4030bp+0,  0x1.0024754db41d5p-54},
    {0x1.306fe0a31b715p+0,  0x1.6f46ad23182e4p-55},
    {0x1.33c08b26416ffp+0,  0x1.32721843659a6p-54},
    {0x1.371a7373aa9cbp+0, -0x1.63aeabf42eae2p-54},
    {0x1.3a7db34e59ff7p+0, -0x1.5e436d661f5e3p-56},
    {0x1.3dea64c123422p+0,  0x1.ada0911f09ebcp-55},
    {0x1.4160a21f72e2ap+0, -0x1.ef3691c309278p-58},
    {0x1.44e086061892dp+0,   0x1.89b7a04ef80dp-59},
    { 0x1.486a2b5c13cdp+0,   0x1.3c1a3b69062fp-56},
    {0x1.4bfdad5362a27p+0,  0x1.d4397afec42e2p-56},
    {0x1.4f9b2769d2ca7p+0, -0x1.4b309d25957e3p-54},
    {0x1.5342b569d4f82p+0, -0x1.07abe1db13cadp-55},
    {0x1.56f4736b527dap+0,  0x1.9bb2c011d93adp-54},
    {0x1.5ab07dd485429p+0,  0x1.6324c054647adp-54},
    {0x1.5e76f15ad2148p+0,  0x1.ba6f93080e65ep-54},
    {0x1.6247eb03a5585p+0, -0x1.383c17e40b497p-54},
    {0x1.6623882552225p+0, -0x1.bb60987591c34p-54},
    {0x1.6a09e667f3bcdp+0, -0x1.bdd3413b26456p-54},
    {0x1.6dfb23c651a2fp+0, -0x1.bbe3a683c88abp-57},
    {0x1.71f75e8ec5f74p+0, -0x1.16e4786887a99p-55},
    {0x1.75feb564267c9p+0, -0x1.0245957316dd3p-54},
    {0x1.7a11473eb0187p+0, -0x1.41577ee04992fp-55},
    {0x1.7e2f336cf4e62p+0,  0x1.05d02ba15797ep-56},
    {0x1.82589994cce13p+0, -0x1.d4c1dd41532d8p-54},
    {0x1.868d99b4492edp+0, -0x1.fc6f89bd4f6bap-54},
    {0x1.8ace5422aa0dbp+0,  0x1.6e9f156864b27p-54},
    {0x1.8f1ae99157736p+0,  0x1.5cc13a2e3976cp-55},
    {0x1.93737b0cdc5e5p+0, -0x1.75fc781b57ebcp-57},
    { 0x1.97d829fde4e5p+0, -0x1.d185b7c1b85d1p-54},
    { 0x1.9c49182a3f09p+0,  0x1.c7c46b071f2bep-56},
    {0x1.a0c667b5de565p+0, -0x1.359495d1cd533p-54},
    {0x1.a5503b23e255dp+0, -0x1.d2f6edb8d41e1p-54},
    {0x1.a9e6b5579fdbfp+0,  0x1.0fac90ef7fd31p-54},
    {0x1.ae89f995ad3adp+0,  0x1.7a1cd345dcc81p-54},
    {0x1.b33a2b84f15fbp+0, -0x1.2805e3084d708p-57},
    {0x1.b7f76f2fb5e47p+0, -0x1.5584f7e54ac3bp-56},
    {0x1.bcc1e904bc1d2p+0,  0x1.23dd07a2d9e84p-55},
    {0x1.c199bdd85529cp+0,  0x1.11065895048ddp-55},
    {0x1.c67f12e57d14bp+0,  0x1.2884dff483cadp-54},
    {0x1.cb720dcef9069p+0,  0x1.503cbd1e949dbp-56},
    {0x1.d072d4a07897cp+0, -0x1.cbc3743797a9cp-54},
    {0x1.d5818dcfba487p+0,  0x1.2ed02d75b3707p-55},
    {0x1.da9e603db3285p+0,  0x1.c2300696db532p-54},
    {0x1.dfc97337b9b5fp+0, -0x1.1a5cd4f184b5cp-54},
    {0x1.e502ee78b3ff6p+0,  0x1.39e8980a9cc8fp-55},
    {0x1.ea4afa2a490dap+0, -0x1.e9c23179c2893p-54},
    {0x1.efa1bee615a27p+0,   0x1.dc7f486a4b6bp-54},
    { 0x1.f50765b6e454p+0,  0x1.9d3e12dd8a18bp-54},
    {0x1.fa7c1819e90d8p+0,  0x1.74853f3a5931ep-55},
};

/* For 0 <= i < 64, T2[i] = (h,l) such that h+l is the best double-double
   approximation of 2^(i/2^12). The approximation error is bounded as follows:
   |h + l - 2^(i/2^12)| < 2^-107. */
static const double T2[][2] = {
    {              0x1p+0,                 0x0p+0},
    {0x1.000b175effdc7p+0,  0x1.ae8e38c59c72ap-54},
    {0x1.00162f3904052p+0, -0x1.7b5d0d58ea8f4p-58},
    {0x1.0021478e11ce6p+0,  0x1.4115cb6b16a8ep-54},
    {0x1.002c605e2e8cfp+0, -0x1.d7c96f201bb2fp-55},
    {0x1.003779a95f959p+0,  0x1.84711d4c35e9fp-54},
    {0x1.0042936faa3d8p+0, -0x1.0484245243777p-55},
    { 0x1.004dadb113dap+0, -0x1.4b237da2025f9p-54},
    {0x1.0058c86da1c0ap+0, -0x1.5e00e62d6b30dp-56},
    {0x1.0063e3a559473p+0,  0x1.a1d6cedbb9481p-54},
    {0x1.006eff583fc3dp+0, -0x1.4acf197a00142p-54},
    {0x1.007a1b865a8cap+0, -0x1.eaf2ea42391a5p-57},
    {0x1.0085382faef83p+0,  0x1.da93f90835f75p-56},
    {0x1.00905554425d4p+0, -0x1.6a79084ab093cp-55},
    {0x1.009b72f41a12bp+0,  0x1.86364f8fbe8f8p-54},
    {0x1.00a6910f3b6fdp+0, -0x1.82e8e14e3110ep-55},
    {0x1.00b1afa5abcbfp+0, -0x1.4f6b2a7609f71p-55},
    {0x1.00bcceb7707ecp+0, -0x1.e1a258ea8f71bp-56},
    {0x1.00c7ee448ee02p+0,  0x1.4362ca5bc26f1p-56},
    {0x1.00d30e4d0c483p+0,  0x1.095a56c919d02p-54},
    {0x1.00de2ed0ee0f5p+0, -0x1.406ac4e81a645p-57},
    { 0x1.00e94fd0398ep+0,  0x1.b5a6902767e09p-54},
    {0x1.00f4714af41d3p+0, -0x1.91b2060859321p-54},
    {0x1.00ff93412315cp+0,  0x1.427068ab22306p-55},
    {0x1.010ab5b2cbd11p+0,  0x1.c1d0660524e08p-54},
    {0x1.0115d89ff3a8bp+0, -0x1.e7bdfb3204be8p-54},
    {0x1.0120fc089ff63p+0,  0x1.843aa8b9cbbc6p-55},
    {0x1.012c1fecd613bp+0, -0x1.34104ee7edae9p-56},
    {0x1.0137444c9b5b5p+0, -0x1.2b6aeb6176892p-56},
    {0x1.01426927f5278p+0,  0x1.a8cd33b8a1bb3p-56},
    {0x1.014d8e7ee8d2fp+0,  0x1.2edc08e5da99ap-56},
    {0x1.0158b4517bb88p+0,  0x1.57ba2dc7e0c73p-55},
    {0x1.0163da9fb3335p+0,  0x1.b61299ab8cdb7p-54},
    {0x1.016f0169949edp+0, -0x1.90565902c5f44p-54},
    {0x1.017a28af25567p+0,  0x1.70fc41c5c2d53p-55},
    {0x1.018550706ab62p+0,  0x1.4b9a6e145d76cp-54},
    {0x1.019078ad6a19fp+0, -0x1.008eff5142bf9p-56},
    {0x1.019ba16628de2p+0, -0x1.77669f033c7dep-54},
    {0x1.01a6ca9aac5f3p+0, -0x1.09bb78eeead0ap-54},
    {0x1.01b1f44af9f9ep+0,  0x1.371231477ece5p-54},
    {0x1.01bd1e77170b4p+0,  0x1.5e7626621eb5bp-56},
    {0x1.01c8491f08f08p+0, -0x1.bc72b100828a5p-54},
    { 0x1.01d37442d507p+0, -0x1.ce39cbbab8bbep-57},
    {0x1.01de9fe280ac8p+0,  0x1.16996709da2e2p-55},
    {0x1.01e9cbfe113efp+0, -0x1.c11f5239bf535p-55},
    {0x1.01f4f8958c1c6p+0,  0x1.e1d4eb5edc6b3p-55},
    {0x1.020025a8f6a35p+0, -0x1.afb99946ee3fp-54},
    {0x1.020b533856324p+0, -0x1.8f06d8a148a32p-54},
    {0x1.02168143b0281p+0, -0x1.2bf310fc54eb6p-55},
    {0x1.0221afcb09e3ep+0, -0x1.c95a035eb4175p-54},
    {0x1.022cdece68c4fp+0, -0x1.491793e46834dp-54},
    {0x1.02380e4dd22adp+0, -0x1.3e8d0d9c49091p-56},
    {0x1.02433e494b755p+0, -0x1.314aa16278aa3p-54},
    {0x1.024e6ec0da046p+0,  0x1.48daf888e9651p-55},
    {0x1.02599fb483385p+0,  0x1.56dc8046821f4p-55},
    {0x1.0264d1244c719p+0,  0x1.45b42356b9d47p-54},
    {0x1.027003103b10ep+0, -0x1.082ef51b61d7ep-56},
    {0x1.027b357854772p+0,  0x1.2106ed0920a34p-56},
    {0x1.0286685c9e059p+0, -0x1.fd4cf26ea5d0fp-54},
    {0x1.02919bbd1d1d8p+0, -0x1.09f8775e78084p-54},
    {0x1.029ccf99d720ap+0,  0x1.64cbba902ca27p-58},
    {0x1.02a803f2d170dp+0,  0x1.4383ef231d207p-54},
    {0x1.02b338c811703p+0,  0x1.4a47a505b3a47p-54},
    {0x1.02be6e199c811p+0,  0x1.e47120223467fp-54},
};

/* The following is a degree-4 polynomial generated by Sollya for exp(x)
   over [-0.000130273,0.000130273] with absolute error < 2^-74.346. */
static const double Q_1[] = {0x1p0,                 /* degree 0 */
                             0x1p0,                 /* degree 1 */
                             0x1p-1,                /* degree 2 */
                             0x1.5555555995d37p-3,  /* degree 3 */
                             0x1.55555558489dcp-5   /* degree 4 */
};

// Approximation for the fast path of exp(z) for z=zh+zl,
// with |z| < 0.000130273 < 2^-12.88 and |zl| < 2^-42.6
// (assuming x^y does not overflow or underflow)
static inline void q_1 (double *hi, double *lo, double zh, double zl) {
  double z = zh + zl;
  double q = __builtin_fma (Q_1[4], zh, Q_1[3]);

  q = __builtin_fma (q, z, Q_1[2]);

  fast_two_sum (hi, lo, Q_1[1], q * z);

  d_mul (hi, lo, zh, zl, *hi, *lo);

  fast_sum (hi, lo, Q_1[0], *hi, *lo);
}

/* __builtin_roundeven was introduced in gcc 10:
   https://gcc.gnu.org/gcc-10/changes.html,
   and in clang 17 */
#if (defined(__GNUC__) && __GNUC__ >= 10) || (defined(__clang__) && __clang_major__ >= 17)
#define HAS_BUILTIN_ROUNDEVEN
#endif

#if !defined(HAS_BUILTIN_ROUNDEVEN) && (defined(__GNUC__) || defined(__clang__)) && (defined(__AVX__) || defined(__SSE4_1__))
inline double __builtin_roundeven(double x){
   double ix;
#if defined __AVX__
   __asm__("vroundsd $0x8,%1,%1,%0":"=x"(ix):"x"(x));
#else /* __SSE4_1__ */
   __asm__("roundsd $0x8,%1,%0":"=x"(ix):"x"(x));
#endif
   return ix;
}
#define HAS_BUILTIN_ROUNDEVEN
#endif

#ifndef HAS_BUILTIN_ROUNDEVEN
#include <math.h>
/* round x to nearest integer, breaking ties to even */
static double
__builtin_roundeven (double x)
{
  double y = round (x); /* nearest, away from 0 */
  if (fabs (y - x) == 0.5)
  {
    /* if y is odd, we should return y-1 if x>0, and y+1 if x<0 */
    union { double f; uint64_t n; } u, v;
    u.f = y;
    v.f = (x > 0) ? y - 1.0 : y + 1.0;
    if (__builtin_ctz (v.n) > __builtin_ctz (u.n))
      y = v.f;
  }
  return y;
}
#endif

/*
  Approximation of exp(x)
  (the code in pow.c has x = xh + xl as input, we simplified it since here
  we have xl=0, we kept the corresponding line in comment for reference).

  exp(x) is approximated by hi + lo.

  For the error analysis, we only consider the case where x^y does not
  overflow or underflow. We get:

  (hi + lo) / exp(x) = 1 + eps with |eps| < 2^-74.139

  At output, we also have 0.99985 < hi+lo < 1.99995 and |lo/hi| < 2^-41.4.
*/

// static inline void exp_1 (double *hi, double *lo, double xh double xl) {
static inline void exp_1 (double *hi, double *lo, double x) {

#define INVLOG2 0x1.71547652b82fep+12 /* |INVLOG2-2^12/log(2)| < 2^-43.4 */
  // double k = __builtin_roundeven (xh * INVLOG2);
  double k = __builtin_roundeven (x * INVLOG2);

  double kh, kl;
#define LOG2H 0x1.62e42fefa39efp-13
#define LOG2L 0x1.abc9e3b39803fp-68
  s_mul (&kh, &kl, k, LOG2H, LOG2L);

  double yh, yl;
  // fast_two_sum (&yh, &yl, xh - kh, xl);
  yh = x - kh;
  // yl -= kl;
  yl = -kl;

  int64_t K = k; /* Note: k is an integer, this is just a conversion. */
  int64_t M = (K >> 12) + 0x3ff;
  int64_t i2 = (K >> 6) & 0x3f;
  int64_t i1 = K & 0x3f;

  double t1h = T1[i2][0], t1l = T1[i2][1], t2h = T2[i1][0], t2l = T2[i1][1];
  d_mul (hi, lo, t2h, t2l, t1h, t1l);

  double qh, ql;
  q_1 (&qh, &ql, yh, yl);

  d_mul (hi, lo, *hi, *lo, qh, ql);
  f64_u _d;

  _d.u = M << 52;
  *hi *= _d.f;
  *lo *= _d.f;
}

/****************** end of code copied from pow.[ch] *************************/

/* The following is a degree-7 polynomial generated by Sollya for exp(z)
   over [-0.000130273,0.000130273] with absolute error < 2^-113.218
   (see file exp_accurate.sollya). Since we use this code only for
   |x| > 0.125 in expm1(x), the corresponding relative error for expm1
   is about 2^-113.218/|expm1(-0.125)| which is about 2^-110. */
static const double Q_2[] = {
  0x1p0, // degree 0, Q_2[0]
  0x1p0, // degree 1, Q_2[1]
  0x1p-1, // degree 2, Q_2[2]
  0x1.5555555555555p-3, 0x1.55555555c4d26p-57, // degree 3, Q_2[3], Q_2[4]
  0x1.5555555555555p-5, // degree 4, Q_2[5]
  0x1.1111111111111p-7, // degree 5, Q_2[6]
  0x1.6c16c3fbb4213p-10, // degree 6, Q_2[7]
  0x1.a01a023ede0d7p-13, // degree 7, Q_2[8]
};

// Approximation for the accurate path of exp(z) for z=zh+zl,
// with |z| < 0.000130273 < 2^-12.88 and |zl| < 2^-42.6
// (assuming x^y does not overflow or underflow)
static inline void q_2 (double *hi, double *lo, double zh, double zl) {
  /* Let q[0]..q[7] be the coefficients of degree 0..7 of Q_2.
     The ulp of q[7]*z^7 is at most 2^-155, thus we can compute q[7]*z^7
     in double precision only.
     The ulp of q[6]*z^6 is at most 2^-139, thus we can compute q[6]*z^6
     in double precision only.
     The ulp of q[5]*z^5 is at most 2^-124, thus we can compute q[5]*z^5
     in double precision only. */
  double z = zh + zl;
  double q = __builtin_fma (Q_2[8], zh, Q_2[7]);

  q = __builtin_fma (q, z, Q_2[6]);

  q = __builtin_fma (q, z, Q_2[5]);

  // multiply q by z and add Q_2[3] + Q_2[4]
  a_mul (hi, lo, q, z);
  double t;
  fast_two_sum (hi, &t, Q_2[3], *hi);
  *lo += t + Q_2[4];

  // multiply hi+lo by zh+zl and add Q_2[2]
  d_mul (hi, lo, *hi, *lo, zh, zl);
  fast_two_sum (hi, &t, Q_2[2], *hi);
  *lo += t;

  // multiply hi+lo by zh+zl and add Q_2[1]
  d_mul (hi, lo, *hi, *lo, zh, zl);
  fast_two_sum (hi, &t, Q_2[1], *hi);
  *lo += t;

  // multiply hi+lo by zh+zl and add Q_2[0]
  d_mul (hi, lo, *hi, *lo, zh, zl);
  fast_two_sum (hi, &t, Q_2[0], *hi);
  *lo += t;
}

// same as exp_1, but with q_1 replaced by q_2 (more accurate)
static inline void exp_2 (double *hi, double *lo, double x) {

#define INVLOG2 0x1.71547652b82fep+12 /* |INVLOG2-2^12/log(2)| < 2^-43.4 */
  double k = __builtin_roundeven (x * INVLOG2);
  // since x <= 0x1.62e42fefa39fp+9 we have k <= 4194305

  double yh, yl;
  /* LOG2H and LOG2L are chosen such that k*LOGH and k*LOGL are exact,
     with |LOG2Hacc+LOG2Macc+LOG2Lacc - log(2)/2^12)| < 2^-132.546 */
#define LOG2Hacc 0x1.62e42ffp-13
#define LOG2Macc -0x1.718432ap-47
#define LOG2Lacc -0x1.b0e2633fe0685p-79
  yh = __builtin_fma (-k, LOG2Hacc, x); // exact
  yl = -k * LOG2Macc;                   // exact
  if (__builtin_fabs (yh) >= __builtin_fabs (yl))
    fast_two_sum (&yh, &yl, yh, yl);
  else
    fast_two_sum (&yh, &yl, yl, yh);
  yl -= k * LOG2Lacc;

  int64_t K = k;
  int64_t M = (K >> 12) + 0x3ff;
  int64_t i2 = (K >> 6) & 0x3f;
  int64_t i1 = K & 0x3f;

  double t1h = T1[i2][0], t1l = T1[i2][1], t2h = T2[i1][0], t2l = T2[i1][1];
  d_mul (hi, lo, t2h, t2l, t1h, t1l);

  double qh, ql;
  q_2 (&qh, &ql, yh, yl);

  d_mul (hi, lo, *hi, *lo, qh, ql);
  f64_u _d;

  _d.u = M << 52;
  *hi *= _d.f;
  *lo *= _d.f;
}

typedef union {double f; uint64_t u;} b64u64_u;

/* The following is a degree-11 polynomial generated by Sollya
   (file expm1_fast.sollya),
   which approximates expm1(x) with relative error bounded by 2^-67.183
   for |x| <= 0.125. */
static const double P[] = {
  0x0,                   // degree 0 (unused)
  0x1p0,                 // degree 1
  0x1p-1,                // degree 2
  0x1.5555555555555p-3,  // degree 3
  0x1.5555555555553p-5,  // degree 4
  0x1.1111111111bbcp-7,  // degree 5
  0x1.6c16c16c1f8a2p-10, // degree 6
  0x1.a01a0183a908bp-13, // degree 7
  0x1.a01a00383b80dp-16, // degree 8
  0x1.71e02a5f3b87p-19,  // degree 9
  0x1.27fcd07571d4ep-22, // degree 10
  0x1.969ce6c7ee119p-26, // degree 11
};

/* |x| <= 0.125, put in h + l a double-double approximation of expm1(x),
   and return the maximal corresponding absolute error.
   We also have |x| > 0x1.6a09e667f3bccp-53.
   With xmin=RR("0x1.6a09e667f3bccp-53",16), the routine
   expm1_fast_tiny_all(xmin,0.125,2^-64.29) in expm1.sage returns
   4.43378275371923e-20 < 2^-64.29, and
   expm1_fast_tiny_all(-0.125,-xmin,2^-64.13) returns
   4.94311016020191e-20 < 2^-64.13, which proves the relative
   error is bounded by 2^-64.13. */
static double
expm1_fast_tiny (double *h, double *l, double x)
{
  /* The maximal value of |P[4]*x^4/expm1(x)| over [-0.125,0.125]
     is less than 2^-13.495, thus we can compute the coefficients of degree
     4 or higher using double precision only. */
  double x2 = x * x, x4 = x2 * x2;
  double c10 = __builtin_fma (P[11], x, P[10]);
  double c8 = __builtin_fma (P[9], x, P[8]);
  double c6 = __builtin_fma (P[7], x, P[6]);
  double c4 = __builtin_fma (P[5], x, P[4]);
  c8 = __builtin_fma (c10, x2, c8);
  c4 = __builtin_fma (c6, x2, c4);
  c4 = __builtin_fma (c8, x4, c4);
  double t;
  // multiply c4 by x and add P[3]
  a_mul (h, l, c4, x);
  fast_two_sum (h, &t, P[3], *h);
  *l += t;
  // multiply (h,l) by x and add P[2]
  s_mul (h, l, x, *h, *l);
  fast_two_sum (h, &t, P[2], *h);
  *l += t;
  // multiply (h,l) by x and add P[1]
  s_mul (h, l, x, *h, *l);
  fast_two_sum (h, &t, P[1], *h);
  *l += t;
  // multiply (h,l) by x
  s_mul (h, l, x, *h, *l);
  return 0x1.d4p-65 * *h; // 2^-64.13 < 0x1.d4p-65
}

/* Given -0x1.2b708872320e2p+5 < x < -0x1.6a09e667f3bccp-53 or
   0x1.6a09e667f3bccp-53 < x < 0x1.62e42fefa39fp+9, put in h + l a
   double-double approximation of expm1(x), and return the maximal
   corresponding absolute error.
   The input tiny is true iff |x| <= 0.125. */
static double
expm1_fast (double *h, double *l, double x, int tiny)
{
  if (tiny) // |x| <= 0.125
    return expm1_fast_tiny (h, l, x);

  /* now -0x1.2b708872320e2p+5 < x < 0.125 or
     0.125  < x < 0x1.62e42fefa39fp+9: we approximate exp(x)
     and subtract 1 */

  exp_1 (h, l, x);
  /* h + l = exp(x) * (1 + eps) with |eps| < 2^-74.139 */
  double err1 = 0x1.d1p-75 * *h; // 2^-74.139 < 0x1.d1p-75
  double u;
  if (x >= 0) // implies h >= 1 and the fast_two_sum pre-condition holds
    fast_two_sum (h, &u, *h, -1.0);
  else
    fast_two_sum (h, &u, -1.0, *h); // x < 0 thus h <= 1
  *l += u;
  /* the error in the above fast_two_sum is bounded by 2^-105*|h|,
     with the new value of h */
  return err1 + 0x1p-105 * *h; // h >= 0
}

/* The following is a degree-16 polynomial generated by Sollya
   (file expm1_accurate.sollya),
   which approximates expm1(x) with relative error bounded by 2^-109.536
   for |x| <= 0.125. */
static const double Q[] = {
  0x1p0,                                         // degree 1: Q[0]
  0x1p-1,                                        // degree 2: Q[1]
  0x1.5555555555555p-3, 0x1.55555555554abp-57,   // degree 3: Q[2], Q[3]
  0x1.5555555555555p-5, 0x1.5555555529b52p-59,   // degree 4: Q[4], Q[5]
  0x1.1111111111111p-7, 0x1.111110fd7800cp-63,   // degree 5: Q[6], Q[7]
  0x1.6c16c16c16c17p-10, -0x1.f49f228e81422p-65, // degree 6: Q[8], Q[9]
  0x1.a01a01a01a01ap-13, 0x1.a1a3748b2ap-73,     // degree 7: Q[10], Q[11]
  0x1.a01a01a01a01ap-16,                         // degree 8: Q[12]
  0x1.71de3a556c733p-19,                         // degree 9: Q[13]
  0x1.27e4fb7789f9fp-22,                         // degree 10: Q[14]
  0x1.ae64567f5755ep-26,                         // degree 11: Q[15]
  0x1.1eed8efedba9bp-29,                         // degree 12: Q[16]
  0x1.612460b437492p-33,                         // degree 13: Q[17]
  0x1.93976857d992ap-37,                         // degree 14: Q[18]
  0x1.ae966f43fe1c7p-41,                         // degree 15: Q[19]
  0x1.ac8bc1457bf6dp-45,                         // degree 16: Q[20]
};

/* Accurate path for 0x1.6a09e667f3bccp-53 < |x| <= 0.125. */
static double
expm1_accurate_tiny (double x)
{
  /* exceptions below have between 46 and 57 identical bits after the
     round bit */
#define EXCEPTIONS 47
  static const double exceptions[EXCEPTIONS][3] = {
   {0x1.0b5d6cc46b3f8p-28, 0x1.0b5d6ccd251f9p-28, 0x1.fffffffffffffp-82},
   {-0x1.0f9b1c5ad2f3p-22, -0x1.0f9b1a1a7f6e3p-22, 0x1.fffffffffffffp-76},
   {-0x1.19e53fcd490dp-23, -0x1.19e53e96dffa9p-23, 0x1.fffffffffffffp-77},
   {0x1.1a4d6f93a29efp-24, 0x1.1a4d702f49f7dp-24, -0x1.fffffffffffffp-78},
   {-0x1.1a9dc8f6df10ap-47, -0x1.1a9dc8f6df0f7p-47, 0x1.fffffffffffffp-101},
   {0x1.2cf34db4807cdp-14, 0x1.2cf6114f2054bp-14, 0x1.fffffffffffffp-68},
   {-0x1.3988e1409212fp-51, -0x1.3988e1409212dp-51, -0x1.fffffffffffffp-105},
   {0x1.44c3d7c85bcf1p-14, 0x1.44c70fce6daabp-14, 0x1.fffffffffffffp-68},
   {-0x1.47b50a2a84ea8p-43, -0x1.47b50a2a84d05p-43, 0x1.fffffffffffffp-97},
   {0x1.51fce10251a48p-16, 0x1.51fdc02094ef7p-16, 0x1.fffffffffffffp-70},
   {-0x1.64808871369c2p-30, -0x1.6480886d55b0bp-30, 0x1.fffffffffffffp-84},
   {-0x1.6e9b2675a667ep-44, -0x1.6e9b2675a6577p-44, -0x1.fffffffffffffp-98},
   {-0x1.8154be277353ep-46, -0x1.8154be27734f5p-46, -0x1.fffffffffffffp-100},
   {0x1.8387d84827defp-38, 0x1.8387d8482c743p-38, 0x1.fffffffffffffp-92},
   {-0x1.8a8597b7c4b28p-23, -0x1.8a859557c5383p-23, -0x1.fffffffffffffp-77},
   {-0x1.92a19fd3ece36p-42, -0x1.92a19fd3ec943p-42, -0x1.fffffffffffffp-96},
   {-0x1.964a682912f4p-31, -0x1.964a68268e23fp-31, 0x1.fffffffffffffp-85},
   {-0x1.99ccc999fff07p-48, -0x1.99ccc999ffef3p-48, 0x1.fffffffffffffp-102},
   {0x1.a31972381bd0cp-23, 0x1.a31974e638221p-23, 0x1.fffffffffffffp-77},
   {-0x1.a8f783d749a8fp-4, -0x1.93aa1590d1e64p-4, -0x1.924af54b72c83p-108},
   {-0x1.ab86cb1743b75p-4, -0x1.95f8998ae5a65p-4, -0x1.84ea6e52b401ep-113},
   {-0x1.abb3b16c80ac4p-32, -0x1.abb3b16b1b63dp-32, 0x1.fffffffffffffp-86},
   {-0x1.b31e4dcde1e8ap-40, -0x1.b31e4dcde076dp-40, -0x1.fffffffffffffp-94},
   {-0x1.b935b38a6abadp-18, -0x1.b935547d3666dp-18, 0x1.fffffffffffffp-72},
   {0x1.bddfe561dbef3p-27, 0x1.bddfe5926531bp-27, -0x1.fffffffffffffp-81},
   {0x1.be2caeebfc83bp-4, 0x1.d761d8637563p-4, 0x1.a3cd02c39fb3ep-106},
   {-0x1.be9eacd95738dp-4, -0x1.a721c6d62e063p-4, 0x1.361e65cd9241p-107},
   {-0x1.c3263f6db7b48p-4, -0x1.ab30fc87097fap-4, 0x1.c47590934b57ap-106},
   {0x1.c58a7e3c93897p-4, 0x1.df9a92a4eb774p-4, -0x1.92f5f627f559fp-108},
   {-0x1.ca36132b4416p-4, -0x1.b182df1ecadb4p-4, 0x1.d689576ff9da1p-107},
   {-0x1.d097524a42e42p-4, -0x1.b7361f5082622p-4, -0x1.233581a73fd4ap-105},
   {-0x1.d4bb2250fc188p-19, -0x1.d4baecad344bfp-19, 0x1.fffffffffffffp-73},
   {-0x1.daf693d64fadap-4, -0x1.c075a87afb8a8p-4, -0x1.2aa7e4ef70195p-109},
   {-0x1.ddf3947c72332p-4, -0x1.c31ea77b4d57dp-4, -0x1.802e09c28d484p-106},
   {-0x1.dfeb80fca1157p-4, -0x1.c4def84730a0ep-4, 0x1.f5a1ebf9018f4p-107},
   {0x1.e0d50de7cdcecp-4, 0x1.fe31412377851p-4, 0x1.f881b1e44c357p-106},
   {-0x1.e6a0cc21f2c9fp-4, -0x1.cad5246110345p-4, -0x1.ffffffffffffep-58},
   {-0x1.e6b201f0d01f4p-4, -0x1.cae46c9e30824p-4, -0x1.182bebf9c627dp-104},
   {0x1.e923c188ea79bp-4, 0x1.03c5a420857cfp-3, 0x1.e63455fa8abf5p-113},
   {-0x1.e997e57006edcp-4, -0x1.cd76f688575e4p-4, -0x1.64015ad7add95p-107},
   {-0x1.ea5a8f57b2fc6p-4, -0x1.ce23adee8eaaap-4, -0x1.e4ae27bc3f8adp-105},
   {0x1.f359f8f048583p-13, 0x1.f369315ef3e8bp-13, -0x1.fffffffffffffp-67},
   {-0x1.f9c22c39aa1f4p-4, -0x1.dbc7c68016605p-4, 0x1.ffffffffffffdp-58},
   {-0x1.fab2bc8ad912p-4, -0x1.dc9c5f1ae8c7ap-4, -0x1.63b24a65b7a68p-105},
   {-0x1.ff504f1b8677cp-4, -0x1.e0afde3e0da82p-4, -0x1.197f08e0f1202p-107},
   {0x1p-52, 0x1.0000000000001p-52, -0x1.fffffffffffffp-106},
   {-0x1.bb67ae8584cabp-52, -0x1.bb67ae8584ca9p-52, -0x1.ffffffffffffep-106},
  };
  for (int i = 0; i < EXCEPTIONS; i++)
    if (x == exceptions[i][0])
      return exceptions[i][1] + exceptions[i][2];
#undef EXCEPTIONS

  double h, l, t;
  double x2 = x * x, x4 = x2 * x2;
  double c15 = __builtin_fma (Q[20], x, Q[19]);
  double c13 = __builtin_fma (Q[18], x, Q[17]);
  double c11 = __builtin_fma (Q[16], x, Q[15]);
  c13 = __builtin_fma (c15, x2, c13);
  // add Q[14]*x+c11*x2+c13*x4 to Q[13]
  fast_two_sum (&h, &l, Q[13], Q[14] * x + c11 * x2 + c13 * x4);
  // multiply h+l by x and add Q[12]
  s_mul (&h, &l, x, h, l);
  fast_two_sum (&h, &t, Q[12], h);
  l += t;
  // multiply h+l by x and add Q[10]+Q[11]
  s_mul (&h, &l, x, h, l);
  fast_two_sum (&h, &t, Q[10], h);
  l += t + Q[11];
  // multiply h+l by x and add Q[8]+Q[9]
  s_mul (&h, &l, x, h, l);
  fast_two_sum (&h, &t, Q[8], h);
  l += t + Q[9];
  // multiply h+l by x and add Q[6]+Q[7]
  s_mul (&h, &l, x, h, l);
  fast_two_sum (&h, &t, Q[6], h);
  l += t + Q[7];
  // multiply h+l by x and add Q[4]+Q[5]
  s_mul (&h, &l, x, h, l);
  fast_two_sum (&h, &t, Q[4], h);
  l += t + Q[5];
  // multiply h+l by x and add Q[2]+Q[3]
  s_mul (&h, &l, x, h, l);
  fast_two_sum (&h, &t, Q[2], h);
  l += t + Q[3];
  // multiply h+l by x and add Q[1]
  s_mul (&h, &l, x, h, l);
  fast_two_sum (&h, &t, Q[1], h);
  l += t;
  // multiply h+l by x
  s_mul (&h, &l, x, h, l);
  // multiply h+l by x
  s_mul (&h, &l, x, h, l);
  // add Q[0]*x = x
  fast_two_sum (&h, &t, x, h);
  l += t;
  return h + l;
}

static double expm1_accurate (double x)
{
  b64u64_u t = {.f = x};
  uint64_t ux = t.u, ax = ux & 0x7ffffffffffffffflu;

  if (ax <= 0x3fc0000000000000lu) // |x| <= 0.125
    return expm1_accurate_tiny (x);

  /* exceptions below have between 47 and 57 identical bits after the
     round bit */
#define EXCEPTIONS 96
  static const double exceptions[EXCEPTIONS][3] = {
   {-0x1.add1dce7cd5bcp-2, -0x1.5f0357a4cf6c6p-2, 0x1.398091600cd41p-105},
   {0x1.aca7ae8da5a7bp+0, 0x1.157d4acd7e557p+2, -0x1.fffffffffffffp-52},
   {0x1.d6336a88077aap+0, 0x1.51a8dff540ff7p+2, 0x1.78f1982b593afp-105},
   {-0x1.1397add4538acp-1, -0x1.aa3b02b31d93ep-2, 0x1.0f8991f2f9ffp-104},
   {0x1.00091a4a0dae5p+2, 0x1.ad0726fd1ccb3p+5, -0x1.ffffffffffffep-49},
   {0x1.273c188aa7b14p+2, 0x1.8f295a96ec6ebp+6, -0x1.fffffffffffffp-48},
   {0x1.0d73e6af47f36p+2, 0x1.097ad3a32b788p+6, -0x1.0f723a11eccf7p-98},
   {0x1.83d4bcdebb3f4p+2, 0x1.ab50b409c8aeep+8, 0x1.16719fcede453p-103},
   {0x1.a9da5c0e731eap+2, 0x1.836c253ffa44ep+9, -0x1.c092bc48c07c3p-96},
   {-0x1.02b72fbea16ep-2, -0x1.c93d7beeed889p-3, -0x1.5fedfe9ad705ap-107},
   {0x1.016e82ceda359p+1, 0x1.9e37fb31fd5fcp+2, 0x1.3709b5649e622p-102},
   {0x1.76e7e5d7b6eacp+3, 0x1.de7bd6751029ap+16, 0x1.d1765ed0dbee1p-89},
   {0x1.62f71c4656b61p-1, 0x1.0012ecb039c9cp+0, 0x1.01dc6b104a893p-105},
   {-0x1.ea16274b0109bp-3, -0x1.b3dbaf5230568p-3, -0x1.21f261d234775p-106},
   {-0x1.343d5853ab1bap-3, -0x1.1e2a26c6cbcffp-3, -0x1.43c9bfa598339p-109},
   {0x1.08f51434652c3p+4, 0x1.daac439b157e5p+23, 0x1.c6823badae774p-84},
   {0x1.1d5c2daebe367p+4, 0x1.a8c02e174c315p+25, -0x1.de0fc9395bbd4p-83},
   {0x1.634b2dd7eb0a3p+4, 0x1.0684c2e7b00fcp+32, 0x1.1e355e7edc3c3p-72},
   {-0x1.789d025948efap-2, -0x1.3b1ee1f952dcdp-2, 0x1.ffffffffffffdp-56},
   {-0x1.dc2b5df1f7d3dp-1, -0x1.35fe01788d71cp-1, 0x1.dc83c7a84cf5fp-108},
   {0x1.2ee70220fb1c5p+5, 0x1.8aa92bc84ff91p+54, 0x1.3c264141f8e99p-54},
   {0x1.1a0408712e00ap-2, 0x1.44acc499153ccp-2, 0x1.758d621c3b9dep-106},
   {0x1.8172a0e02f90ep-2, 0x1.d404e97601d65p-2, -0x1.ffffffffffffap-56},
   {-0x1.22e24fa3d5cf9p-1, -0x1.bbd1d708f42adp-2, 0x1.dd168cf650e6cp-107},
   {0x1.fde31a71ddba9p-3, 0x1.217c79b0566b5p-2, 0x1.066766014f376p-106},
   {0x1.067b7708b71b9p-3, 0x1.180b891078de3p-3, 0x1.0702eef85fd7fp-105},
   {-0x1.0ce0a43467d9fp-3, -0x1.f7f3398ad73bdp-4, 0x1.82fd5865e6c2p-106},
   {0x1.57ecb12f42f09p-3, 0x1.767d3ffbe9a82p-3, 0x1.699b9b978bb19p-104},
   {0x1.a20bbf6967c11p-3, 0x1.cfc51d8639b8fp-3, 0x1.0a73acd1aa664p-105},
   {-0x1.cddf723d3e52fp-3, -0x1.9d7ec7df33dbcp-3, -0x1.4c519851f4cf7p-106},
   {0x1.28d1885215445p-3, 0x1.3f67cb950a619p-3, 0x1.fffffffffffeep-57},
   {-0x1.bc21f0ba4ae83p-3, -0x1.8f4678512ce5ep-3, 0x1.e4d526363d49p-107},
   {-0x1.58db2327d4e2bp-3, -0x1.3d627e6f47f5ap-3, 0x1.629cd619c9258p-104},
   {-0x1.f31bfe026a32ep-2, -0x1.8b0b6b63cdd01p-2, 0x1.37751462f58edp-106},
   {0x1.6587e74ac8c65p+0, 0x1.854e6c6c05b13p+1, -0x1.2c5c28c210ab2p-103},
   {0x1.005ae04256babp-1, 0x1.4cbb1357e7a3dp-1, 0x1.10f83e22a66fcp-106},
   {0x1.a3a7add74f25ap-2, 0x1.0359f11a22a8dp-1, 0x1.1318aadf4a74fp-104},
   {0x1.04ac36c54a838p-3, 0x1.15fd2bf20273bp-3, 0x1.1bfee381a9d05p-104},
   {0x1.8eaa8cb0d0f38p-3, 0x1.b81d26b109c8p-3, 0x1.d8bfad7e2b897p-104},
   {0x1.bcab27d05abdep-2, 0x1.166ce703b05e9p-1, 0x1.dfe7b252154edp-106},
   {0x1.22a9d3042f3bcp-3, 0x1.384d130e4e667p-3, 0x1.ffffffffffff7p-57},
   {-0x1.74c12f94c4363p-3, -0x1.54cccbecda5e9p-3, 0x1.7a2c5613c769ep-104},
   {0x1.27f4980d511ffp-2, 0x1.5728eea2bbdddp-2, -0x1.45704d346de3ap-108},
   {0x1.8bbe2fb45c151p-2, 0x1.e3186ba9d4d49p-2, 0x1.fffffffffffffp-56},
   {-0x1.43f9e7a1919fep-3, -0x1.2ba68646726cfp-3, -0x1.ffffffffffff7p-57},
   {-0x1.938f3a33191ffp-2, -0x1.4d87fe71b1badp-2, 0x1.3df9b0557c98ap-104},
   {-0x1.474d4de7c14bbp-2, -0x1.182619ac31282p-2, 0x1.31d9b27789f7p-104},
   {0x1.81f63829c4e0ap-3, 0x1.a8ba20a181bb4p-3, 0x1.dcb7557eb5134p-105},
   {0x1.c195b6198a1d2p-3, 0x1.f6c011cfb84c9p-3, 0x1.d430b1fc276fdp-105},
   {0x1.cd3848dec31b2p-3, 0x1.02a4321c73106p-2, 0x1.7bef5eff06922p-105},
   {0x1.d707029bb59d9p-2, 0x1.2b092fef66a7ep-1, 0x1.e27333f64e0b9p-106},
   {0x1.f6e4c3ced7c72p-3, 0x1.1d0232e560f38p-2, 0x1.644b7f5399dfp-107},
   {0x1.fab8ff1fa05f1p-1, 0x1.b0be8b6d80446p+0, 0x1.68f6a55a21bcap-103},
   {0x1.accfbe46b4efp-1, 0x1.4f85c9783dce1p+0, -0x1.4f3d25cc4570bp-107},
   {-0x1.b8144d498cc5bp-3, -0x1.8c024d0aa27b3p-3, 0x1.0702d4719540dp-108},
   {0x1.4e88c5accfda5p-3, 0x1.6b68447b2f2fdp-3, -0x1.2cb8fff116072p-110},
   {0x1.7d7fc2e4f5fccp-3, 0x1.a3583db6ebf94p-3, 0x1.1d43ccdb25e74p-106},
   {0x1.1c38132777b26p-2, 0x1.4794729ba52adp-2, -0x1.5abd16c70b908p-106},
   {-0x1.119aae6072d39p-2, -0x1.e033b6b48a1fbp-3, -0x1.5924cbdac85bbp-107},
   {0x1.0727af5fee8f6p-1, 0x1.5806551a5d846p-1, 0x1.05a62956baf36p-104},
   {0x1.e9375280398bbp+0, 0x1.70a159662f223p+2, 0x1.b23b8a707f1cdp-103},
   {-0x1.7725e67d94824p-3, -0x1.56cb40d337d1dp-3, 0x1.382199a20e409p-104},
   {-0x1.3b89bb1b787cdp-3, -0x1.246e9c005c036p-3, -0x1.4c6ed3df37a13p-107},
   {-0x1.8aeb636f3ce35p-3, -0x1.672d45e082548p-3, 0x1.82b6b66e03876p-110},
   {-0x1.d3f3799439415p-3, -0x1.a256cd99e8292p-3, -0x1.07e9d5d1fbf62p-108},
   {0x1.704f3cd72bc67p-3, 0x1.93816aeb7ae9ep-3, 0x1.4036d6435d65p-104},
   {-0x1.bf6548c02eec9p-3, -0x1.91e66bbf6409fp-3, -0x1.0f5e23310d5e8p-104},
   {0x1.8a31fa78c68afp-3, 0x1.b2b0210d13ed3p-3, 0x1.a132d283b2aa6p-107},
   {-0x1.290ea09e36479p-3, -0x1.1484b3cd038fp-3, -0x1.09bec3f4113eep-111},
   {0x1.a065fefae814fp-3, 0x1.cdc010de032e4p-3, 0x1.ff41a6986a845p-106},
   {0x1.4032183482ed7p-3, 0x1.5a955756bd3e1p-3, 0x1.848861561a9bp-106},
   {0x1.223469ea438e4p-3, 0x1.37c5c66727519p-3, 0x1.08e9c703aa567p-104},
   {-0x1.f193dbe5f18bbp-3, -0x1.b9be811308062p-3, 0x1.8e3bbdc8e1fd7p-106},
   {-0x1.82b5dfaf59b4cp-2, -0x1.4213802eb28ffp-2, 0x1.ffffffffffffdp-56},
   {0x1.d086543694c5ap-1, 0x1.7a417a07cafafp+0, 0x1.62cf4b32655bdp-106},
   {0x1.37c9e03e6099dp+0, 0x1.30a9340959f93p+1, 0x1.bd43cb2443c73p-103},
   {-0x1.2a9cad9998262p+0, -0x1.60870b262bf1cp-1, -0x1.e757fe830d60ep-109},
   {-0x1.ac71ace4f979p-3, -0x1.8298cb774359bp-3, 0x1.b5ca0b5a8abeep-105},
   {0x1.677e81300d1f4p-3, 0x1.88f9dbf1b4e57p-3, 0x1.8a69b407735c4p-105},
   {0x1.4297ec53f6b7fp-1, 0x1.c16640ad39959p-1, 0x1.ffffffffffffdp-55},
   {-0x1.0a54d87783d6fp+0, -0x1.4b1887d4d477cp-1, 0x1.d81f352752164p-108},
   {0x1.1f0da93354198p+7, 0x1.0bd73b73fc74cp+207, 0x1.588526e93304cp+103},
   {0x1.7a60ee15e3e9dp+6, 0x1.62e4dc3bbf53fp+136, 0x1.ae7c8ed9b6bcbp+30},
   {0x1.0bc04af1b09f5p+9, 0x1.7b1d97c902985p+772, 0x1.551dfecc05bd4p+667},
   {0x1.9e7b643238a14p+8, 0x1.f5da7fe652978p+597, 0x1.0429700e71228p+494},
   {0x1.3f37fb551e418p+9, 0x1.0ccdc52c38712p+921, 0x1.376e0ad4f23b9p+818},
   {0x1.2da9e5e6af0bp+8, 0x1.27d6fe867d6f6p+435, 0x1.0a1d500c39996p+330},
   {0x1.556c678d5e976p+7, 0x1.37e7ac4e7f9b3p+246, 0x1.01a99afd82b06p+142},
   {0x1.519fd95037e31p+6, 0x1.b53c2f00bb322p+121, 0x1.fa70a379db2fbp+18},
   {0x1.54cd1fea7663ap+7, 0x1.c90810d354618p+245, 0x1.2925a9627fb2cp+136},
   {0x1.6474c604cc0d7p+6, 0x1.7a8f65ad009bdp+128, -0x1.0b611958ec877p+21},
   {0x1.d6479eba7c971p+8, 0x1.62a88613629b6p+678, -0x1.3f69a2085428cp+569},
   {0x1.7945e34b18a9ap+7, 0x1.1b0e4936a8c9bp+272, -0x1.f983e7e9b16f1p+167},
   {0x1.c44ce0d716a1ap+4, 0x1.b890ca8636ae2p+40, -0x1.bde9e7639f651p-68},
   {0x1.c7206c1b753e4p+8, 0x1.8670de0b68cadp+656, -0x1.7599cebd802f7p+549},
   {0x1.89d56a0c38e6fp+5, 0x1.0410c95b580b9p+71, -0x1.7d1a66d4c94f2p-40},
  };
  for (int i = 0; i < EXCEPTIONS; i++)
    if (x == exceptions[i][0])
      return exceptions[i][1] + exceptions[i][2];
#undef EXCEPTIONS

  /* now -0x1.2b708872320e2p+5 < x < 0.125 or
     0.125  < x < 0x1.62e42fefa39fp+9: we approximate exp(x)
     and subtract 1 */

  double h, l;
  exp_2 (&h, &l, x);
  double u;
  if (x >= 0) // implies h >= 1 and the fast_two_sum pre-condition holds
    fast_two_sum (&h, &u, h, -1.0);
  else
    fast_two_sum (&h, &u, -1.0, h); // x < 0 thus h <= 1
  l += u;
  /* the error in the above fast_two_sum is bounded by 2^-105*|h|,
     with the new value of h */
  return h + l;
}

double
cr_expm1 (double x)
{
  b64u64_u t = {.f = x};
  uint64_t ux = t.u, ax = ux & 0x7ffffffffffffffflu;

  if (__builtin_expect (ux >= 0xc042b708872320e2, 0))
  {
    // x = -NaN or x <= -0x1.2b708872320e2p+5
    if ((ux >> 52) == 0xfff) // -NaN or -Inf
      return (ux > 0xfff0000000000000lu) ? x : -1.0;
    // for x <= -0x1.2b708872320e2p+5, expm1(x) rounds to -1 to nearest
    return -1.0 + 0x1p-54;
  }
  else if (__builtin_expect (ax >= 0x40862e42fefa39f0, 0))
  {
    // x = +NaN or x >= 0x1.62e42fefa39fp+9
    if ((ux >> 52) == 0x7ff) // +NaN
      return x;
    // for x >= 0x1.62e42fefa39fp+9, expm1(x) rounds to +Inf to nearest
    return 0x1.fffffffffffffp+1023 * x;
  }
  else if (ax <= 0x3ca6a09e667f3bcc) // |x| <= 0x1.6a09e667f3bccp-53
    /* then expm1(x) rounds to x (to nearest), with Taylor expansion
       x + x^2/2 + ... */
  {
    if (ax <= 0x3c96a09e667f3bcc)
      /* |x| <= 0x1.6a09e667f3bccp-54: x^2 < 1/2 ulp(x), we have to deal with
         -0 apart since fma (-0, -0, -0) is (+0) + (-0) which evaluates to +0
         for some rounding modes */
      return (x == 0) ? x : __builtin_fma (x, x, x);
    else
      /* 0x1p-53 <= |x| <= 0x1.6a09e667f3bccp-53: x/4 is exactly
         representable, and x^2/4 < 1/2 ulp(x) */
      return __builtin_fma (x, x * 0.25, x);
  }

  /* -0x1.2b708872320e2p+5 < x < -0x1.6a09e667f3bccp-53 or
     0x1.6a09e667f3bccp-53 < x < 0x1.62e42fefa39fp+9 */

  double err, h, l;
  err = expm1_fast (&h, &l, x, ax <= 0x3fc0000000000000lu);
  double left = h + (l - err);
  double right = h + (l + err);
  if (left == right)
    return left;

  return expm1_accurate (x);
}
