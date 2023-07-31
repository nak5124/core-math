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

// FIXME: remove skipped in check_worst_uni.c

// #define TRACE 0x1.5178df87adcfap-4
#define TRACE 0x1.0b5d6cc46b3f8p-28

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
  /* The maximal value of |P_fast_tiny[4]*x^4/expm1(x)| over [-0.125,0.125]
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
   corresponding absolute error. */
static double
expm1_fast (double *h, double *l, double x)
{
  b64u64_u t = {.f = x};
  // FIXME: pass ax from cr_expm1()
  uint64_t ux = t.u, ax = ux & 0x7ffffffffffffffflu;

  if (ax <= 0x3fc0000000000000lu) // |x| <= 0.125
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
  double h, l;
  double x2 = x * x, x4 = x2 * x2;
  double c15 = __builtin_fma (Q[20], x, Q[19]);
  double c13 = __builtin_fma (Q[18], x, Q[17]);
  double c11 = __builtin_fma (Q[16], x, Q[15]);
  double c9 = __builtin_fma (Q[14], x, Q[13]);
  c13 = __builtin_fma (c15, x2, c13);
  c9 = __builtin_fma (c11, x2, c9);
  c9 = __builtin_fma (c13, x4, c9);
  double t;
  // multiply c9 by x and add Q[12]
  a_mul (&h, &l, c9, x);
  fast_two_sum (&h, &t, Q[12], h);
  l += t + Q[11];
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
  if (x == TRACE) printf ("h=%la l=%la\n", h, l);
  // multiply h+l by x
  s_mul (&h, &l, x, h, l);
  if (x == TRACE) printf ("h=%la l=%la\n", h, l);
  // add Q[0]*x = x
  fast_two_sum (&h, &t, x, h);
  l += t;
  if (x == TRACE) printf ("h=%la l=%la\n", h, l);
  return h + l;
}

static double expm1_accurate (double x)
{
  b64u64_u t = {.f = x};
  // FIXME: pass ax from cr_expm1()
  uint64_t ux = t.u, ax = ux & 0x7ffffffffffffffflu;

  if (ax <= 0x3fc0000000000000lu) // |x| <= 0.125
    return expm1_accurate_tiny (x);

  return -2;
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
    return 0x1.fffffffffffffp+1023 + 0x1.fffffffffffffp+1023;
  }
  else if (ax <= 0x3ca6a09e667f3bcc) // |x| <= 0x1.6a09e667f3bccp-53
    /* then expm1(x) rounds to x (to nearest), with Taylor expansion
       x + x^2/2 + ... */
  {
    if (ax < 0x3ca0000000000000)
    {
      /* |x| < 0x1p-53: x^2 < 1/2 ulp(x), we have to deal with -0 apart
         since fma (-0, -0, -0) is (+0) + (-0) which evaluates to +0
         for some rounding modes */
      return (x == 0) ? x : __builtin_fma (x, x, x);
    }
    else
      /* 0x1p-53 <= |x| <= 0x1.6a09e667f3bccp-53: x/4 is exactly
         representable, and x^2/4 < 1/2 ulp(x) */
      return __builtin_fma (x, x * 0.25, x);
  }

  /* -0x1.2b708872320e2p+5 < x < -0x1.6a09e667f3bccp-53 or
     0x1.6a09e667f3bccp-53 < x < 0x1.62e42fefa39fp+9 */

  double err, h, l;
  err = expm1_fast (&h, &l, x);
  double left = h + (l - err);
  double right = h + (l + err);
  if (x==TRACE) printf ("x=%la h=%la l=%la err=%la left=%la right=%la\n", x, h, l, err, left, right);
  if (left == right)
    return left;

  if (x==TRACE) printf ("fast path failed\n");
  
  return expm1_accurate (x);
}
