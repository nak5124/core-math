/* Correctly rounded base-10 exponential function for binary64 values.

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

/* References:
   [5] Towards a correctly-rounded and fast power function in binary64
       arithmetic, Tom Hubrecht, Claude-Pierre Jeannerod, Paul Zimmermann,
       ARITH 2023 - 30th IEEE Symposium on Computer Arithmetic, 2023.
       Detailed version (with full proofs) available at
       https://inria.hal.science/hal-04159652.
 */

#include <stdint.h>

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
// We can ignore al * bl when assuming al <= ulp(ah) and bl <= ulp(bh)
static inline void d_mul(double *hi, double *lo, double ah, double al,
                         double bh, double bl) {
  double s, t;

  a_mul(hi, &s, ah, bh);
  t = __builtin_fma(al, bh, s);
  *lo = __builtin_fma(ah, bl, t);
}

static inline void fast_two_sum(double *hi, double *lo, double a, double b) {
  double e;

  *hi = a + b;
  e = *hi - a;
  *lo = b - e;
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
   over [-0.000130273,0.000130273]
   with absolute error < 2^-74.346 (see sollya/Q_1.sollya). */
static const double Q_1[] = {0x1p0,                 /* degree 0 */
                             0x1p0,                 /* degree 1 */
                             0x1p-1,                /* degree 2 */
                             0x1.5555555995d37p-3,  /* degree 3 */
                             0x1.55555558489dcp-5   /* degree 4 */
};

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

/* Given (zh,zl) such that |zh+zl| < 0.000130273 and |zl| < 2^-42.7260,
   this routine puts in qh+ql an approximation of exp(zh+zl) such that

   | (qh+ql) / exp(zh+zl) - 1 | < 2^-74.169053

   See Lemma 6 from reference [5].
*/
static inline void q_1 (double *qh, double *ql, double zh, double zl) {
  double z = zh + zl;
  double q = __builtin_fma (Q_1[4], zh, Q_1[3]); /* q3+q4*z */

  q = __builtin_fma (q, z, Q_1[2]); /* q2+q3*z+q4*z^2 */

  fast_two_sum (qh, ql, Q_1[1], q * z);

  d_mul (qh, ql, zh, zl, *qh, *ql);

  fast_sum (qh, ql, Q_1[0], *qh, *ql);
}

#define NAN (0.0/0.0)

/* Code adapted from pow.c, by removing the argument s (sign is always
   positive here).
   When called from exp10_fast(), the condition
   |rl/rh| is fulfilled because |rl| <= ulp(rh).
   Also the condition |rl| < 2^-14.4187 is fulfilled because |rl| <= 2^-43.
   We also have -0x1.74385446d71c3p+9 <= rh <= 0x1.62e42fefa39efp+9.

   Given RHO1 <= rh <= RHO2, |rl/rh| < 2^-23.8899 and |rl| < 2^-14.4187,
   this routine computes an approximation eh+el of exp(rh+rl) such that:

   | (eh+el) / exp(rh+rl) - 1 | < 2^-74.16.

   Moreover |el/eh| <= 2^-41.7.

   See Lemma 7 from reference [5].

   The result eh+el is multiplied by s (which is +1 or -1).
*/
static inline void
exp_1 (double *eh, double *el, double rh, double rl) {

#define RHO0 -0x1.74910ee4e8a27p+9
#define RHO1 -0x1.577453f1799a6p+9
#define RHO2 0x1.62e42e709a95bp+9
#define RHO3 0x1.62e4316ea5df9p+9

  if (__builtin_expect(rh > RHO2, 0)) {
    /* since rh <= 0x1.62e42fefa39efp+9 when called from exp10_fast(),
       we can't have rh > RHO3 */
    *eh = *el = NAN; // delegate to the accurate step
    return;
  }

  if (__builtin_expect(rh < RHO1, 0)) {
    // since -0x1.74385446d71c3p+9 <= rh, we can't have rh < RHO0
    // RHO0 <= rh < RHO1: delegate to the accurate step
    *eh = *el = NAN;
    return;
  }

#define INVLOG2 0x1.71547652b82fep+12
  double k = __builtin_roundeven (rh * INVLOG2);

  double kh, kl;
#define LOG2H 0x1.62e42fefa39efp-13
#define LOG2L 0x1.abc9e3b39803fp-68
  s_mul (&kh, &kl, k, LOG2H, LOG2L);

  double yh, yl;
  fast_two_sum (&yh, &yl, rh - kh, rl);
  yl -= kl;

  int64_t K = k; /* Note: k is an integer, this is just a conversion. */
  int64_t M = (K >> 12) + 0x3ff;
  int64_t i2 = (K >> 6) & 0x3f;
  int64_t i1 = K & 0x3f;

  double t1h = T1[i2][0], t1l = T1[i2][1], t2h = T2[i1][0], t2l = T2[i1][1];
  d_mul (eh, el, t2h, t2l, t1h, t1l);

  double qh, ql;
  q_1 (&qh, &ql, yh, yl);

  d_mul (eh, el, *eh, *el, qh, ql);
  f64_u _d;

  /* we should have 1 < M < 2047 here, since we filtered out
     potential underflow/overflow cases at the beginning of this function */

  _d.u = M << 52;
  *eh *= _d.f;
  *el *= _d.f;
}

/****************** end of code copied from pow.[ch] *************************/

typedef union {double f; uint64_t u;} b64u64_u;

/* Put into h+l a double-double approximation of 10^x, and return a bound
   on the absolute error.
   Assumes -0x1.434e6420f4374p+8 < x < -0x1.bcb7b1526e50cp-55
   or 0x1.bcb7b1526e50cp-55 < x < 0x1.34413509f79ffp+8.

   FIXME: instead of multiplying x by log(10), giving rh+rl, and then in exp_1
   we  compute k = round(rh*2^12/log(2)), we could compute k directly as
   round(x*2^12*log(10)/log(2)) as in exp10_accurate(). But this would require
   to perform the error analysis again.
*/
static double
exp10_fast (double *h, double *l, double x)
{
  double rh, rl;
  /* first multiply x by an approximation of log(10):
     | LOG10H + LOG10L - log(10) | < 2^-106.3 */
#define LOG10H 0x1.26bb1bbb55516p+1
#define LOG10L -0x1.f48ad494ea3e9p-53
  s_mul (&rh, &rl, x, LOG10H, LOG10L);
  /* The rounding error from s_mul is bounded by ulp(rl).
     Since |x| < 0x1.434e6420f4374p+8, we have |rh| < 744.5
     thus |rl| <= ulp(rh) <= 2^-43, and ulp(rl) <= 2^-95.
     The approximation error from LOG10H+LOG10L is bounded by
     |x|*2^-106.3 < 2^-97.96.
     Moreover we have -0x1.74385446d71c3p+9 <= rh <= 0x1.62e42fefa39efp+9.
     Thus:
     | rh + rl - log(10)*x | < 2^-95 + 2^-97.96 < 2^-94.82 */
  exp_1 (h, l, rh, rl);
  /* We have from exp_1():

     | h + l - exp(rh + rl) | < 2^-74.16 * |h + l|

     and since |l/h| <= 2^-41.7:

     | h + l - exp(rh + rl) | < 2^-74.16 * |h| * (1 + 2^-41.7)
                              < 2^-74.159 * |h|                   (1)

     Since rh + rl = log(10)*x - eps with |eps| < 2^-94.82,
     10^x = exp(rh+rl) * exp(eps) thus:
     | h + l - 10^x | < 2^-74.159 * |h| + |exp(rh+rl)| * |1 - exp(eps)|
     and from (1) we get |exp(rh + rl)| < |h| + |l| + 2^-74.159 * |h|
                                        < (1 + 2^-41.7 + 2^-74.159) * |h|
                                        < 1.01 * |h|, we obtain:
     | h + l - 10^x | < 2^-74.159 * |h| + 1.01 * |1 - exp(eps)| * |h|
                      < 2^-74.158 * |h| */
  return 0x1.cbp-75 * *h; /* 2^-74.158 < 0x1.cbp-75 */
}

/* assumes -0x1.434e6420f4374p+8 < x < -0x1.bcb7b1526e50ep-56
   or 0x1.bcb7b1526e50ep-56 < x < 0x1.34413509f79ffp+8 */
static double
exp10_accurate (double x)
{
#define EXCEPTIONS 139
  static const double exceptions[EXCEPTIONS][3] = {
    {-0x1.21f6fc63d1e5p+8, 0x1.b0e12d8b3e4bp-964, 0x1.f55bdebc1959cp-1067},
    {-0x1.e095699d01895p+7, 0x1.b3e5f2b6612dfp-799, 0x1.fffffffffffffp-853},
    {-0x1.b5e974e4bcc37p+7, 0x1.900b0d01e2f71p-728, -0x1.ff977cb91c0cfp-833},
    {-0x1.77d933c1a88e1p+7, 0x1.a8639e89f5e46p-625, -0x1.3ee0ffcfa8744p-731},
    {-0x1.da5b10d8689fdp+6, 0x1.0a267404067fbp-394, -0x1.333ba46e440b3p-499},
    {-0x1.aac4ae5e3fa54p+5, 0x1.ba1f963d6e5c6p-178, -0x1.75ee4aff24726p-284},
    {-0x1.97e160d9b2526p+5, 0x1.8c8b39cbd7772p-170, -0x1.5ee60d91a567ep-276},
    {-0x1.1b8269857f23p+5, 0x1.35d07db329e0fp-118, -0x1.ffffffffffffep-172},
    {-0x1.a66820b6bbfbcp+4, 0x1.3b1b5689138a6p-88, 0x1.b27566dbd9f58p-195},
    {-0x1.a4de864908099p+4, 0x1.892503a8d2d08p-88, 0x1.8f4a8601f2b25p-192},
    {-0x1.8555fbe43c9fap+4, 0x1.1f3203f445e66p-81, 0x1.2081a50ae0576p-184},
    {-0x1.4cd4af2fca2b4p+4, 0x1.dce2dfaaf67c9p-70, 0x1.fffffffffffffp-124},
    {-0x1.0e841fe04c639p+4, 0x1.c8c3d2c1783ddp-57, -0x1.8b41b5d7ee523p-163},
    {-0x1.701d8d4d5f3efp+3, 0x1.b95e9b901a3ddp-39, 0x1.5222f11b19512p-144},
    {-0x1.18c0d46ea66c8p+3, 0x1.cf0578c8e8b8fp-30, -0x1.0737ba64e313bp-134},
    {-0x1.aa5575135e2d3p+2, 0x1.d43c2f5ee8437p-23, -0x1.c78c412807822p-128},
    {-0x1.7c988946c9573p+1, 0x1.16b2db95cf769p-10, 0x1.31dffbfdd2567p-115},
    {-0x1.3d75c9fa059e2p+1, 0x1.b1dcf92c18a36p-9, 0x1.a7524a9b64366p-114},
    {-0x1.e7cff3e41caadp+0, 0x1.9751483c97c78p-7, -0x1.ccc422187d9c5p-112},
    {-0x1.69374506458f3p+0, 0x1.3df91b68547c6p-5, 0x1.aef89c358e2bap-109},
    {-0x1.d18176754aac7p-1, 0x1.f8db9b16b1432p-4, 0x1.15e197166aa6fp-109},
    {-0x1.1416c72a588a6p-1, 0x1.27d838f22d0ap-2, -0x1.b2a409fb9c0bcp-121},
    {-0x1.c03419f51b93ep-2, 0x1.75c468428b7cp-2, 0x1.5be35b21b3043p-109},
    {-0x1.c360cdde773f7p-3, 0x1.343a134e572cep-1, 0x1.2411f7a782ef7p-105},
    {-0x1.a9cf11e5adbc5p-4, 0x1.9301d2902eba7p-1, -0x1.e29cc5e9bb647p-109},
    {-0x1.4c763938f02acp-4, 0x1.a8b846065f18p-1, -0x1.48aa7b519ca4cp-105},
    {-0x1.485b2ac5f7ec9p-5, 0x1.d2dc5e913a8fbp-1, 0x1.fffffffffffffp-55},
    {-0x1.5cccb7c4ac953p-6, 0x1.e781cf5117abbp-1, -0x1.4f2d6afdd582p-106},
    {-0x1.5b25114a07a72p-6, 0x1.e79ed6c92103cp-1, -0x1.0b8d27ba7b623p-107},
    {-0x1.5935097c1e6a6p-6, 0x1.e7c0d61c044e2p-1, 0x1.33aee3f43fb47p-105},
    {-0x1.3b95082297ea7p-6, 0x1.e9c9ca3499b0dp-1, -0x1.faf2495ae45b4p-108},
    {-0x1.39e42a1447eap-6, 0x1.e9e796e7e56ffp-1, 0x1.3468b204e0b39p-105},
    {-0x1.da040c968a4eep-7, 0x1.ef3a101770015p-1, 0x1.fffffffffffffp-55},
    {-0x1.b44e17164ce91p-7, 0x1.f08a757b7706ep-1, 0x1.b65e92030dcebp-106},
    {-0x1.26c5d9aa6e074p-7, 0x1.f580f7788ee3bp-1, -0x1.7e6fe164bf03fp-106},
    {-0x1.f6f96f005fd47p-8, 0x1.f708106038d25p-1, -0x1.fffffffffffffp-55},
    {-0x1.e435674dbad28p-8, 0x1.f75cff6484afp-1, 0x1.94d795665fbfep-107},
    {-0x1.be699eeec2c4ep-8, 0x1.f8083b79eaacfp-1, 0x1.8c838de5b4e09p-104},
    {-0x1.9ef920bed31b3p-8, 0x1.f896d78214048p-1, 0x1p-54},
    {-0x1.99ca2c1ae2284p-9, 0x1.fc53cfff2f3a2p-1, 0x1.249748240eac9p-104},
    {-0x1.4cadba1297213p-9, 0x1.fd04369659d16p-1, 0x1.36dbd33474872p-105},
    {-0x1.f6fd8cd6c9011p-10, 0x1.fdbe302f48bfep-1, -0x1.33ebd53aff143p-104},
    {-0x1.ebb11d32c9493p-10, 0x1.fdcb23b597a6fp-1, -0x1.df2e048ab83f2p-111},
    {-0x1.95334e650bbb9p-10, 0x1.fe2e5331a7495p-1, 0x1.db71cc20349a6p-105},
    {-0x1.6ee0b282f77ecp-10, 0x1.fe5a4bd73731bp-1, -0x1.ffffffffffffep-55},
    {-0x1.4c7a2be09b10ep-11, 0x1.ff40c0192f17p-1, 0x1.19f3dfca93d21p-107},
    {-0x1.e12494018e44cp-12, 0x1.ff7596cdd0f69p-1, -0x1.841150c917e2ep-106},
    {-0x1.a1b18d3a28957p-12, 0x1.ff87d54fafc19p-1, -0x1.f243e08383ddbp-108},
    {-0x1.63df14c04ab23p-12, 0x1.ff999cab1babbp-1, -0x1.02e5cd8b9642ep-107},
    {-0x1.126b256835fc5p-12, 0x1.ffb10a2f910ebp-1, 0x1.f7c9b471942d7p-105},
    {-0x1.7a7f33cc3fd0bp-13, 0x1.ffc98a8e780fbp-1, 0x1.ddaec6371e813p-106},
    {-0x1.21c4325902ac5p-13, 0x1.ffd64e4f403edp-1, -0x1.7cdc831614cf1p-111},
    {-0x1.69924ee4c9dbap-14, 0x1.ffe5fc447bee3p-1, 0x1.e9ee3902b414dp-105},
    {-0x1.3ad0c5feeab03p-14, 0x1.ffe9596287cc7p-1, 0x1.e537bc0e8b1bdp-104},
    {-0x1.17362e953393bp-14, 0x1.ffebe92022beap-1, 0x1.26bbec2bfbd34p-107},
    {-0x1.cc32651a2c138p-15, 0x1.ffef71b2b17b4p-1, 0x1.1340870d74334p-104},
    {-0x1.a0584cc4fdb48p-15, 0x1.fff1058a696f5p-1, 0x1.1896e76b1b0c4p-105},
    {-0x1.6506061aae6f7p-15, 0x1.fff327da56828p-1, -0x1.2ff1bf0e8acd8p-107},
    {-0x1.45ddb10382e3fp-15, 0x1.fff446cc3068dp-1, -0x1.437237d73a494p-109},
    {-0x1.fdebe98f7bd2bp-19, 0x1.fffeda7783bedp-1, 0x1.fffffffffffffp-55},
    {-0x1.1fa316819f30fp-22, 0x1.ffffeb4d87f8fp-1, -0x1.fffffffffffffp-55},
    {-0x1.f38d3497c750ap-23, 0x1.ffffee06f60b7p-1, 0x1.fffffffffffffp-55},
    {-0x1.0892e5bab3afbp-29, 0x1.ffffffd9ecbd9p-1, -0x1.fffffffffffffp-55},
    {-0x1.4d89c4fdd2bcbp-54, 0x1.fffffffffffffp-1, -0x1.fffffffffffffp-55},
    {0x1.bcb7b1526e50ep-55, 0x1.0000000000001p+0, -0x1.fffffffffffffp-54},
    {0x1.69554012f9a17p-51, 0x1.0000000000007p+0, -0x1.fffffffffffffp-54},
    {0x1.3fa407733f49bp-50, 0x1.000000000000bp+0, 0x1.fffffffffffffp-54},
    {0x1.df760b2cdeed3p-49, 0x1.0000000000023p+0, -0x1.fffffffffffffp-54},
    {0x1.16d12aac2e244p-46, 0x1.00000000000a1p+0, -0x1.fffffffffffffp-54},
    {0x1.0fc218301ba95p-33, 0x1.0000000138df9p+0, -0x1.fffffffffffffp-54},
    {0x1.58f219ee7f4ecp-33, 0x1.000000018d223p+0, -0x1.fffffffffffffp-54},
    {0x1.b572206147992p-30, 0x1.0000000fbd065p+0, 0x1.fffffffffffffp-54},
    {0x1.f7a045885e2p-27, 0x1.00000090f48a5p+0, -0x1.fffffffffffffp-54},
    {0x1.ca3d76530c29bp-17, 0x1.00020f93a20a5p+0, -0x1.fffffffffffffp-54},
    {0x1.28808b27db928p-14, 0x1.000aab1b954dbp+0, 0x1.135ebda15efc8p-107},
    {0x1.a4e071be911a1p-14, 0x1.000f24dda87bbp+0, 0x1.77f4c61f053ffp-106},
    {0x1.33dcea095e035p-12, 0x1.002c51ecc06a7p+0, 0x1.0c8757bf7abfep-106},
    {0x1.c8d5d475181c2p-11, 0x1.00839ea838256p+0, 0x1.2d336fe64a2b9p-110},
    {0x1.de9217d71bdbfp-11, 0x1.0089e3794caf8p+0, 0x1.01e0a65f9dff7p-105},
    {0x1.7c3ddd23ac8cap-10, 0x1.00db40291e4f5p+0, 0x1.fffffffffffffp-54},
    {0x1.78d2f9978e8bep-9, 0x1.01b3461fe368ep+0, 0x1.94f3a6eb57b23p-109},
    {0x1.ec65645edc394p-8, 0x1.0477a5e0c9213p+0, -0x1.a72aa52fb78bcp-107},
    {0x1.4425c74e7eccdp-7, 0x1.05e5e2a071633p+0, -0x1.faa4410533583p-106},
    {0x1.90d7373b3a546p-7, 0x1.07502eb858f1fp+0, 0x1.fffffffffffffp-54},
    {0x1.03e1f649d39bap-6, 0x1.0985db3c7a556p+0, 0x1.a5b2198aeb0abp-105},
    {0x1.7e3c84f2cb9b5p-6, 0x1.0e20cd9bd42a9p+0, 0x1.fffffffffffffp-54},
    {0x1.10be9f41719bp-5, 0x1.1465bd5d90697p+0, -0x1.97b7b0ec13b05p-105},
    {0x1.25765968ecd68p-5, 0x1.1602fbb34d8aep+0, -0x1.d9938e11bf839p-106},
    {0x1.9aa6fd4d21a47p-5, 0x1.1f525a9f09337p+0, 0x1.30dbd33a0ca8cp-108},
    {0x1.e7b525705edefp-5, 0x1.259ccebe317bfp+0, -0x1.1f35737b5e004p-109},
    {0x1.1a45d0d5b0a79p-3, 0x1.5f9da82fdb9f5p+0, 0x1.ad336373da19dp-105},
    {0x1.e0df7a9a954ccp-3, 0x1.b795554e9eacp+0, -0x1.4c92a870f8e27p-104},
    {0x1.c414aa8bd83b1p-2, 0x1.61bfb59b7b127p+1, 0x1.d706652e05873p-105},
    {0x1.d7d271ab4eeb4p-2, 0x1.71ce472eb84c8p+1, -0x1.0fa6863d41876p-117},
    {0x1.f1f32b8b01af8p-2, 0x1.882ef00e579a7p+1, -0x1.b0840e77e9a02p-103},
    {0x1.f53be0fe16695p-2, 0x1.8b17140ac1f1bp+1, -0x1.a29f4bce2218ap-103},
    {0x1.1fe5f30572361p-1, 0x1.d33582532e5bp+1, 0x1.b2dc31ea1fc54p-106},
    {0x1.8eab1f62d8e8dp-1, 0x1.8070cd731f577p+2, -0x1.aafc6c221a723p-104},
    {0x1.e107654fb6916p-1, 0x1.1664ceca10f9ap+3, -0x1.d9dce9d4d6d4cp-103},
    {0x1p+0, 0x1.4p+3, 0x0p+0},
    {0x1.154b3463f846p+0, 0x1.838d08649702dp+3, 0x1.553c239d8c386p-103},
    {0x1.1daf94cf0bd01p+0, 0x1.a1efc3ce340d2p+3, 0x1.0564e2abb55fbp-102},
    {0x1.75f49c6ad3badp+0, 0x1.ce41d8fa665fap+4, -0x1.539027cc91e93p-114},
    {0x1.a3c782d4f54fcp+0, 0x1.5d05adfa6a8a1p+5, 0x1.c0afe5567f632p-101},
    {0x1.cc30b915ec8c4p+0, 0x1.f60165d5bc3e1p+5, 0x1.165e71f25cbfep-102},
    {0x1p+1, 0x1.9p+6, 0x0p+0},
    {0x1.8p+1, 0x1.f4p+9, 0x0p+0},
    {0x1.ee9674267e65fp+1, 0x1.c8ed39b9d8a37p+12, -0x1.d6a2e7398745cp-97},
    {0x1p+2, 0x1.388p+13, 0x0p+0},
    {0x1.2d5494eb1dd13p+2, 0x1.8f189a48a6a87p+15, -0x1.fffffffffffffp-39},
    {0x1.4p+2, 0x1.86ap+16, 0x0p+0},
    {0x1.60266cc425d1p+2, 0x1.367d049688476p+18, -0x1.d1269fa7b2c2ep-87},
    {0x1.8p+2, 0x1.e848p+19, 0x0p+0},
    {0x1.cp+2, 0x1.312dp+23, 0x0p+0},
    {0x1p+3, 0x1.7d784p+26, 0x0p+0},
    {0x1.2p+3, 0x1.dcd65p+29, 0x0p+0},
    {0x1.4p+3, 0x1.2a05f2p+33, 0x0p+0},
    {0x1.6p+3, 0x1.74876e8p+36, 0x0p+0},
    {0x1.8p+3, 0x1.d1a94a2p+39, 0x0p+0},
    {0x1.ap+3, 0x1.2309ce54p+43, 0x0p+0},
    {0x1.cp+3, 0x1.6bcc41e9p+46, 0x0p+0},
    {0x1.ep+3, 0x1.c6bf52634p+49, 0x0p+0},
    {0x1p+4, 0x1.1c37937e08p+53, 0x0p+0},
    {0x1.1p+4, 0x1.6345785d8ap+56, 0x0p+0},
    {0x1.2p+4, 0x1.bc16d674ec8p+59, 0x0p+0},
    {0x1.3p+4, 0x1.158e460913dp+63, 0x0p+0},
    {0x1.4p+4, 0x1.5af1d78b58c4p+66, 0x0p+0},
    {0x1.5p+4, 0x1.b1ae4d6e2ef5p+69, 0x0p+0},
    {0x1.6p+4, 0x1.0f0cf064dd592p+73, 0x0p+0},
    {0x1.7p+4, 0x1.52d02c7e14af6p+76, 0x1p+23},
    {0x1.89063309f3004p+4, 0x1.83fe31fe6a9d6p+81, 0x1.48902817c6b41p-24},
    {0x1.2a59b82b6fc5ep+6, 0x1.b5f92b2a65f97p+247, -0x1.65b160bda474fp+141},
    {0x1.9bc658ed53988p+6, 0x1.f601dfa307562p+341, 0x1.e0445076e7a83p+236},
    {0x1.f5b1d27635637p+6, 0x1.914044aa9cfc8p+416, 0x1.5cf990ba3c0f9p+312},
    {0x1.f9b1d27635637p+6, 0x1.f59055d5443bap+419, 0x1.b437f4e8cb138p+315},
    {0x1.464c8348af94ap+7, 0x1.f5b55de961a8ep+541, 0x1.1211abe3b7492p+437},
    {0x1.cbe37694f4d1p+7, 0x1.d01b2ef68a124p+763, 0x1.91f3dd2aff9d7p+657},
    {0x1.cde37694f4d1p+7, 0x1.2210fd5a164b7p+767, -0x1.fffffffffffffp+713},
    {0x1.2999c72e3120dp+8, 0x1.86361271301e9p+988, -0x1.6b6e03f7f22bp+884},
  };
  int a, b, c;
  for (a = 0, b = EXCEPTIONS; a + 1 != b;)
  {
    /* Invariant: if x is an exceptional case, we have
       exceptions[a][0] <= x and (x < exceptions[b][0] or b = EXCEPTIONS) */
    c = (a + b) / 2;
    if (exceptions[c][0] <= x)
      a = c;
    else
      b = c;
  }
  if (x == exceptions[a][0])
    return exceptions[a][1] + exceptions[a][2];
  double eh, el;
#define INVLOG2_10 0x1.a934f0979a371p+13 // 2^12*log(10)/log(2)
  double k = __builtin_roundeven (x * INVLOG2_10); // -4399104 <= k <= 4194304
  if (__builtin_expect (k == 4194304, 0))
    k = 4194303; // ensures M < 2047 below

  double yh, yl;
  // LOG2H+LOG2_10M+LOG2L approximates log(2)/log(10)/2^12
#define LOG2_10H 0x1.34413508p-14
#define LOG2_10M 0x1.f79fef3p-46
#define LOG2_10L 0x1.1f12b35816f92p-78
  yh = __builtin_fma (-k, LOG2_10H, x); // exact
  fast_two_sum (&yh, &yl, yh, -k * LOG2_10M);
  yl = __builtin_fma (-k, LOG2_10L, yl);
  // now multiply yh+yl by log(10)
  s_mul (&eh, &el, yh, LOG10H, LOG10L);
  yh = eh;
  yl = __builtin_fma (yl, LOG10H, el);

  int64_t K = k; /* Note: k is an integer, this is just a conversion. */
  int64_t M = (K >> 12) + 0x3ff; // -51 <= M <= 2046
  int64_t i2 = (K >> 6) & 0x3f;
  int64_t i1 = K & 0x3f;

  double t1h = T1[i2][0], t1l = T1[i2][1], t2h = T2[i1][0], t2l = T2[i1][1];
  d_mul (&eh, &el, t2h, t2l, t1h, t1l);

  double qh, ql;
  q_2 (&qh, &ql, yh, yl);

  d_mul (&eh, &el, eh, el, qh, ql);
  f64_u _d;
  
  if (__builtin_expect (M >= 1, 0)) // normal case
  {
    _d.u = M << 52;
    return __builtin_fma (el, _d.f, _d.f * eh);
  }

  // -51 <= M <= 0: subnormal case
  // M=0 corresponds to 2^-1023, M=-51 corresponds to 2^-1074
  // scale eh, el by 2^52 to get out of the subnormal range
  _d.u = (M + 52) << 52;
  fast_two_sum (&yh, &yl, _d.f * eh, _d.f * el);
  double res = yh * 0x1p-52;
  double corr = __builtin_fma (-res, 0x1p52, yh);
  corr += yl;
  /* now add corr*0x1p-52 to res */
  return __builtin_fma (corr, 0x1p-52, res);
}

double cr_exp10 (double x)
{
  b64u64_u t = {.f = x};
  uint64_t ax = t.u & (~0ul>>1);
  if (__builtin_expect (ax >= 0x40734413509f79fful, 0))
    // x = NaN or |x| >= 0x1.34413509f79ffp+8
  {
    if (ax > 0x7ff0000000000000ul)
      return x; // NaN
    if (x >= 0x1.34413509f79ffp+8) /* 10^x > 2^1024*(1-2^-54) */
      return 0x1p1023 + 0x1p1023;
    if (x <= -0x1.439b746e36b53p+8) /* 10^x < 2^-1075 */
      return 0x1p-1074 * 0.5;
    if (x <= -0x1.434e6420f4374p+8) /* 2^-1075 < 10^x < 2^-1074 */
      return 0x3p-1074 * 0.25;
  }
  else if (__builtin_expect (ax <= 0x3c7bcb7b1526e50eul, 0))
    // |x| <= 0x1.bcb7b1526e50ep-56
    return 1 + x; /* for |x| <= -0x1.bcb7b1526e50ep-56, exp10(x) rounds to
                     1 to nearest */
  /* now -0x1.434e6420f4374p+8 < x < -0x1.bcb7b1526e50ep-56
     or 0x1.bcb7b1526e50ep-56 < x < 0x1.34413509f79ffp+8 */

  double h, l, err;
  err = exp10_fast (&h, &l, x);
  double left =  h + (l - err);
  double right = h + (l + err);
  if (left == right)
    return left;

  return exp10_accurate (x);
}
