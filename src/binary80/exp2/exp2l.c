/* Correctly rounded exp2l function for binary64 values.

Copyright (c) 2024 Paul Zimmermann

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
#include <assert.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

// anonymous structs, see https://port70.net/~nsz/c/c11/n1570.html#6.7.2.1p19
typedef union {long double f; struct {uint64_t m; uint16_t e;};} b80u80_t;

/* s + t <- a + b, assuming |a| >= |b| */
static inline void
fast_two_sum (long double *s, long double *t, long double a, long double b)
{
  *s = a + b;
  long double e = *s - a;
  *t = b - e;
}

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void
a_mul (long double *hi, long double *lo, long double a, long double b) {
  *hi = a * b;
  *lo = __builtin_fmal (a, b, -*hi);
}

// Returns (ah + al) * (bh + bl) - (al * bl)
static inline void
d_mul (long double *hi, long double *lo, long double ah, long double al,
       long double bh, long double bl) {
  a_mul (hi, lo, ah, bh);
  *lo = __builtin_fmal (ah, bl, *lo);
  *lo = __builtin_fmal (al, bh, *lo);
}

// T2[i] approximates 2^(i/2^5) with relative error < 2^-129.565
static const long double T2[32][2] = {
   {0x1p+0L, 0x0p+0L},
   {0x1.059b0d31585743aep+0L, 0x1.f1523ada32905ffap-66L},
   {0x1.0b5586cf9890f62ap+0L, -0x1.d1b5239ef559f27p-66L},
   {0x1.11301d0125b50a4ep+0L, 0x1.77e35db26319d58cp-65L},
   {0x1.172b83c7d517adcep+0L, -0x1.06e75e29d6b0dbfap-69L},
   {0x1.1d4873168b9aa78p+0L, 0x1.6e00a2643c1ea62ep-66L},
   {0x1.2387a6e75623866cp+0L, 0x1.fadb1c15cb593b04p-68L},
   {0x1.29e9df51fdee12c2p+0L, 0x1.7457d6892a8ef2a2p-66L},
   {0x1.306fe0a31b7152dep+0L, 0x1.1ab48c60b90bdbdap-65L},
   {0x1.371a7373aa9caa72p+0L, -0x1.755fa17570cf0384p-65L},
   {0x1.3dea64c12342235cp+0L, -0x1.7dbb83d8511808bap-65L},
   {0x1.44e086061892d032p+0L, -0x1.9217ec41fcc08562p-65L},
   {0x1.4bfdad5362a271d4p+0L, 0x1.cbd7f621710701b2p-67L},
   {0x1.5342b569d4f81dfp+0L, 0x1.507893b0d4c7e9ccp-65L},
   {0x1.5ab07dd48542958cp+0L, 0x1.2602a323d668bb12p-65L},
   {0x1.6247eb03a5584b2p+0L, -0x1.e0bf205a4b7a89c6p-65L},
   {0x1.6a09e667f3bcc908p+0L, 0x1.65f626cdd52afa7cp-65L},
   {0x1.71f75e8ec5f73dd2p+0L, 0x1.b879778566b65a1ap-67L},
   {0x1.7a11473eb0186d7ep+0L, -0x1.5dfb81264bc14218p-65L},
   {0x1.82589994cce128acp+0L, 0x1.f115f56694021ed6p-65L},
   {0x1.8ace5422aa0db5bap+0L, 0x1.f156864b26ecf9bcp-66L},
   {0x1.93737b0cdc5e4f46p+0L, -0x1.fc781b57ebba5a08p-65L},
   {0x1.9c49182a3f0901c8p+0L, -0x1.dca7c706a0d3912ap-67L},
   {0x1.a5503b23e255c8b4p+0L, 0x1.2248e57c3de40286p-67L},
   {0x1.ae89f995ad3ad5e8p+0L, 0x1.cd345dcc8169fefp-66L},
   {0x1.b7f76f2fb5e46eaap+0L, 0x1.ec206ad4f14d5322p-66L},
   {0x1.c199bdd85529c222p+0L, 0x1.9625412374ccf288p-69L},
   {0x1.cb720dcef906915p+0L, 0x1.e5e8f4a4edbb0ecap-67L},
   {0x1.d5818dcfba48725ep+0L, -0x1.7e9452647c8d582ap-66L},
   {0x1.dfc97337b9b5eb96p+0L, 0x1.195873da5236e44cp-65L},
   {0x1.ea4afa2a490d9858p+0L, 0x1.ee7431ebb6603f0ep-65L},
   {0x1.f50765b6e4540674p+0L, 0x1.f096ec50c575ff32p-65L},
};

// T1[i] approximates 2^(i/2^10) with relative error < 2^-129.048
static const long double T1[32][2] = {
   {0x1p+0L, 0x0p+0L},
   {0x1.002c605e2e8cec5p+0L, 0x1.b486ff22688e8042p-66L},
   {0x1.0058c86da1c09ea2p+0L, -0x1.cc5ad661a130c72ep-73L},
   {0x1.0085382faef831dap+0L, 0x1.27f2106beea70f16p-65L},
   {0x1.00b1afa5abcbed62p+0L, -0x1.aca9d827dc46d578p-65L},
   {0x1.00de2ed0ee0f4f6p+0L, -0x1.ab13a069914e78d8p-67L},
   {0x1.010ab5b2cbd11708p+0L, -0x1.7ccfd6d8fbc56654p-65L},
   {0x1.0137444c9b5b4ed4p+0L, 0x1.2a293d12edc0f6d8p-65L},
   {0x1.0163da9fb33356d8p+0L, 0x1.299ab8cdb737e9p-66L},
   {0x1.019078ad6a19efp+0L, -0x1.1dfea2857f2adcfap-65L},
   {0x1.01bd1e77170b415ep+0L, 0x1.d899887ad6abfd84p-66L},
   {0x1.01e9cbfe113eec7ep+0L, -0x1.f5239bf535594f58p-67L},
   {0x1.02168143b0280da8p+0L, 0x1.9de0756294cca9f6p-68L},
   {0x1.02433e494b754b3ap+0L, 0x1.aaf4ec3aae71c11ep-65L},
   {0x1.027003103b10def8p+0L, -0x1.77a8db0ebeced796p-67L},
   {0x1.029ccf99d720a05ap+0L, -0x1.9a22b7e9aec548fp-65L},
   {0x1.02c9a3e778060ee6p+0L, 0x1.ef95949ef4537bd2p-65L},
   {0x1.02f67ffa765e5c8ep+0L, 0x1.278b1213c0c9e1b6p-66L},
   {0x1.032363d42b0277fap+0L, 0x1.46b0f6b00b29401ep-65L},
   {0x1.03504f75ef0716fp+0L, 0x1.77472bccd623cb4ap-65L},
   {0x1.037d42e11bbcc0acp+0L, -0x1.7ee11521ee5bb3bp-65L},
   {0x1.03aa3e170aafd83ap+0L, -0x1.49e0dc1269659b0ep-65L},
   {0x1.03d7411915a8a6ep+0L, -0x1.ff8c2457133e5c34p-65L},
   {0x1.04044be896ab6678p+0L, -0x1.e2913831fef18048p-65L},
   {0x1.04315e86e7f84bd8p+0L, -0x1.8e0cbbe4b703226p-65L},
   {0x1.045e78f5640b9136p+0L, -0x1.0de542c45976151ep-66L},
   {0x1.048b9b35659d809p+0L, 0x1.cd53d5e8b6609244p-65L},
   {0x1.04b8c54847a27e18p+0L, 0x1.9b69feee140b2d6cp-66L},
   {0x1.04e5f72f654b1298p+0L, 0x1.bc9d50684640c7dap-66L},
   {0x1.051330ec1a03f5e6p+0L, 0x1.45f11ce522be682ep-65L},
   {0x1.0540727fc176195p+0L, 0x1.a8eda3f31093fe7cp-65L},
   {0x1.056dbbebb786b20ep+0L, -0x1.bc5e64449ba34522p-66L},
};

// T0[i] approximates 2^(i/2^15) with relative error < 2^-129.004
static const long double T0[32][2] = {
   {0x1p+0L, 0x0p+0L},
   {0x1.000162e525ee0548p+0L, -0x1.5775054cd5adbfb2p-65L},
   {0x1.0002c5cc37da9492p+0L, -0x1.7b3d1e5b9cb8c262p-67L},
   {0x1.000428b535c857eep+0L, -0x1.64c2b3ef9bd797e4p-67L},
   {0x1.00058ba01fb9f96ep+0L, -0x1.26a6569cfedd0784p-65L},
   {0x1.0006ee8cf5b22326p+0L, 0x1.77c33a014414bc8ep-66L},
   {0x1.0008517bb7b37f32p+0L, 0x1.a5127d0b5ff94c8cp-68L},
   {0x1.0009b46c65c0b7aep+0L, -0x1.cbd67bc2e9bcfbf6p-67L},
   {0x1.000b175effdc76bap+0L, 0x1.c718b38e549cb934p-67L},
   {0x1.000c7a538609667cp+0L, -0x1.a6bef4105b137bf2p-70L},
   {0x1.000ddd49f84a311cp+0L, -0x1.7ca37fadb538a1d8p-65L},
   {0x1.000f404256a180c4p+0L, -0x1.77da3c7a168d87dap-71L},
   {0x1.0010a33ca111ffa6p+0L, -0x1.bb7ff655871c632cp-67L},
   {0x1.00120638d79e57f4p+0L, -0x1.a2a9629bed7b0238p-69L},
   {0x1.00136936fa4933e6p+0L, -0x1.06c95b8aba5aab5ep-65L},
   {0x1.0014cc3709153db6p+0L, -0x1.d04108b0bf2a604p-65L},
   {0x1.00162f3904051fa2p+0L, -0x1.ae86ac75479c344p-65L},
   {0x1.0017923ceb1b83ecp+0L, -0x1.d7d5ad2426d98758p-67L},
   {0x1.0018f542be5b14dap+0L, 0x1.68bbbf240fe795acp-65L},
   {0x1.001a584a7dc67cb8p+0L, -0x1.1f7bc1b6df8284a4p-65L},
   {0x1.001bbb54296065dp+0L, -0x1.b7ebdcc748e85934p-65L},
   {0x1.001d1e5fc12b7a72p+0L, 0x1.59de63237804a4cep-65L},
   {0x1.001e816d452a64f6p+0L, 0x1.342315e8f1e6f0fap-65L},
   {0x1.001fe47cb55fcfb4p+0L, -0x1.a7064b4959898e28p-65L},
   {0x1.0021478e11ce6504p+0L, 0x1.5cb6b16a8e0ad03cp-66L},
   {0x1.0022aaa15a78cf4ap+0L, -0x1.03a60c77b646fde4p-66L},
   {0x1.00240db68f61b8e6p+0L, 0x1.7649ad42d581bc88p-65L},
   {0x1.002570cdb08bcc42p+0L, 0x1.5119c9d215fbae7p-66L},
   {0x1.0026d3e6bdf9b3c8p+0L, -0x1.752352535fcc167ep-65L},
   {0x1.00283701b7ae19e4p+0L, -0x1.19822944d4228146p-70L},
   {0x1.00299a1e9daba90ap+0L, 0x1.2baca861d8c8d1f4p-65L},
   {0x1.002afd3d6ff50bbp+0L, 0x1.ca8a335347ceeba2p-65L},
};

// put in h+l an approximation of 2^x for |x| < 2^-16, with relative error
// bounded by 2^-83.747 (see routine analyze_P in exp2l.sage)
static void
P (long double *h, long double *l, long double x)
{
  /* the following degree-4 polynomial generated by exp2.sollya has absolute
     error bounded by 2^-83.748 for |x| < 2^-16 */
  static const long double p[] = {1.0L, 0x1.62e42fefa39ef358p-1L,
                            0x1.ebfbdff82c58ea86p-3L, 0x1.c6b08d6835c26dep-5L,
                            0x1.3b2ab70cf131bd7ep-7L};
  long double y = __builtin_fmal (p[4], x, p[3]);
  y = __builtin_fmal (y, x, p[2]);
  fast_two_sum (h, l, p[1], y * x);
  long double t;
  a_mul (h, &t, *h, x);
  t = __builtin_fmal (*l, x, t);
  fast_two_sum (h, l, p[0], *h);
  *l += t;
}

#define TRACE 0x8.0017f89c43248p+6L

/* Assume -16446 < x < -0x1.71547652b82fe176p-65
   or 0x1.71547652b82fe176p-64 < x < 16384.
   Return h + l approximating 2^x with relative error < 2^-83.727
   or h = l = NaN.
*/
static void
fast_path (long double *h, long double *l, long double x)
{
  // if (x == TRACE) printf ("enter fast_path x=%La\n", x);
  b80u80_t v = {.f = x};
  
  int32_t k = __builtin_roundl (0x1p15L * x); // -16445*2^15 <= k <= 16383*2^15
  // if (x == TRACE) printf ("k=%d\n", k);
  long double r = __builtin_fmal (k, -0x1p-15L, x);
  // if (x == TRACE) printf ("r=%La\n", r);
  int32_t i = (k + 538869760) & 32767;
  // if (x == TRACE) printf ("i=%d\n", i);
  int32_t e = (k - i) >> 15;
  // if (x == TRACE) printf ("e=%d\n", e);
  int32_t i0 = i & 0x1f, i1 = (i >> 5) & 0x1f, i2 = i >> 10;
  // if (x == TRACE) printf ("i2=%d i1=%d i0=%d\n", i2, i1, i0);
  // k = e*2^15 + i2*2^10 + i1*2^5 + i0
  // x = k*2^-15 + r with |r| < 2^-16
  // 2^x = 2^e * 2^(i2/2^5) * 2^(i1/2^10) * 2^(i0/2^15) * 2^r
  P (h, l, r); // relative error bounded by 2^-83.747
  // if (x == TRACE) printf ("P: h=%La l=%La\n", *h, *l);
  long double hh, ll;
  d_mul (&hh, &ll, T2[i2][0], T2[i2][1], T1[i1][0], T1[i1][1]);
  /* The d_mul() call decomposes into with ah=T2[i2][0], al=T2[i2][1],
     bh=T1[i1][0], bl=T1[i1][1]:
     a_mul (hh, ll0, ah, bh)
     ll1 = __builtin_fmal (ah, bl, ll0)
     ll = __builtin_fmal (al, bh, ll1)
     Since ah <= 0x1.f50765b6e4540674p+0 and bh <= 0x1.056dbbebb786b20ep+0,
     we have ah*bh < 2, thus |ll0| <= ulp(1) = 2^-63.
     Since |ah| < 2 and |bl| < 2^-64, we have |ll1| <= |ah*bl+ll0| < 2^-62,
     thus the error in the first fma() is bounded by ulp(2^-62-eps) = 2^-126.
     Since |al| < 2^-64 and |bh| < 2, we have |ll| <= |al*bh+ll1| < 2^-62+2^-63
     thus the error in the 2nd fma() is bounded by ulp(2^-62+2^-63) = 2^-125.
     The ignored term al*bl is bounded by 2^-128, which gives a maximal
     absolute error for this d_mul() call of:
     2^-126 + 2^-125 + 2^-128 < 2^-124.299.
     The same bound holds for the relative error since each of ah+al and
     bh+bl is >= 1, thus hh+ll >= 1. */
  d_mul (&hh, &ll, hh, ll, T0[i0][0], T0[i0][1]);
  // if (x == TRACE) printf ("hh=%La ll=%La\n", hh, ll);
  /* The d_mul() call decomposes into with ah=hh_in, al=ll_in,
     bh=T0[i0][0], bl=T0[i0][1]:
     a_mul (hh, ll0, ah, bh)
     ll1 = __builtin_fmal (ah, bl, ll0)
     ll = __builtin_fmal (al, bh, ll1)
     Since ah <= 0x1.ffa74ea381efc218p+0 (this can be proven by multiplying
     T2[31][0] and T1[31][0] with rounding up) and bh<=0x1.002afd3d6ff50bbp+0,
     we have ah*bh < 2, thus |ll0| <= ulp(1) = 2^-63.
     Since |ah| < 2 and |bl| < 2^-64, we have |ll1| <= |ah*bl+ll0| < 2^-62,
     thus the error in the first fma() is bounded by ulp(2^-62-eps) = 2^-126.
     Since |al| < 2^-62+2^-63 (bound for ll in the first d_mul() call) and
     |bh| < 2, we have |ll| <= |al*bh+ll1| < 2^-60,
     thus the error in the 2nd fma() is bounded by ulp(2^-60-eps) = 2^-124.
     The ignored term al*bl is bounded by 2^-125, which gives a maximal
     absolute error for this d_mul() call of:
     2^-126 + 2^-124 + 2^-125 < 2^-123.192.
     The same bound holds for the relative error since each of ah+al and
     bh+bl is >= 1, thus hh+ll >= 1. */

  /* After the two d_mul() calls, the relative error is bounded by,
     taking also into account the errors in T2, T1, T0:
     |hh + ll - 2^(i2/2^5+i1/2^10+i0/2^15)| < (1 + err) * (hh + ll)
     where err < (1 + 2^-124.299) * (1 + 2^-123.192) * (1 + 2^-129.565) *
     (1 + 2^-129.048) * (1 + 2^-129.004) - 1 < 2^-122.596. */

  d_mul (h, l, *h, *l, hh, ll);
  /* The d_mul() call decomposes into with ah=h_in, al=l_in, bh=hh, bl=ll:
     a_mul (h, l0, ah, bh)
     l1 = __builtin_fmal (ah, bl, l0)
     l = __builtin_fmal (al, bh, l1)
     Since ah <= 2^(2^-16) <= 1.00002 and bh <= 1.99996 (this can be seen
     by considering the i2=i1=i0=31 above), we have ah*bh < 2, thus
     |l0| <= ulp(1) = 2^-63.
     Since |ah| < 2 and |bl| < 2^-60 (see above), we have |l1| <= |ah*bl+l0|
     < 2^-58, thus the error in the first fma() is bounded by ulp(2^-58-eps)
     = 2^-122.
     Since |al| < 2^-63 (see analyze_P()) and |bh| < 2, we have
     |l| <= |al*bh+l1| < 2^-62+2^-58
     thus the error in the 2nd fma() is bounded by ulp(2^-62+2^-58) = 2^-121.
     The ignored term al*bl is bounded by 2^-123, which gives a maximal
     absolute error for this d_mul() call of:
     2^-122 + 2^-121 + 2^-123 < 2^-120.192.
     For the relative error, we might have h_in = 2^(-2^-16), which gives a
     bound of 2^-120.192/2^(-2^-16) < 2^-120.191. */

  /* The relative errors are:
   * that on h_in + l_in, bounded by 2^-83.747
   * that on hh + ll, bounded by 2^-122.596
   * that from this last d_mul() call, bounded by 2^-120.192
   This gives a final relative error:
   |h + l - 2^(i2/2^5+i1/2^10+i0/2^15+r)| < (1 + 2^-83.747) * (1 + 2^-122.596)
   * (1 + 2^-120.191) - 1 < 2^-83.746, with 0.99998 <= h + l < 2. */

  if (__builtin_expect (e >= -16355, 1))
  {
    /* Multiply h, l by 2^e. Since e >= -16355, we have 2^x>=0.99998*2^-16355
       thus if l*2^e is in the subnormal range, we have an additional absolute
       error of at most 2^-16445, which corresponds to an additional relative
       error < 2^-16445/(0.99998*2^-16355) < 2^-89.999. This gives a final
       bound of (1 + 2^-83.746) * (1 + 2^-89.999) - 1 < 2^-83.727.
       No overflow is possible here since x < 16384. */
    v.f = *h;
    v.e += e;
    *h = v.f;
    b80u80_t w = {.f = *l};
    w.e += e;
    *l = w.f;
  }
  else
  {
    v.e = 32767;                                      
    v.m = 0xc000000000000000ul;
    *h = *l = v.f; // +qnan
  }
}

long double
cr_exp2l (long double x)
{
  b80u80_t v = {.f = x};
  uint16_t e = v.e & 0x7fff;

  // printf ("x=%La v.e=%u\n", x, v.e);

  // check NaN, Inf, overflow, underflow
  // overflow for x >= 16384, i.e., 16397 <= e <= 32767
  // the smallest subnormal is 2^-16445
  if (__builtin_expect (e >= 16397, 0))
  {
    if (e == 0x7fff)
    { // NaN or Inf: 2^x = x for x = NaN or +Inf, 2^-Inf = 0
      if (v.e == 0xffff && v.m == 0x8000000000000000ul) // -Inf
        return 0x0p0L;
      return x;
    }
    if (x >= 0x1p+14L) // x >= 16384
      return 0x1p16383L + 0x1p16383L;
    // now x < 0
    if (x <= -0x1.00f8p+14L) // x <= -16446
      return 0x1p-16445L * 0.5L;
  }

  // case of tiny inputs
  // for 0 <= x <= 0x1.71547652b82fe176p-64, 2^x rounds to 1 to nearest
  // for -0x1.71547652b82fe176p-65 <= x <= 0, 2^x rounds to 1 to nearest
  if (__builtin_expect (e <= 16319, 0))
  {
    if (0 <= x && x <= 0x1.71547652b82fe176p-64L)
      return __builtin_fmal (x, x, 0x1p0L);
    if (-0x1.71547652b82fe176p-65L <= x && x < 0)
      return __builtin_fmal (x, -x, 0x1p0L);
  }

  long double h, l;
  fast_path (&h, &l, x);
  static const long double err = 0x1.36p-84; // 2^-83.727 < err
  long double left = h +  __builtin_fmal (h, -err, l);
  long double right = h + __builtin_fmal (h, err, l);
  // if (x == TRACE) printf ("left=%La right=%La\n", left, right);
  if (__builtin_expect (left == right, 1))
    return left;

  return -1;
}
