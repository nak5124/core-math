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

// Veltkamp's splitting: split x into xh + xl such that
// x = xh + xl exactly
// xh fits in 32 bits and |xh| <= 2^e if 2^(e-1) <= |x| < 2^e
// xl fits in 32 bits and |xl| < 2^(e-32)
static inline void
split (long double *xh, long double *xl, long double x)
{
  static const long double C = 0x1.00000001p+32L;
  long double gamma = C * x;
  long double delta = x - gamma;
  *xh = gamma + delta;
  *xl = x - *xh;
}

// Dekker's algorithm: rh + rl = u * v
// Reference: Algorithm Mul12 from https://ens-lyon.hal.science/ensl-01529804 pages 21-22
// See also Handbook of Floating-Point Arithmetic, 2nd edition, Veltkamp splitting (Algorith 4.9)
// and Dekker's product (Algorithm 4.10)
// The Handbook only mentions rounding to nearest, but exhaustive tests up to precision 10
// seem to indicate it also works for directed roundings.
// This is confirmed by "Note on the Veltkamp/Dekker Algorithms with Directed Roundings",
// Paul Zimmermann, February 2024 (assuming no underflow/overflow).
static inline void
a_mul (long double *rh, long double *rl, long double u, long double v)
{
  long double u1, u2, v1, v2;
  split (&u1, &u2, u);
  split (&v1, &v2, v);
  *rh = u * v;
  *rl = (((u1 * v1 - *rh) + u1 * v2) + u2 * v1) + u2 * v2;
}

// Return in hi+lo a 128-bit approximation of (ah + al) * (bh + bl)
static inline void
d_mul (long double *hi, long double *lo, long double ah, long double al,
       long double bh, long double bl) {
  a_mul (hi, lo, ah, bh); // exact
  *lo += ah * bl;
  *lo += al * bh;
}

/* Return in hi+lo a 96-bit approximation of (ah + al) * (bh + bl), assuming
   1 <= ah+al, bh+bl < 2. */
static inline void
d_mul1 (long double *hi, long double *lo, long double ah, long double al,
        long double bh, long double bl) {
  static const long double C = 0x1.8p+32l;
  long double ahh = (C + ah) - C, bhh = (C + bh) - C;
  long double ahl = ah - ahh, bhl = bh - bhh;
  *hi = ahh * bhh; // exact since ahh and bhh have at most 32 significant bits
  long double t1 = ahh * (bhl + bl);
  long double t2 = (ahl + al) * bhh;
  long double t3 = (ahl + al) * (bhl + bl);
  *lo = t1 + (t2 + t3);
}

// Same as d_mul1, but assumes ah and bh fit into 32 bits
static inline void
d_mul2 (long double *hi, long double *lo, long double ah, long double al,
        long double bh, long double bl) {
  *hi = ah * bh; // exact
  long double t1 = ah * bl, t2 = al * bh, t3 = al * bl;
  *lo = (t1 + t2) + t3;
}

// Same as d_mul1, but assumes bh fits into 32 bits
static inline void
d_mul3 (long double *hi, long double *lo, long double ah, long double al,
        long double bh, long double bl) {
  static const long double C = 0x1.8p+32l; // ulp(C) = 2^-31
  long double ahh = (C + ah) - C, ahl = ah - ahh;
  *hi = ahh * bh; // exact
  long double t1 = ahh * bl;
  long double t2 = (ahl + al) * bh;
  long double t3 = (ahl + al) * bl;
  *lo = (t1 + t3) + t2;
}

// T2fast[i] approximates 2^(i/2^5) with absolute error < 2^-97.150
static const long double T2fast[32][2] = {
   {0x1p+0L, 0x0p+0L},
   {0x1.059b0d32p+0L, -0x1.4f5178a30756e292p-33L},
   {0x1.0b5586dp+0L, -0x1.9dbc2759d1b5239ep-34L},
   {0x1.11301d02p+0L, -0x1.b495eb62881ca24ep-33L},
   {0x1.172b83c8p+0L, -0x1.5742919041b9d78ap-35L},
   {0x1.1d487316p+0L, 0x1.17354f00b7005132p-33L},
   {0x1.2387a6e8p+0L, -0x1.53b8f327c0a49c7ep-33L},
   {0x1.29e9df52p+0L, -0x1.08f69ed175052edap-39L},
   {0x1.306fe0a4p+0L, -0x1.c91d5a42e54b73ap-33L},
   {0x1.371a7374p+0L, -0x1.558d563aeabf42eap-34L},
   {0x1.3dea64c2p+0L, -0x1.b97bb9497dbb83d8p-33L},
   {0x1.44e08606p+0L, 0x1.892d03136f409dfp-36L},
   {0x1.4bfdad54p+0L, -0x1.3abb1c578d0a0278p-33L},
   {0x1.5342b56ap+0L, -0x1.583f107abe1db13cp-35L},
   {0x1.5ab07dd4p+0L, 0x1.0a852b192602a324p-33L},
   {0x1.6247eb04p+0L, -0x1.6a9ed383c17e40b4p-34L},
   {0x1.6a09e668p+0L, -0x1.8866dee9a09d9322p-37L},
   {0x1.71f75e8ep+0L, 0x1.8bee7ba46e1e5de2p-33L},
   {0x1.7a11473ep+0L, 0x1.6030dafaa2047edap-33L},
   {0x1.82589994p+0L, 0x1.99c25159f115f566p-33L},
   {0x1.8ace5422p+0L, 0x1.541b6b74f8ab4326p-33L},
   {0x1.93737b0cp+0L, 0x1.b8bc9e8a0387e4a8p-33L},
   {0x1.9c49182ap+0L, 0x1.f8480e3e235838fap-35L},
   {0x1.a5503b24p+0L, -0x1.daa374bdbb6e3508p-36L},
   {0x1.ae89f996p+0L, -0x1.4b14a85e32cba234p-34L},
   {0x1.b7f76f3p+0L, -0x1.286e455613df952cp-34L},
   {0x1.c199bdd8p+0L, 0x1.54a7088832c4a824p-34L},
   {0x1.cb720dcep+0L, 0x1.f20d22a0797a3d2ap-33L},
   {0x1.d5818ddp+0L, -0x1.16de36897e945264p-34L},
   {0x1.dfc97338p+0L, -0x1.192851a5cd4f184cp-34L},
   {0x1.ea4afa2ap+0L, 0x1.24366163dce863d8p-34L},
   {0x1.f50765b6p+0L, 0x1.c8a80ce9f096ec5p-33L},
};

// T1fast[i] approximates 2^(i/2^10) with absolute error < 2^-97.024
static const long double T1fast[32][2] = {
   {0x1p+0L, 0x0p+0L},
   {0x1.002c605ep+0L, 0x1.74676283690dfe44p-35L},
   {0x1.0058c86ep+0L, -0x1.78fd85780398b5acp-34L},
   {0x1.0085383p+0L, -0x1.441f3895b01bdf28p-34L},
   {0x1.00b1afa6p+0L, -0x1.50d04a7b5953b05p-34L},
   {0x1.00de2edp+0L, 0x1.dc1e9ebf953b17e6p-33L},
   {0x1.010ab5b2p+0L, 0x1.97a22e0e83302928p-33L},
   {0x1.0137444cp+0L, 0x1.36b69da92a293d12p-33L},
   {0x1.0163daap+0L, -0x1.3332a49ed6654732p-34L},
   {0x1.019078aep+0L, -0x1.2bcc22011dfea286p-33L},
   {0x1.01bd1e78p+0L, -0x1.d1e97d4313b33bc2p-33L},
   {0x1.01e9cbfep+0L, 0x1.13eec7dc15b8c816p-36L},
   {0x1.02168144p+0L, -0x1.3f5fc95f9887e2a8p-34L},
   {0x1.02433e4ap+0L, -0x1.6915698a550b13c6p-33L},
   {0x1.0270031p+0L, 0x1.d886f7be885724f2p-35L},
   {0x1.029ccf9ap+0L, -0x1.46fafd36688adfa6p-35L},
   {0x1.02c9a3e8p+0L, -0x1.0ff3e232106a6b62p-33L},
   {0x1.02f67ffap+0L, 0x1.d9797239278b1214p-34L},
   {0x1.032363d4p+0L, 0x1.5813bfd51ac3dacp-35L},
   {0x1.03504f76p+0L, -0x1.0f8e90f445c6a19ap-36L},
   {0x1.037d42e2p+0L, -0x1.c8867ea97ee11522p-33L},
   {0x1.03aa3e18p+0L, -0x1.eaa04f8d49e0dc12p-33L},
   {0x1.03d7411ap+0L, -0x1.d4aeb241ff8c2458p-33L},
   {0x1.04044be8p+0L, 0x1.2d56ccee1d6ec7cep-33L},
   {0x1.04315e86p+0L, 0x1.cff097ae71f3441cp-33L},
   {0x1.045e78f6p+0L, -0x1.37e8dd9486f2a162p-33L},
   {0x1.048b9b36p+0L, -0x1.34c4fede32ac2a18p-33L},
   {0x1.04b8c548p+0L, 0x1.1e89f8619b69feeep-34L},
   {0x1.04e5f73p+0L, -0x1.3569dacf21b157ccp-33L},
   {0x1.051330ecp+0L, 0x1.a03f5e6a2f88e72ap-36L},
   {0x1.0540728p+0L, -0x1.f44f35795c497034p-35L},
   {0x1.056dbbecp+0L, -0x1.21e537c9bc5e6444p-34L},
};

// T0fast[i] approximates 2^(i/2^15) with absolute error < 2^-97.055
static const long double T0fast[32][2] = {
   {0x1p+0L, 0x0p+0L},
   {0x1.000162e6p+0L, -0x1.b423f5715775054cp-33L},
   {0x1.0002c5ccp+0L, 0x1.bed4a48e84c2e1a4p-35L},
   {0x1.000428b6p+0L, -0x1.946f50245930acfcp-33L},
   {0x1.00058bap+0L, 0x1.fb9f96d6cacd4b18p-36L},
   {0x1.0006ee8cp+0L, 0x1.eb64464cbbe19dp-33L},
   {0x1.0008517cp+0L, -0x1.2132033796bb60bep-34L},
   {0x1.0009b46cp+0L, 0x1.9702deb71a14c21ep-34L},
   {0x1.000b175ep+0L, 0x1.ffb8ed7471c62ce4p-33L},
   {0x1.000c7a54p+0L, -0x1.e7da66101a6bef42p-34L},
   {0x1.000ddd4ap+0L, -0x1.ed73b92f946ff5b6p-38L},
   {0x1.000f4042p+0L, 0x1.5a86030ff4412e1cp-34L},
   {0x1.0010a33cp+0L, 0x1.4223ff4b9120026ap-33L},
   {0x1.00120638p+0L, 0x1.af3cafe7e5d569d6p-33L},
   {0x1.00136936p+0L, 0x1.f49267caf936a476p-33L},
   {0x1.0014cc38p+0L, -0x1.edd58495d04108bp-33L},
   {0x1.00162f3ap+0L, -0x1.f7f5c0bdae86ac76p-33L},
   {0x1.0017923cp+0L, 0x1.d63707d78a0a94b6p-33L},
   {0x1.0018f542p+0L, 0x1.7cb629b568bbbf24p-33L},
   {0x1.001a584ap+0L, 0x1.f719f2ddc1087c92p-34L},
   {0x1.001bbb54p+0L, 0x1.4b032e7920508ce2p-35L},
   {0x1.001d1e6p+0L, -0x1.f6a42c6a98867372p-35L},
   {0x1.001e816ep+0L, -0x1.75ab3612cbdcea18p-33L},
   {0x1.001fe47cp+0L, 0x1.6abf9f6658f9b4b6p-33L},
   {0x1.0021478ep+0L, 0x1.1ce6504572dac5aap-36L},
   {0x1.0022aaa2p+0L, -0x1.4b0e616c81d3063cp-33L},
   {0x1.00240db6p+0L, 0x1.1ec371cd7649ad42p-33L},
   {0x1.002570cep+0L, -0x1.3dd0cef6aee6362ep-34L},
   {0x1.0026d3e6p+0L, 0x1.7bf3678e8adcadacp-33L},
   {0x1.00283702p+0L, -0x1.2147987011982294p-34L},
   {0x1.00299a1ep+0L, 0x1.3b5752152baca862p-33L},
   {0x1.002afd3ep+0L, -0x1.2015e89e3575ccacp-33L},
};

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
// bounded by 2^-78.947 (see routine analyze_P in exp2l.sage), and |l| < 2^-63
static void
P (long double *h, long double *l, long double x)
{
  /* the following degree-4 polynomial generated by exp2.sollya has absolute
     error bounded by 2^-83.748 for |x| < 2^-16 */
  static const long double p[] = {1.0L, 0x1.62e42fefa39ef358p-1L,
                            0x1.ebfbdff82c58ea86p-3L, 0x1.c6b08d6835c26dep-5L,
                            0x1.3b2ab70cf131bd7ep-7L};
  long double y = p[4] * x + p[3];
  y = y * x + p[2];
  *h = y * x + p[1];
  *h = *h * x;
  fast_two_sum (h, l, p[0], *h);
}

// put in h+l an approximation of 2^x for |x| < 2^-16, with relative error
// bounded by 2^-125.403 (see routine analyze_Pacc in exp2l.sage), and |l| < 2^-62.999
static void
Pacc (long double *h, long double *l, long double x)
{
  /* the following degree-6 polynomial generated by exp2acc.sollya has absolute
     error bounded by 2^-133.987 for |x| < 2^-16 */
  static const long double p[] = {1.0L, // degree 0
                                  0x1.62e42fefa39ef358p-1L, -0x1.b0e2633fe0676a9cp-67L, // degree 1
                                  0x1.ebfbdff82c58ea86p-3L, 0x1.e2d60dd936b9ba5ep-68L,  // degree 2
                                  0x1.c6b08d704a0bf8b4p-5L, -0x1.8b4ba2fbcf44117p-70L,  // degree 3
                                  0x1.3b2ab6fba4e7729cp-7L,  // degree 4
                                  0x1.5d87fe78ad725bcep-10L, // degree 5
                                  0x1.4309131bde9fabeap-13L, // degree 6
  };
  long double y = p[9] * x + p[8]; // a6*x+a5
  y = y * x + p[7];                // y*x+a4
  y = y * x;                       // y*x
  fast_two_sum (h, l, p[5], y);    // a3h+y
  *l += p[6];                      // +a3l
  // multiply h+l by x
  long double t;
  a_mul (h, &t, *h, x);            // exact
  *l = *l * x + t;
  // add a2h+a2l
  fast_two_sum (h, &t, p[3], *h);
  *l += t + p[4];
  // multiply h+l by x
  a_mul (h, &t, *h, x);            // exact
  *l = *l * x + t;
  // add a1h+a1l
  fast_two_sum (h, &t, p[1], *h);
  *l += t + p[2];
  // multiply h+l by x
  a_mul (h, &t, *h, x);            // exact
  *l = *l * x + t;
  // add a0
  fast_two_sum (h, &t, p[0], *h);
  *l += t;
}

#define TRACE -0xf.f84cbb09f612fcap+10L

/* Assume -16446 < x < -0x1.71547652b82fe176p-65
   or 0x1.71547652b82fe176p-64 < x < 16384.
   Return h + l approximating 2^x with relative error < 2^-78.891
   or h = l = NaN.
*/
static void
fast_path (long double *h, long double *l, long double x)
{
  // if (x == TRACE) printf ("enter fast_path x=%La\n", x);
  b80u80_t v = {.f = x};

  int32_t k = __builtin_roundl (0x1p15L * x); // -16445*2^15 <= k <= 16383*2^15
  // if (x == TRACE) printf ("k=%d\n", k);
  long double r = x - (long double) k * 0x1p-15L;
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
  P (h, l, r); // relative error bounded by 2^-78.947
  // if (x == TRACE) printf ("P: h=%La l=%La\n", *h, *l);
  long double hh, ll;
  d_mul2 (&hh, &ll, T2fast[i2][0], T2fast[i2][1], T1fast[i1][0], T1fast[i1][1]);
  /* We have |T2fast[i2][0]|, |T1fast[i1][0]| < 2,
     |T2fast[i2][1]|, |T1fast[i1][1]| < 2^-32
     and |hh| < 2:
     * the absolute error on T2fast[i2] is bounded by 2^-97.150, which gives
       an absolute error < 2^-97.150*T1fast[i1] < 2^-97.119
     * the absolute error on T1fast[i1] is bounded by 2^-97.024, which gives
       an absolute error < 2^-97.024*T2fast[i2] < 2^-96.055
     * the 2nd term error is bounded by 2^-97.150*2^-97.024 = 2^-194.174
     * the d_mul2() call decomposes into:
       hi = T2fast[i2][0] * T1fast[i1][0] [exact]
       t1 = T2fast[i2][0] * T1fast[i1][1]
       t2 = T2fast[i2][1] * T1fast[i1][0]
       t3 = T2fast[i2][1] * T1fast[i1][1]
       lo = (t1 + t2) + t3
       Since |T2fast[i2][0]| < 2 and |T1fast[i1][1]| < 2^-32, we have |t1| < 2^-31
       and the rounding error on t1 is bounded by ulp(2^-32) = 2^-95.
       Since |T2fast[i2][1]| < 2^-32 and |T1fast[i1][0]| < 2, we have |t2| < 2^-31
       and the rounding error on t2 is bounded by ulp(2^-32) = 2^-95.
       Since |T2fast[i2][1]| < 2^-32 and |T1fast[i1][1]| < 2^-32, we have |t3| < 2^-64
       and the rounding error on t3 is bounded by ulp(2^-65) = 2^-128.
       Then we have |t1+t2| < 2^-30 and the rounding error on t1+t2 is bounded by 2^-94.
       Then (t1+t2)+t3 < 2^-29 and the rounding error on (t1+t2)+t3 is bounded by 2^-93.
       The absolute error is thus bounded by:
       2^-97.119 + 2^-96.055 + 2^-194.174 + 2^-95 + 2^-95 + 2^-94 + 2^-93 < 2^-91.877:
       | hh + ll - 2^(i2/2^5) * 2^(i1/2^10) | < 2^-92.763
       with |hh| < 2 and |ll| < 2^-29.999.
  */
  d_mul3 (&hh, &ll, hh, ll, T0fast[i0][0], T0fast[i0][1]);
  /* We have |hh_in|, |T0fast[i0][0]| < 2, |ll_in| < 2^-30, |T0fast[i0][1]| < 2^-32:
     * the absolute error on T0fast[i0] is bounded by 2^-97.055, which gives
       an absolute error < (2+2^-30)*2^-97.055 < 2^-96.054
     * the absolute error on hh_in+ll_in is bounded by 2^-92.763, which gives
       an absolute error < 2^-92.763*T0fast[i0] < 2^-92.762
     * the 2nd term error is bounded by 2^-92.763*2^-97.055 = 2^-189.818
     * the d_mul3() call decomposes into:
       hi = ahh * T0fast[i0][0] [exact]
       t1 = ahh * T0fast[i0][1]
       t2 = (ahl + al) * T0fast[i0][0]
       t3 = (ahl + al) * T0fast[i0][1]
       lo = (t1 + t3) + t2
       Since |ahh| < 2 and |T0fast[i0][1]| < 2^-32, we have |t1| < 2^-31
       and the rounding error on t1 is bounded by ulp(2^-32) = 2^-95.
       Since |ahl + al| < 2^-31+2^-31 = 2^-30 and |T0fast[i0][0]| < 2, we have |t2| < 2^-29
       and the rounding error on t2 is bounded by ulp(2^-30) = 2^-93.
       Since |ahl + al| < 2^-31+2^-31 = 2^-30 and |T0fast[i0][1]| < 2^-32, we have |t3| < 2^-62
       and the rounding error on t3 is bounded by ulp(2^-63) = 2^-126.
       Then we have |t1+t3| < 2^-31+2^-62 and the rounding error on t1+t3 is bounded
       by ulp(2^-31+2^-62) = 2^-94.
       Then we have |(t1+t3)+t2| < 2^-31+2^-62+2^-29 and the rounding error on (t1+t3)+t2 is bounded
       by ulp(2^-31+2^-62+2^-29) = 2^-92.
       The absolute error is thus bounded by:
       2^-96.054 + 2^-92.762 + 2^-189.818 + 2^-95 + 2^-93 + 2^-126 + 2^-94 + 2^-92 < 2^-90.663:
       | hh + ll - 2^(i2/2^5) * 2^(i1/2^10) * 2^(i0/2^15) | < 2^-90.663
       with |hh| < 2 and |ll| < 2^-28.678.
  */
  d_mul1 (h, l, *h, *l, hh, ll);
  /* At input, we have 0.999989 < h_in + l_in < 1.000011 and 1 <= hh + ll < 1.999958,
     thus at output 0.999989 < h + l < 1.999980. The different errors are:
   * that on h_in + l_in, bounded by 2^-78.947 (relative), which gives an absolute error
     bounded by 2^-78.947*1.000011*1.999980 < 2^-77.946
   * that on hh + ll, bounded by 2^-90.663 (absolute), which gives an absolute error
     bounded by 1.000011*2^-90.663 < 2^-90.662
   * the d_mul1() call decomposes into:
     - hi = ahh * bhh [exact]
     - t1 = ahh * (bhl + bl) where ahh is the upper part of h_in, bhl is the lower part of hh,
       and bl = ll. We have |bhl| < ulp(C) = 2^-31 where C is the magic constant in d_mul1(),
       and |bl| < 2^-28.678 thus |bhl + bl| < 2^-28.414 and the rounding error on bhl+bl is
       bounded by ulp(2^-28.414) = 2^-92. This error is multiplied by |ahh| < 2 thus contributes
       to at most 2^-91. Now |t1| < 2*2^-28.414 thus the rounding error on t1 is bounded by
       ulp(2^-27.414) = 2^-91. The total rounding error on t1 is thus at most 2^-91+2^-91=2^-90.
     - t2 = (ahl + al) * bhh where ahl is the lower part of h_in, al = l, and bhh is the upper
       part of hh. We have |ahl| < ulp(C) = 2^-31 and |al| < 2^-63, thus the rounding error on
       ahl+al is bounded by ulp(2^-31+2^-63) = 2^-94. This error is multiplied by |bhh| < 2 thus
       contributes to at most 2^-93. Now |t2| < (2^-31+2^-63)*2 = 2^-30+2^-62 thus the rounding
       error on t2 is bounded by ulp(2^-30+2^-62) = 2^-93. The total rounding error on t2 is thus
       bounded by 2^-93+2^-93 = 2^-92.
     - t3 = (ahl + al) * (bhl + bl). We have |ahl+al| < 2^-31+2^-63 and |bhl+bl| < 2^-31+2^-28.678,
       thus the rounding error on ahl + al, which is bounded by 2^-94, is multiplied by at most
       2^-31+2^-28.678 and contributes to at most 2^-122.414. The rounding error on bhl + bl, which
       is bounded by 2^-92, is multiplied by at most 2^-31+2^-63 and contributes to at most
       2^-122.999. Now |t3| < (2^-31+2^-63)*(2^-31+2^-28.678) < 2^-59.414 and the rounding error on
       t3 is bounded by ulp(2^-59.414) = 2^-123. The total rounding error on t3 is thus bounded by
       2^-122.414 + 2^-122.999 + 2^-123 < 2^-121.191.
     - lo = t1 + (t2 + t3). We have |t2| < 2^-30+2^-62 < 2^-29.999 and |t3| < 2^-59.414 thus
       |t2+t3| < 2^-29.998 and the rounding error on t2+t3 is bounded by ulp(2^-29.998) = 2^-93.
       Now |t1 + (t2+t3)| < 2^-27.414 + 2^-29.998 < 2^-27.191 thus the rounding error on
       t1 + (t2+t3) is bounded by ulp(2^-27.191) = 2^-91.  The total rounding error on lo is thus
       bounded by 2^-93 + 2^-91 < 2^-90.678.
     The absolute error on h + l is thus bounded by:
     2^-77.946 + 2^-90.662 + 2^-90 + 2^-92 + 2^-121.191 + 2^-90.678 < 2^-77.945.
     Since |h + l| > 0.999989 this yields a relative error < 2^-77.945/0.999989 < 2^-77.944.
   */
  if (__builtin_expect (e >= -16355, 1))
  {
    /* Multiply h, l by 2^e. Since e >= -16355, we have 2^x>=0.99998*2^-16355
       thus if l*2^e is in the subnormal range, we have an additional absolute
       error of at most 2^-16445, which corresponds to an additional relative
       error < 2^-16445/(0.99998*2^-16355) < 2^-89.999. This gives a final
       bound of (1 + 2^-77.944) * (1 + 2^-89.999) - 1 < 2^-77.943.
       No overflow is possible here since x < 16384. */
    // since |h| > 0.5, |h*2^e| > 2^-16356 and is exactly representable
    v.f = *h;
    v.e += e;
    *h = v.f;
    b80u80_t w = {.f = *l};
    if (__builtin_expect ((w.e & 0x7fff) + e > 0, 1))
      {
        w.e += e;
        *l = w.f;
      }
    else
      *l = __builtin_ldexpl (*l, e);
  }
  else
  {
    v.e = 32767;
    v.m = 0xc000000000000000ul;
    *h = *l = v.f; // +qnan
  }
}

static void
accurate_path (long double *h, long double *l, long double x)
{
#define EXCEPTIONS 59
static const long double exceptions[EXCEPTIONS][3] = {
    {-0xb.8aa3b295c17f0bcp-68L, 0x1.fffffffffffffffep-1L, 0x1.fffffffffffffffep-66L},
    {0xb.8aa3b295c17f0bcp-67L, 0x1.0000000000000002p+0L, -0x1.fffffffffffffffep-65L},
    {0xa.194f3c43094f2a2p-64L, 0x1.0000000000000006p+0L, 0x1.fffffffffffffffep-65L},
    {0xc.434dedbf1d96fc1p-63L, 0x1.0000000000000012p+0L, -0x1.fffffffffffffffep-65L},
    {0xb.6fc4ed79fcd7255p-53L, 0x1.0000000000003f6ap+0L, 0x1.fffffffffffffffep-65L},
    {0xf.49f104ab3cc2d94p-52L, 0x1.000000000000a98ep+0L, 0x1.fffffffffffffffep-65L},
    {0x9.f1ecf60af3e5853p-47L, 0x1.00000000000dc966p+0L, 0x1.fffffffffffffffep-65L},
    {0xc.3dc8cf1463af62fp-47L, 0x1.000000000010f85ap+0L, -0x1.fffffffffffffffep-65L},
    {0x9.ad1f062a8ab29ffp-40L, 0x1.0000000006b50272p+0L, 0x1.fffffffffffffffep-65L},
    {0xd.abfd779809f67b6p-38L, 0x1.0000000025e8087ap+0L, -0x1.fffffffffffffffep-65L},
    {0xc.762d7684ae1beeap-37L, 0x1.00000000451a19cep+0L, 0x1.fffffffffffffffep-65L},
    {0xe.0c9e1609da847dbp-37L, 0x1.000000004de7e1e2p+0L, 0x1.fffffffffffffffep-65L},
    {0x9.aab514ef3077eddp-36L, 0x1.000000006b3561fep+0L, -0x1.fffffffffffffffep-65L},
    {0xd.f39d71dc272a58p-29L, 0x1.0000004d5d3d3d86p+0L, -0x1.fffffffffffffffep-65L},
    {0xa.824ad65265e94b6p-25L, 0x1.000003a4626653aap+0L, 0x1.fffffffffffffffep-65L},
    {0xd.0527fc86dd2ec59p-25L, 0x1.000004832f1eead2p+0L, -0x1.fffffffffffffffep-65L},
    {0xd.ca1bcc03e818338p-25L, 0x1.000004c7714ce422p+0L, 0x1.fffffffffffffffep-65L},
    {0xc.5f396165dfc60bap-11L, 0x1.0112fe9112c95b06p+0L, 0x1.fffffffffffffffep-65L},
    {0x1.1cac23cf32997fa6p-6L, 0x1.031a0d2f944dc4d8p+0L, 0x1.fc33e05ac1b1158ap-129L},
    {0x1.248230c2bb787ce4p-16L, 0x1.0000cac0b15d6024p+0L, -0x1.ab58fc5c42eab87p-130L},
    {0x1.2574cfe96b07e51ep-15L, 0x1.000196d25dbbb85p+0L, -0x1.650ba11717cb4bbcp-130L},
    {0x1.270a4a527eb90b6cp-7L, 0x1.019a4aa31b259dccp+0L, -0x1.7e68a9c64a6a7efp-131L},
    {0x1.35e0b2e14748db7cp-7L, 0x1.01aefe25aea5272ap+0L, -0x1.80c0b33e4cf8aac2p-127L},
    {0x1.3ac9a43d4e7d192ep-5L, 0x1.06e901f58091b67ap+0L, 0x1.120ee5fe92e5b42cp-129L},
    {0x1.3f02d33da85d3b6ep-2L, 0x1.3db3eddfcd080064p+0L, 0x1.7075b144578cbff8p-129L},
    {0x1.491705f0ae9f98bep-4L, 0x1.0ea943b7cdc4830cp+0L, -0x1.97b4ec60a25776eep-126L},
    {0x1.4df4919b6022268cp-6L, 0x1.03a47e1e06af54d4p+0L, -0x1.08060332aa1ef138p-128L},
    {0x1.50919d96b5fae21p-5L, 0x1.0765299e343f756ep+0L, 0x1.c4f0626b24f2151cp-127L},
    {0x1.5178a614b366f2fap-5L, 0x1.076a4fcbe306eadp+0L, 0x1.dc18dc836e58cc56p-125L},
    {0x1.529f4845f565b744p-2L, 0x1.41f2cb598284c76ap+0L, 0x1.d2f63b235d1b5822p-128L},
    {0x1.58b0bc0151b40e26p+0L, 0x1.457c21a3a033a3ecp+1L, -0x1.56dfc93184a53a02p-126L},
    {0x1.5afc7d79dedd2a4cp-6L, 0x1.03c92571dc388a4cp+0L, 0x1.78fb4b5ddf1a16ccp-129L},
    {0x1.5ead8ebb36c52e3p-16L, 0x1.0000f312bd341228p+0L, 0x1.ef4c0926ab586534p-132L},
    {0x1.5f5b152690eba5dap-13L, 0x1.00079c717ef7efcp+0L, 0x1.313adf5b534e0502p-127L},
    {0x1.62c2f00546d03898p-2L, 0x1.457c21a3a033a3ecp+0L, -0x1.56dfc93184a53a02p-127L},
    {0x1.658382b8511ee5ccp-10L, 0x1.003dfb508259ecacp+0L, 0x1.aff6ac6986857a6cp-126L},
    {0x1.6ec1e220c34be404p-1L, 0x1.a49af00837c3b46ap+0L, 0x1.55129bf7e816581ap-128L},
    {0x1.6f9ce5a8b3243262p-7L, 0x1.01ff9b337f526032p+0L, 0x1.25f7555adb61477cp-128L},
    {0x1.70fd6310d1b4994cp-6L, 0x1.0407157c0ce85144p+0L, 0x1.0e68d791be9eb2fcp-133L},
    {0x1.7d098c9ba167b4bap+0L, 0x1.6725658526f34c7ap+1L, -0x1.977481b2530f44f6p-127L},
    {0x1.8dae021561102834p-2L, 0x1.4f145246ca66c496p+0L, 0x1.38c74600bb4d06a4p-125L},
    {0x1.a4ed7fbb4a9fb356p-4L, 0x1.12e68526b08d8282p+0L, -0x1.dbb94f6d0a942a3ap-127L},
    {0x1.aaded45884e59364p-12L, 0x1.00127ed001fc8accp+0L, -0x1.0ac20ca1ef316aeep-128L},
    {0x1.ad988d3081bcbb9cp-4L, 0x1.134dd395bd76f908p+0L, 0x1.dc94128e60787ebp-127L},
    {0x1.ae30b1e652dca39ap-12L, 0x1.0012a3a3fccb6446p+0L, 0x1.6106632122af6d9cp-129L},
    {0x1.b3aa5032fa7f12c8p-1L, 0x1.cdbb2250ecf28d18p+0L, 0x1.51f7c471f44bbd42p-125L},
    {0x1.b760f11061a5f202p+0L, 0x1.a49af00837c3b46ap+1L, 0x1.55129bf7e816581ap-127L},
    {0x1.c400323ab65060d8p-4L, 0x1.14598c62848ce032p+0L, 0x1.a574d511f0618ab2p-127L},
    {0x1.cf8852012559841ep-2L, 0x1.5e5a8e406ecbb63ap+0L, 0x1.ab1104fa34c02b38p-131L},
    {0x1.d00a4c793a1d6d4ep-16L, 0x1.000141a6b8f91d42p+0L, -0x1.b86975165f93cd9p-128L},
    {0x1.d2eb2bfd12d6f486p-4L, 0x1.150c5eb3832acc14p+0L, 0x1.2883e8680287fe9ap-128L},
    {0x1.d9d528197d3f8964p+0L, 0x1.cdbb2250ecf28d18p+1L, 0x1.51f7c471f44bbd42p-124L},
    {0x1.db4b22a09e022f6p-13L, 0x1.000a4bcb36ef561p+0L, -0x1.56ab41256e8ece16p-130L},
    {0x1.e2dda3cd8c341298p-11L, 0x1.0029d9b9a11881b8p+0L, -0x1.1422c5751fe6962cp-128L},
    {0x1.e5b7eae7259fcb4cp-5L, 0x1.0abd81e709e4f1a4p+0L, 0x1.6109741735fe354ap-127L},
    {0x1.eaab0d7de0384c5ap-3L, 0x1.2e3f3978515cbfap+0L, 0x1.57a35d3d4f378412p-126L},
    {0x1.eb990e74b7582b7p-5L, 0x1.0adf7c7d0f3e7b3p+0L, 0x1.7449760cad2f03d4p-125L},
    {0x1.ecea940cbe9fc4b2p+1L, 0x1.cdbb2250ecf28d18p+3L, 0x1.51f7c471f44bbd42p-122L},
    {0x1.f426326e859ed2e8p-2L, 0x1.6725658526f34c7ap+0L, -0x1.977481b2530f44f6p-128L},
  };
  for (int i = 0; i < EXCEPTIONS; i++)
  {
    if (x == exceptions[i][0])
    {
      *h = exceptions[i][1];
      *l = exceptions[i][2];
      return;
    }
  }
  
  int32_t k = __builtin_roundl (0x1p15L * x); // -16445*2^15 <= k <= 16383*2^15
  // if (x == TRACE) printf ("k=%d\n", k);
  long double r = x - (long double) k * 0x1p-15L;
  // if (x == TRACE) printf ("r=%La\n", r);
  int32_t i = (k + 538869760) & 32767;
  // if (x == TRACE) printf ("i=%d\n", i);
  int32_t e = (k - i) >> 15;
  // if (x == TRACE) printf ("e=%d\n", e);
  int32_t i0 = i & 0x1f, i1 = (i >> 5) & 0x1f, i2 = i >> 10;
  // if (x == TRACE) printf ("i2=%d i1=%d i0=%d\n", i2, i1, i0);
  Pacc (h, l, r);
  // if (x == TRACE) printf ("P: h=%La l=%La\n", *h, *l);
  long double hh, ll;
  d_mul (&hh, &ll, T2[i2][0], T2[i2][1], T1[i1][0], T1[i1][1]);
  d_mul (&hh, &ll, hh, ll, T0[i0][0], T0[i0][1]);
  d_mul (h, l, *h, *l, hh, ll);
  // normalize h+l
  fast_two_sum (h, l, *h, *l);
  // if (x == TRACE) printf ("x=%La h=%La l=%La e=%d\n", x, *h, *l, e);
  if (e >= -16381)
  {
    /* Since |h| > 0.5, ulp(h) >= 2^-64, thus ulp(h)*2^e >= 2^-16445 which is the smallest
       subnormal, thus 2^e*h is exact. */
    *h = __builtin_ldexpl (*h, e);
    *l = __builtin_ldexpl (*l, e);
  }
  else // near subnormal range
  {
    hh = *h;
    *h = __builtin_ldexpl (*h, e); // might not equal 2^e*h
    // if (x == TRACE) printf ("h=%La\n", *h);
    hh = hh - __builtin_ldexpl (*h, -e); // remaining (truncated) part
    // if (x == TRACE) printf ("hh=%La\n", hh);
    hh += *l;
    // if (x == TRACE) printf ("hh=%La\n", hh);
    *l = __builtin_ldexpl (hh, e);
    // if (x == TRACE) printf ("l=%La\n", *l);
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
  if (__builtin_expect (e <= 16319, 0)) // |x| < 2^-63
  {
    if (0 <= x && x <= 0x1.71547652b82fe176p-64L)
      return __builtin_fmal (x, x, 0x1p0L);
    if (-0x1.71547652b82fe176p-65L <= x && x < 0)
      return __builtin_fmal (x, -x, 0x1p0L);
  }

  // now -16446 < x < -0x1.71547652b82fe176p-65 or 0x1.71547652b82fe176p-64 < x < 16384

  long double h, l;
  fast_path (&h, &l, x);
  static const long double err = 0x1.0bp-78; // 2^-77.943 < err
  //if (x == TRACE) printf ("h=%La l=%La\n", h, l);
  long double left = h +  (l - h * err);
  long double right = h + (l + h * err);
  //if (x == TRACE) printf ("left=%La right=%La\n", left, right);
  if (__builtin_expect (left == right, 1))
    return left;

  //if (x == TRACE) printf ("fast path failed\n");

  accurate_path (&h, &l, x);
  return h + l;
}
