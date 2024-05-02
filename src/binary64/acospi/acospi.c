/* Correctly-rounded acospi of binary64 value.

Copyright (c) 2021-2023 Paul Zimmermann and Tom Hubrecht

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

/* This implementation is derived from the CORE-MATH acos.c code
   (revision b9d8cae), using the formula acospi (x) = acos(x) / pi.
   The parts specific to acospi are marked either with "acospi_specific",
   or enclosed between "acospi_begin" and "acospi_end", to make it easier
   to maintain this code if acos is fixed or improved.
*/

#include <stdint.h>
#include <math.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

/****************** code copied from pow.[ch] ********************************/

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma (a, b, -*hi);
}

// Returns (ah + al) * (bh + bl) - (al * bl)
static inline void d_mul(double *hi, double *lo, double ah, double al,
                         double bh, double bl) {
  a_mul (hi, lo, ah, bh);
  *lo = __builtin_fma (ah, bl, *lo);
  *lo = __builtin_fma (al, bh, *lo);
}

/****************** end of code copied from pow.[ch] *************************/

// acospi_begin
/* ONE_OVER_PIH + ONE_OVER_PIL is a double-double approximation of 1/pi,
   with | ONE_OVER_PIH + ONE_OVER_PIL - 1/pi | < 2^-109.523 */
#define ONE_OVER_PIH 0x1.45f306dc9c883p-2
#define ONE_OVER_PIL -0x1.6b01ec5417056p-56
// acospi_end

/* h + l <- a * b */
static void
dekker (double *h, double *l, double a, double b)
{
  *h = a * b;
#ifdef __FP_FAST_FMA
  *l = __builtin_fma (a, b, -*h);
#else /* use Dekker's algorithm */
#define MAGIC 0x8000001
  double ah = a * MAGIC, bh = b * MAGIC;
  ah = (a - ah) + ah;
  bh = (b - bh) + bh;
  double al = a - ah;
  double bl = b - bh;
  *l = (((ah * bh - *h) + ah * bl) + al * bh) + al * bl;
#endif
}

/* return h=sqrt(x) and l such that h+l approximates sqrt(x) to about
   106 bits of accuracy */
static double
sqrt_dbl_dbl (double x, double *l)
{
  double h = __builtin_sqrt (x);
#ifdef __FP_FAST_FMA
  double e = -__builtin_fma (h, h, -x);
#else
  double u, v;
  dekker (&u, &v, h, h);
  double e = (x - u) - v ;
#endif
  /* e = x - h^2 */
  *l = e / (h + h);
  return h;
}

/* s + t <- a + b, assuming |a| >= |b| */
static void
fast_two_sum (double *s, double *t, double a, double b)
{
  *s = a + b;
  double e = *s - a;
  *t = b - e;
}

/* h_out + l_out = (h_in + l_in) * y */
static void
mul2_1 (double *h_out, double *l_out, double h_in, double l_in, double y)
{
  dekker (h_out, l_out, h_in, y);
  *l_out += l_in * y;
}

/* The following table was generated with
   Gen_P_aux(6,256,threshold=192,verbose=true,extra0=12,extra1=3,k1=192)
   for the first part (0 <= i < 192), and for the second part (192 <= i < 256)
   with
   Gen_P_aux(5,256,threshold=192,verbose=true,extra0=12,extra1=3,k0=192)
   For 0 <= i < 192, each entry contains 8 coefficients corresponding to
   a degree-6 polynomial and its evaluation point:
   * T[i][0] is p[0] (the upper 12 bits of the low part are zero)
   * T[i][1] is p[1] (the upper 3 bits of the low part are zero)
   * T[i][2] is p[2]
   * T[i][3] is p[3]
   * T[i][4] is p[4]
   * T[i][5] is p[5]
   * T[i][6] is p[6]
   * T[i][7] is the evaluation point (xmid).
   For 192 <= i < 256, each entry contains 7 coefficients corresponding to
   a degree-5 polynomial and its evaluation point:
   * T[i][0] is p[0]
   * T[i][1] is p[1]
   * T[i][2] is p[2]
   * T[i][3] is p[3]
   * T[i][4] is p[4]
   * T[i][5] is p[5]
   * T[i][6] is the evaluation point (xmid).
*/
static const double T[256][8] = {
  { 0x1.921fb54442d18p+0, /*0x4.6989e4b05fa3p-56,*/ -0x1p+0, -0xe.3605b41a5d07p-52, -0x2.aaaaaaaa26314p-4, -0x1.d02e84f046fbp-28, -0x1.3330337ffdfc1p-4, -0x2.5ca13405c66b6p-12, 0x0p+0 }, /* i=0 73.489 */
  { 0x1.909fb4b44253fp+0, -0x1.00012001e6037p+0, -0xc.00288071475c8p-12, -0x2.aac5ab293276cp-4, -0x9.005461cb2ecep-12, -0x1.3355447aa843dp-4, -0x7.809a8297c2bbcp-12, 0x1.7fffffffebc6fp-8 }, /* i=1 74.13 */
  { 0x1.8f9fb2a990bfbp+0, -0x1.0003200ea64c8p+0, -0x1.400bb85b9a23fp-8, -0x2.aaf5ae7b579bep-4, -0xf.0186b75746918p-12, -0x1.33914cbc551f1p-4, -0xc.8295e5760bb98p-12, 0x2.80000000140e4p-8 }, /* i=2 74.124 */
  { 0x1.8e9fae1ec5e86p+0, -0x1.00062038483f3p+0, -0x1.c02029ec97416p-8, -0x2.ab3db95355f9p-4, -0x1.504305dafd9f7p-8, -0x1.33eb677d14e3dp-4, -0x1.186edf2c1538bp-8, 0x3.800000002e7b8p-8 }, /* i=3 74.117 */
  { 0x1.8d9fa613b89eap+0, -0x1.000a2099d022ap+0, -0x2.40445ec276a82p-8, -0x2.ab9dd2bb5d952p-4, -0x1.b08e7d9b5ed11p-8, -0x1.3463a5fffc8c4p-4, -0x1.68eb18bc28196p-8, 0x4.7fffffffc6fbcp-8 }, /* i=4 74.108 */
  { 0x1.8c9f99881e5bcp+0, -0x1.000f215747cefp+0, -0x2.c07cda70f2514p-8, -0x2.ac16041872d36p-4, -0x2.1104415368d9p-8, -0x1.34fa1f513523fp-4, -0x1.b9aa6a457192fp-8, 0x5.8000000039eep-8 }, /* i=5 74.099 */
  { 0x1.8b9f877b8725bp+0, -0x1.0015229dc21c5p+0, -0x3.40ce2285cdfe2p-8, -0x2.aca6592c394f6p-4, -0x2.71adc7b80261cp-8, -0x1.35aef04dfe9d2p-4, -0x2.0ac0bc5ec1148p-8, 0x6.7fffffffc9304p-8 }, /* i=6 74.09 */
  { 0x1.8a9f6eed4e15fp+0, -0x1.001c24a35f6c6p+0, -0x3.c13cbefda1e8ap-8, -0x2.ad4ee01738444p-4, -0x2.d29491715f22cp-8, -0x1.36823bad3f08bp-4, -0x2.5c3a32c94c508p-8, 0x7.8000000009a98p-8 }, /* i=7 74.078 */
  { 0x1.899f4edc9635ap+0, -0x1.002427a752d4cp+0, -0x4.41cd3ab2eccd8p-8, -0x2.ae0fa95b8dcp-4, -0x3.33c22ae8a1b92p-8, -0x1.37742a0b7a1e4p-4, -0x2.ae29dc08ce6d2p-8, 0x8.7ffffffff7918p-8 }, /* i=8 74.064 */
  { 0x1.889f26483e1ecp+0, -0x1.002d2bf1e8623p+0, -0x4.c28423d24ffe8p-8, -0x2.aee8c7e01c9eep-4, -0x3.95402e27b217cp-8, -0x1.3884e9f63c6e9p-4, -0x3.009d1698267ecp-8, 0x9.7fffffff3f388p-8 }, /* i=9 74.049 */
  { 0x1.879ef42ed5a5bp+0, -0x1.003731d48c3fap+0, -0x5.43660c4e68278p-8, -0x2.afda50f431414p-4, -0x3.f71844b60aeb4p-8, -0x1.39b4affcb63bp-4, -0x3.539f5bba69818p-8, 0xa.7fffffff76be8p-8 }, /* i=10 74.036 */
  { 0x1.869eb78e991c6p+0, -0x1.004239a9d28aep+0, -0x5.c4778a518ae48p-8, -0x2.b0e45c5398ecep-4, -0x4.5954297214a34p-8, -0x1.3b03b6be3ce92p-4, -0x3.a74d8446ef28ap-8, 0xb.7fffffffb7a98p-8 }, /* i=11 74.021 */
  { 0x1.859e6f6564ab3p+0, -0x1.004e43d5805a6p+0, -0x6.45bd38b44adbcp-8, -0x2.b207042b3956cp-4, -0x4.bbfdaa8ac5514p-8, -0x1.3c723f003ff56p-4, -0x3.fbae296f78cd8p-8, 0xc.800000001912p-8 }, /* i=12 74.004 */
  { 0x1.849e1ab0aca76p+0, -0x1.005b50c4957e2p+0, -0x6.c73bb7725c3f8p-8, -0x2.b342651e1e994p-4, -0x5.1f1eab5ff321cp-8, -0x1.3e008fbf7c925p-4, -0x4.50d3670de1bep-8, 0xd.8000000023c4p-8 }, /* i=13 73.987 */
  { 0x1.839db86d7338ep+0, -0x1.006960ed574a5p+0, -0x7.48f7ac21cd944p-8, -0x2.b4969e4b085ecp-4, -0x5.82c1267cb2a44p-8, -0x1.3faef648f515ap-4, -0x4.a6cf5d1c9e298p-8, 0xe.7fffffffec158p-8 }, /* i=14 73.969 */
  { 0x1.829d47983f084p+0, -0x1.007874cf5c462p+0, -0x7.caf5c26acc3e4p-8, -0x2.b603d152745a8p-4, -0x5.e6ef2f93f1878p-8, -0x1.417dc652943adp-4, -0x4.fdb49ecb2ee9cp-8, 0xf.80000000188fp-8 }, /* i=15 73.947 */
  { 0x1.819cc72d142dbp+0, -0x1.00888cf398a56p+0, -0x8.4d3aac7f6b99p-8, -0x2.b78a225d25584p-4, -0x6.4bb2f588a1cdp-8, -0x1.436d5a16421b5p-4, -0x5.558bce13bdd5p-8, 0x1.07fffffff9c8fp-4 }, /* i=16 73.926 */
  { 0x1.809c3627667cep+0, -0x1.0099a9ec6c1b3p+0, -0x8.cfcb2398098bp-8, -0x2.b929b82336418p-4, -0x6.b116c470aedb4p-8, -0x1.457e126dfa8e5p-4, -0x5.ae725e241fdap-8, 0x1.18000000005e6p-4 }, /* i=17 73.901 */
  { 0x1.7f9b938215cf2p+0, -0x1.00abcc55b001dp+0, -0x9.52abe86bf7968p-8, -0x2.bae2bbf3a27ep-4, -0x7.172507b730604p-8, -0x1.47b056f437d6ep-4, -0x6.08797fdf25adp-8, 0x1.27fffffffc329p-4 }, /* i=18 73.881 */
  { 0x1.7e9ade375e48ep+0, -0x1.00bef4d4c742ap+0, -0x9.d5e1c3b1b0b2p-8, -0x2.bcb559bc684d8p-4, -0x7.7de84c3c1cbd8p-8, -0x1.4a049625c398ep-4, -0x6.63b04475fb8dcp-8, 0x1.3800000005926p-4 }, /* i=19 73.859 */
  { 0x1.7d9a1540d593ap+0, -0x1.00d32418ae534p+0, -0xa.59718699fe43p-8, -0x2.bea1c0131ee6ep-4, -0x7.e56b427ac8b98p-8, -0x1.4c7b4585fc60fp-4, -0x6.c02b232597bc8p-8, 0x1.47fffffff9be3p-4 }, /* i=20 73.839 */
  { 0x1.7c99379758aa9p+0, -0x1.00e85ada0d548p+0, -0xa.dd600b54a1fe8p-8, -0x2.c0a8203e34c9cp-4, -0x8.4db8c0c794ac8p-8, -0x1.4f14e1c4ac981p-4, -0x7.1dfe2bcaa789p-8, 0x1.580000000393cp-4 }, /* i=21 73.813 */
  { 0x1.7b9844330bddep+0, -0x1.00fe99db49c9ap+0, -0xb.61b2358d6977p-8, -0x2.c2c8ae3e9b126p-4, -0x8.b6dbc5930227p-8, -0x1.51d1eee6a4e68p-4, -0x7.7d366b13065a8p-8, 0x1.68000000008bep-4 }, /* i=22 73.788 */
  { 0x1.7a973a0b48f84p+0, -0x1.0115e1e89abf8p+0, -0xb.e66cf2f442008p-8, -0x2.c503a0da260bap-4, -0x9.20df79b5bba78p-8, -0x1.54b2f873a0e64p-4, -0x7.ddefcbc5733acp-8, 0x1.7800000001812p-4 }, /* i=23 73.758 */
  { 0x1.79961816995f6p+0, -0x1.012e33d81d053p+0, -0xc.6b953bc11a3a8p-8, -0x2.c75931a66c2d6p-4, -0x9.8bcf32dc703b8p-8, -0x1.57b891a2a8474p-4, -0x8.4032bdd1257dp-8, 0x1.880000000b248p-4 }, /* i=24 73.731 */
  { 0x1.7894dd4aad755p+0, -0x1.01479089e89e7p+0, -0xc.f130133b373a8p-8, -0x2.c9c99d1443048p-4, -0x9.f7b675dfcd44p-8, -0x1.5ae3558a8e26bp-4, -0x8.a42ad40c037ap-8, 0x1.9800000000ff4p-4 }, /* i=25 73.705 */
  { 0x1.7793889c4dcbdp+0, -0x1.0161f8e827e39p+0, -0xd.77428845eb048p-8, -0x2.cc55227be899p-4, -0xa.64a0f962883a8p-8, -0x1.5e33e75926448p-4, -0x9.09dcfb7973f4p-8, 0x1.a7fffffff17b3p-4 }, /* i=26 73.674 */
  { 0x1.769218ff52aeap+0, -0x1.017d6de72f22bp+0, -0xd.fdd1b5ec4de8p-8, -0x2.cefc0429c02fp-4, -0xa.d29aa84b1257p-8, -0x1.61aaf285d96ecp-4, -0x9.715ac9d491178p-8, 0x1.b7fffffffb717p-4 }, /* i=27 73.643 */
  { 0x1.75908d669e9ecp+0, -0x1.0199f08594facp+0, -0xe.84e2c3edcf2cp-8, -0x2.d1be876bab71ap-4, -0xb.41afa45e93158p-8, -0x1.65492b0f68076p-4, -0x9.dacb9f52cbd58p-8, 0x1.c800000003818p-4 }, /* i=28 73.613 */
  { 0x1.748ee4c40f636p+0, -0x1.01b781cc4cd07p+0, -0xf.0c7ae75234ebp-8, -0x2.d49cf49f2623cp-4, -0xb.b1ec490698e68p-8, -0x1.690f4db6cd4p-4, -0xa.4643cb605d618p-8, 0x1.d7fffffffdaabp-4 }, /* i=29 73.586 */
  { 0x1.738d1e087244bp+0, -0x1.01d622cec2209p+0, -0xf.949f62fea16e8p-8, -0x2.d79797400a2c6p-4, -0xc.235d2e10a8bap-8, -0x1.6cfe2040ea4a9p-4, -0xa.b3d85c425f238p-8, 0x1.e800000002cdp-4 }, /* i=30 73.547 */
  { 0x1.728b38237fc97p+0, -0x1.01f5d4aaf412dp+0, -0x1.01d55884965c6p-4, -0x2.daaebdf7f10cp-4, -0xc.960f2a874efcp-8, -0x1.711671b94ba69p-4, -0xb.239e9f0c1c8b8p-8, 0x1.f7ffffffed1edp-4 }, /* i=31 73.511 */
  { 0x1.71893203c6f8fp+0, -0x1.0216988994331p+0, -0x1.0a6a2b7997277p-4, -0x2.dde2baae7fb84p-4, -0xd.0a0f57a961478p-8, -0x1.75591abb938bbp-4, -0xb.95a8b052283d8p-8, 0x2.07fffffff8b82p-4 }, /* i=32 73.485 */
  { 0x1.70870a96adc0dp+0, -0x1.02386f9e23cdp+0, -0x1.1308c60fda193p-4, -0x2.e133e29a31444p-4, -0xd.7f6b13df54cbp-8, -0x1.79c6fdbcfdd48p-4, -0xc.0a285b98d4f48p-8, 0x2.1800000012728p-4 }, /* i=33 73.451 */
  { 0x1.6f84c0c86028fp+0, -0x1.025b5b271497p+0, -0x1.1bb1804d0488ep-4, -0x2.e4a28e521a2f2p-4, -0xd.f63005f9011ep-8, -0x1.7e61075b239fp-4, -0xc.812ccfb30e728p-8, 0x2.2800000035fc4p-4 }, /* i=34 73.411 */
  { 0x1.6e825383cd1afp+0, -0x1.027f5c6de8d69p+0, -0x1.2464b3456dc13p-4, -0x2.e82f19e049f4ep-4, -0xe.6e6c2055b0758p-8, -0x1.83282eb1ebe64p-4, -0xc.fac5acd3a4c8p-8, 0x2.37fffffff27dep-4 }, /* i=35 73.377 */
  { 0x1.6d7fc1b280c44p+0, -0x1.02a474c7599bp+0, -0x1.2d22b9277965ep-4, -0x2.ebd9e4d577614p-4, -0xe.e82da4310dacp-8, -0x1.881d75af61cedp-4, -0xd.773458d78cfc8p-8, 0x2.480000000e48p-4 }, /* i=36 73.334 */
  { 0x1.6c7d0a3cc04edp+0, -0x1.02caa59375263p+0, -0x1.35ebed44f8951p-4, -0x2.efa3525c9ebc2p-4, -0xf.6383253880d1p-8, -0x1.8d41e972e5acp-4, -0xd.f6634e8b35108p-8, 0x2.580000001cc4cp-4 }, /* i=37 73.302 */
  { 0x1.6b7a2c0965d62p+0, -0x1.02f1f03dc7808p+0, -0x1.3ec0ac1eecacp-4, -0x2.f38bc9505e116p-4, -0xf.e07b8cdad4bb8p-8, -0x1.9296a2aab38d5p-4, -0xe.789cd4b1e4b7p-8, 0x2.6800000004e6p-4 }, /* i=38 73.264 */
  { 0x1.6a7725fdd5bc9p+0, -0x1.031a563d81032p+0, -0x1.47a15370afc09p-4, -0x2.f793b450e2d7ep-4, -0x1.05f261e2d7b9fp-4, -0x1.981cc5fa65a18p-4, -0xe.fdfddd3bfefbp-8, 0x2.77fffffff2bbep-4 }, /* i=39 73.227 */
  { 0x1.6973f6fdf680fp+0, -0x1.0343d9159dd8dp+0, -0x1.508e423b48b34p-4, -0x2.fbbb81dabb2f4p-4, -0x1.0df9279aafcc4p-4, -0x1.9dd58465923cp-4, -0xf.8696701d4ec38p-8, 0x2.880000001f244p-4 }, /* i=40 73.184 */
  { 0x1.68709dec32aa8p+0, -0x1.036e7a550d25cp+0, -0x1.5987d8d0a96fep-4, -0x3.0003a45e70268p-4, -0x1.161d0a111844cp-4, -0x1.a3c21bbece4cdp-4, -0x1.012818df20cc4p-4, 0x2.97fffffff5e62p-4 }, /* i=41 73.151 */
  { 0x1.676d19a94e1bcp+0, -0x1.039a3b96e0e74p+0, -0x1.628e78e0bf137p-4, -0x3.046c9259d9ff4p-4, -0x1.1e5f0fb7dd5ddp-4, -0x1.a9e3d71bc4af6p-4, -0x1.0a1f6dd2593c5p-4, 0x2.a7fffffff9bb2p-4 }, /* i=42 73.105 */
  { 0x1.666969147661cp+0, -0x1.03c71e8275c5fp+0, -0x1.6ba28584c398ep-4, -0x3.08f6c67188262p-4, -0x1.26c0459a58df2p-4, -0x1.b03c0f5207f82p-4, -0x1.135226a841bfcp-4, 0x2.b800000001b18p-4 }, /* i=43 73.071 */
  { 0x1.65658b0b297aep+0, -0x1.03f524cba3263p+0, -0x1.74c4634c4bf68p-4, -0x3.0da2bf8be7a8p-4, -0x1.2f41bfa442c8ep-4, -0x1.b6cc2b758d196p-4, -0x1.1cc12f7dd6dbcp-4, 0x2.c800000003fd4p-4 }, /* i=44 73.023 */
  { 0x1.64617e692a3p+0, -0x1.04245032ea9fp+0, -0x1.7df4784a2f89bp-4, -0x3.127100ed4e5f4p-4, -0x1.37e498e85441ep-4, -0x1.bd95a160f6a9dp-4, -0x1.26704ce77d56bp-4, 0x2.d800000007654p-4 }, /* i=45 72.984 */
  { 0x1.635d420875812p+0, -0x1.0454a285a8d14p+0, -0x1.87332c21b99f1p-4, -0x3.1762125526c3ap-4, -0x1.40a9f3eceb242p-4, -0x1.c499f6433ff6p-4, -0x1.30627a5c5d56cp-4, 0x2.e800000000d74p-4 }, /* i=46 72.936 */
  { 0x1.6258d4c132bacp+0, -0x1.04861d9e48dbdp+0, -0x1.9080e81462ee2p-4, -0x3.1c76801c66a64p-4, -0x1.4992fafb1c3e8p-4, -0x1.cbdabf336ef02p-4, -0x1.3a97b45c4c021p-4, 0x2.f80000000207p-4 }, /* i=47 72.903 */
  { 0x1.61543569aaa93p+0, -0x1.04b8c364782cep+0, -0x1.99de170fa5909p-4, -0x3.21aedb55207dap-4, -0x1.52a0e06c1e2d6p-4, -0x1.d359a1cdc2e3ap-4, -0x1.4515a21e5983cp-4, 0x3.07fffffff4dbep-4 }, /* i=48 72.859 */
  { 0x1.604f62d634e06p+0, -0x1.04ec95cd5deeep+0, -0x1.a34b25bb8a91fp-4, -0x3.270bb9eb8f572p-4, -0x1.5bd4df01d79e5p-4, -0x1.db1854d5fd617p-4, -0x1.4fdd284dce08ap-4, 0x3.17ffffffef9f6p-4 }, /* i=49 72.808 */
  { 0x1.5f4a5bd92e0a3p+0, -0x1.052196dbd285ap+0, -0x1.acc8828947cc5p-4, -0x3.2c8db6c84bc88p-4, -0x1.65303a3a304dbp-4, -0x1.e318a0e3e743ep-4, -0x1.5af453ba13ccdp-4, 0x3.27fffffff7d2ap-4 }, /* i=50 72.764 */
  { 0x1.5e451f42e9f7cp+0, -0x1.0557c8a099a78p+0, -0x1.b6569dc26b1efp-4, -0x3.323571f3f8c7cp-4, -0x1.6eb43eae7c501p-4, -0x1.eb5c6115dd46ap-4, -0x1.665b1206461b2p-4, 0x3.38000000043a4p-4 }, /* i=51 72.718 */
  { 0x1.5d3fabe1a56dep+0, -0x1.058f2d3a9e659p+0, -0x1.bff5e99873377p-4, -0x3.380390bc636fp-4, -0x1.7862426dfbd5p-4, -0x1.f3e583cf00be9p-4, -0x1.721512657c562p-4, 0x3.47ffffffff82ap-4 }, /* i=52 72.682 */
  { 0x1.5c3a0081755adp+0, -0x1.05c7c6d731b24p+0, -0x1.c9a6da34f05bep-4, -0x3.3df8bddb32204p-4, -0x1.823ba560faea5p-4, -0x1.fcb60b7cd79eep-4, -0x1.7e27d5c87e887p-4, 0x3.57ffffffefa82p-4 }, /* i=53 72.629 */
  { 0x1.5b341bec33a5bp+0, -0x1.060197b24b7cap+0, -0x1.d369e5ca30f48p-4, -0x3.4415a99e33bdap-4, -0x1.8c41d1b19e43p-4, -0x2.05d00f66cf7e6p-4, -0x1.8a94b92c5f705p-4, 0x3.680000003ed6cp-4 }, /* i=54 72.587 */
  { 0x1.5a2dfce98c38p+0, -0x1.063ca216c6ap+0, -0x1.dd3f84a331e7ep-4, -0x3.4a5b0a108adfp-4, -0x1.96763c3284854p-4, -0x2.0f35bc883f4e6p-4, -0x1.975f97a8dedd6p-4, 0x3.78000000088f8p-4 }, /* i=55 72.535 */
  { 0x1.5927a23ebd7cdp+0, -0x1.0678e85eafecdp+0, -0x1.e7283136e849ep-4, -0x3.50c99b27567b8p-4, -0x1.a0da64cf7b789p-4, -0x2.18e956757d6bap-4, -0x1.a48e869e2cc49p-4, 0x3.880000000f778p-4 }, /* i=56 72.494 */
  { 0x1.58210aaebb6d3p+0, -0x1.06b66cf38158dp+0, -0x1.f12468385dc26p-4, -0x3.57621eedc7cf2p-4, -0x1.ab6fd7015f33ap-4, -0x2.22ed384ca44fap-4, -0x1.b2239b0d502c5p-4, 0x3.97ffffffaa52ap-4 }, /* i=57 72.448 */
  { 0x1.571a34f9eb46ap+0, -0x1.06f5324e76f25p+0, -0x1.fb34a8ab3949bp-4, -0x3.5e255db58f592p-4, -0x1.b6382a459425bp-4, -0x2.2d43d5b1b6e88p-4, -0x1.c023bcb6ff9fbp-4, 0x3.a7fffffffa8c2p-4 }, /* i=58 72.395 */
  { 0x1.56131fde742b5p+0, -0x1.07353af8c2f65p+0, -0x2.055973f327a5ap-4, -0x3.651426456badap-4, -0x1.c1350299056c3p-4, -0x2.37efbbd4d44a8p-4, -0x1.ce923cdd7680ap-4, 0x3.b7fffffe0c4c6p-4 }, /* i=59 72.349 */
  { 0x1.550bca175781dp+0, -0x1.0776898c0fbffp+0, -0x2.0f934defd9c5p-4, -0x3.6c2f4e11dea6p-4, -0x1.cc681101b077fp-4, -0x2.42f39290a4488p-4, -0x1.dd759c6cb4961p-4, 0x3.c7fffffff10b2p-4 }, /* i=60 72.289 */
  { 0x1.5404325db6048p+0, -0x1.07b920b27c2f3p+0, -0x2.19e2bd041e3f4p-4, -0x3.7377b1693990ep-4, -0x1.d7d3140bbb35p-4, -0x2.4e521d80d6a76p-4, -0x1.eccf208ae3d17p-4, 0x3.d7fffffffba76p-4 }, /* i=61 72.251 */
  { 0x1.52fc576782753p+0, -0x1.07fd03273bf8ap+0, -0x2.24484a370e876p-4, -0x3.7aee33b3555d8p-4, -0x1.e377d85f7b4dp-4, -0x2.5a0e3d41785p-4, -0x1.fca7b6a10ee9p-4, 0x3.e7fffffffdee6p-4 }, /* i=62 72.194 */
  { 0x1.51f437e808d47p+0, -0x1.084233b6c7882p+0, -0x2.2ec48143c2d8ep-4, -0x3.8293bfa711b9cp-4, -0x1.ef583954ae3fep-4, -0x2.662af0abb66d8p-4, -0x2.0d0210b36ce7cp-4, 0x3.f800000005bbcp-4 }, /* i=63 72.143 */
  { 0x1.50ebd28fc514dp+0, -0x1.0888b53f39c35p+0, -0x2.3957f0b08508cp-4, -0x3.8a694786de49p-4, -0x1.fb76218aae08cp-4, -0x2.72ab5629f9c36p-4, -0x2.1de2c67cff35p-4, 0x4.07ffffffd041cp-4 }, /* i=64 72.101 */
  { 0x1.4fe3260c3d6b5p+0, -0x1.08d08ab0b0357p+0, -0x2.440329e69c616p-4, -0x3.926fc55fb45bep-4, -0x2.07d38b8a38954p-4, -0x2.7f92ad1c173dap-4, -0x2.2f515f4ba4222p-4, 0x4.17fffffffdc94p-4 }, /* i=65 72.042 */
  { 0x1.4eda3108130a3p+0, -0x1.0919b70d9f765p+0, -0x2.4ec6c148a62b8p-4, -0x3.9aa83b49014a4p-4, -0x2.1472826cf4f2cp-4, -0x2.8ce4574be423cp-4, -0x2.41568bd3e2aeap-4, 0x4.27fffffff5254p-4 }, /* i=66 71.987 */
  { 0x1.4dd0f22acc672p+0, -0x1.09643d6b3d714p+0, -0x2.59a34e4c85964p-4, -0x3.a313b3a95b46cp-4, -0x2.21552293e187ep-4, -0x2.9aa3da7bddc58p-4, -0x2.53f008869a3b6p-4, 0x4.38000000044bp-4 }, /* i=67 71.939 */
  { 0x1.4cc76818d8ee5p+0, -0x1.09b020f1df239p+0, -0x2.64996b93cee0cp-4, -0x3.abb3417c46418p-4, -0x2.2e7d9a52c7d62p-4, -0x2.a8d4e1ff53b7cp-4, -0x2.672ca7adf46c4p-4, 0x4.4800000002238p-4 }, /* i=68 71.893 */
  { 0x1.4bbd917370a87p+0, -0x1.09fd64dd635e9p+0, -0x2.6fa9b7067c38p-4, -0x3.b488009c01be6p-4, -0x2.3bee2abce5f48p-4, -0x2.b77b4072d6b0ap-4, -0x2.7b0eb38661e62p-4, 0x4.5800000017b98p-4 }, /* i=69 71.839 */
  { 0x1.4ab36cd88d8ebp+0, -0x1.0a4c0c7d99782p+0, -0x2.7ad4d1ed70ae2p-4, -0x3.bd93160d91cccp-4, -0x2.49a92865b8ac6p-4, -0x2.c69af183b9498p-4, -0x2.8f9e4924ebb48p-4, 0x4.67fffffff6c3cp-4 }, /* i=70 71.784 */
  { 0x1.49a8f8e2c6ba4p+0, -0x1.0a9c1b36b409cp+0, -0x2.861b610f21aa6p-4, -0x3.c6d5b050f98a8p-4, -0x2.57b0fc35fad62p-4, -0x2.d6381bd45aa62p-4, -0x2.a4e40fc475512p-4, 0x4.77ffffffda2ccp-4 }, /* i=71 71.733 */
  { 0x1.489e342943062p+0, -0x1.0aed9481b8618p+0, -0x2.917e0ccc376d6p-4, -0x3.d05107b409152p-4, -0x2.66082448af0dep-4, -0x2.e65712f8d46a4p-4, -0x2.bae5089c6da46p-4, 0x4.8800000016538p-4 }, /* i=72 71.68 */
  { 0x1.47931d3fb83abp+0, -0x1.0b407beced595p+0, -0x2.9cfd813c98424p-4, -0x3.da065ea82bb82p-4, -0x2.74b134ccaa99p-4, -0x2.f6fc598a01988p-4, -0x2.d1aa90e84261p-4, 0x4.97ffffffe300cp-4 }, /* i=73 71.626 */
  { 0x1.4687b2b62a177p+0, -0x1.0b94d51c619ddp+0, -0x2.a89a6e5037e5ep-4, -0x3.e3f7021e2dcc8p-4, -0x2.83aed8fa03c04p-4, -0x3.082ca35eb837ep-4, -0x2.e93c3caa9a7bep-4, 0x4.a7fffffff63a4p-4 }, /* i=74 71.572 */
  { 0x1.457bf318fe923p+0, -0x1.0beaa3ca5b872p+0, -0x2.b45587ed4f60ap-4, -0x3.ee2449e2c79e2p-4, -0x2.9303d409fa8cp-4, -0x3.19ecd7d7941bcp-4, -0x3.01a40a66a664ap-4, 0x4.b7fffffffc21cp-4 }, /* i=75 71.516 */
  { 0x1.446fdcf0ce1e2p+0, -0x1.0c41ebc7e1e01p+0, -0x2.c02f86129b0b2p-4, -0x3.f88f99019fd1cp-4, -0x2.a2b30241d09fcp-4, -0x3.2c4214554864cp-4, -0x3.1aea47ad70d62p-4, 0x4.c800000019e58p-4 }, /* i=76 71.456 */
  { 0x1.43636ec2c3571p+0, -0x1.0c9ab0fd19b47p+0, -0x2.cc2924f4615f8p-4, -0x4.033a5e26dea2cp-4, -0x2.b2bf59fe3500cp-4, -0x3.3f31aec103b98p-4, -0x3.351a0aa633ab2p-4, 0x4.d7fffff9db5ecp-4 }, /* i=77 71.403 */
  { 0x1.4256a70e6ccc2p+0, -0x1.0cf4f76a81efdp+0, -0x2.d8432537f523p-4, -0x4.0e26141f60588p-4, -0x2.c32becf6c0daap-4, -0x3.52c1386b0d8e8p-4, -0x3.503d2946e793cp-4, 0x4.e80000000cdfp-4 }, /* i=78 71.35 */
  { 0x1.41498451c1bcep+0, -0x1.0d50c3282299fp+0, -0x2.e47e4be95cdaap-4, -0x4.1954421dc75ep-4, -0x2.d3fbe92dd691ap-4, -0x3.66f6809e24286p-4, -0x3.6c5c4e3a6b944p-4, 0x4.f800000004578p-4 }, /* i=79 71.295 */
  { 0x1.403c0504d76d7p+0, -0x1.0dae186785edcp+0, -0x2.f0db62d1a7cf6p-4, -0x4.24c67c59ee7fp-4, -0x2.e5329a6dfb784p-4, -0x3.7bd797ed2f294p-4, -0x3.8983b00daadap-4, 0x5.0800000199e7p-4 }, /* i=80 71.237 */
  { 0x1.3f2e279c592fep+0, -0x1.0e0cfb7372c78p+0, -0x2.fd5b3880dbc82p-4, -0x4.307e646f7075cp-4, -0x2.f6d36b70aea06p-4, -0x3.916ad33981b3p-4, -0x3.a7bfa906ae634p-4, 0x5.18000000095dp-4 }, /* i=81 71.18 */
  { 0x1.3e1fea87ffd4ap+0, -0x1.0e6d70b107882p+0, -0x3.09fea087bb348p-4, -0x4.3c7da9ebc1628p-4, -0x3.08e1e754819e2p-4, -0x3.a7b6cf29245aep-4, -0x3.c71ae566ea698p-4, 0x5.27ffffffbb014p-4 }, /* i=82 71.12 */
  { 0x1.3d114c3369c94p+0, -0x1.0ecf7ca00785cp+0, -0x3.16c6739722a12p-4, -0x4.48c60ac810ffcp-4, -0x3.1b61bafa35122p-4, -0x3.bec273a26fc8cp-4, -0x3.e7a51f6231fc4p-4, 0x5.37ffffffde72cp-4 }, /* i=83 71.063 */
  { 0x1.3c024b05a180cp+0, -0x1.0f3323dba2f4p+0, -0x3.23b38fb078f7cp-4, -0x4.555953f76ddb8p-4, -0x3.2e56b69126ce8p-4, -0x3.d694f7923f352p-4, -0x4.09676021c0728p-4, 0x5.48000000074fp-4 }, /* i=84 71.002 */
  { 0x1.3af2e56101fd3p+0, -0x1.0f986b1b23d44p+0, -0x3.30c6d853735b2p-4, -0x4.623961f71a398p-4, -0x3.41c4cf2554448p-4, -0x3.ef35e4e719dd8p-4, -0x4.2c755ddc14358p-4, 0x5.580000002823p-4 }, /* i=85 70.946 */
  { 0x1.39e319a321df4p+0, -0x1.0fff57329df85p+0, -0x3.3e0136ad3956cp-4, -0x4.6f68216502c48p-4, -0x3.55b02050f76bcp-4, -0x4.08ad1cbded514p-4, -0x4.50d7be5b90d38p-4, 0x5.6800000020cc8p-4 }, /* i=86 70.882 */
  { 0x1.38d2e624afd24p+0, -0x1.1067ed13aa36ap+0, -0x3.4b6399ca02406p-4, -0x4.7ce78f9dac06cp-4, -0x3.6a1cedf19bb68p-4, -0x4.2302dbd1f71a4p-4, -0x4.76a5abb44b40cp-4, 0x5.780000001f4ap-4 }, /* i=87 70.825 */
  { 0x1.37c249396449fp+0, -0x1.10d231ce1f749p+0, -0x3.58eef6c7812ccp-4, -0x4.8ab9bb6032258p-4, -0x3.7f0fa6082268cp-4, -0x4.3e3fbf2c6f53cp-4, -0x4.9de849b57e75cp-4, 0x5.87ffffffb490cp-4 }, /* i=88 70.768 */
  { 0x1.36b1412fbd2e5p+0, -0x1.113e2a90e6ebp+0, -0x3.66a4490beec4ap-4, -0x4.98e0c57c7f44cp-4, -0x3.948ce29f00208p-4, -0x4.5a6cc917739bcp-4, -0x4.c6b4a74d0d5b8p-4, 0x5.980000000ab1p-4 }, /* i=89 70.708 */
  { 0x1.359fcc5120892p+0, -0x1.11abdcaaae8afp+0, -0x3.74849279f721ap-4, -0x4.a75ee1841ea28p-4, -0x3.aa996bd2b0cep-4, -0x4.7793665d69794p-4, -0x4.f11c144b80da8p-4, 0x5.a80000000449p-4 }, /* i=90 70.646 */
  { 0x1.348de8e17d2d3p+0, -0x1.121b4d8ad51a5p+0, -0x3.8290dbad43868p-4, -0x4.b6365689ff39p-4, -0x3.c13a39f68ae66p-4, -0x4.95bd73d66e9bcp-4, -0x5.1d31e9d0c6914p-4, 0x5.b7fffffff01b4p-4 }, /* i=91 70.583 */
  { 0x1.337b951f42c49p+0, -0x1.128c82c23b213p+0, -0x3.90ca3434e2096p-4, -0x4.c5697fe6799bcp-4, -0x3.d87477d1965b6p-4, -0x4.b4f544487ed5p-4, -0x5.4b0bb3e448ef8p-4, 0x5.c80000000f398p-4 }, /* i=92 70.523 */
  { 0x1.3268cf4349b4ep+0, -0x1.12ff82042105cp+0, -0x3.9f31b2d0c403p-4, -0x4.d4face053f2ap-4, -0x3.f04d84fce06d2p-4, -0x4.d545a6a06977cp-4, -0x5.7abbe7f3fbe74p-4, 0x5.d8000000028f8p-4 }, /* i=93 70.461 */
  { 0x1.3155958099f49p+0, -0x1.137451271a74fp+0, -0x3.adc875b3644dep-4, -0x4.e4ecc73ecb938p-4, -0x4.08caf85f4a00cp-4, -0x4.f6b9ec8e92fd4p-4, -0x5.ac59d94946f1cp-4, 0x5.e8000000181ap-4 }, /* i=94 70.395 */
  { 0x1.3041e6045d564p+0, -0x1.13eaf625f791dp+0, -0x3.bc8fa2c3b4dfcp-4, -0x4.f54208b991538p-4, -0x4.21f2a2cb87d7cp-4, -0x5.195df17dcc82p-4, -0x5.dfffa701a4b9p-4, 0x5.f800000001d18p-4 }, /* i=95 70.336 */
  { 0x1.2f2dbef5a889dp+0, -0x1.14637720c7a0cp+0, -0x3.cb8867e3e2346p-4, -0x5.05fd4758358b4p-4, -0x4.3bca91c784684p-4, -0x5.3d3e2203bf1ccp-4, -0x6.15c38e6fcbaa4p-4, 0x6.07ffffffe5824p-4 }, /* i=96 70.273 */
  { 0x1.2e191e755a64ep+0, -0x1.14ddda5ddae36p+0, -0x3.dab3fb39b7a22p-4, -0x5.172150b1d406cp-4, -0x4.565912722c624p-4, -0x5.626783bdbf14p-4, -0x6.4dbf8b3d74f64p-4, 0x6.1800000011878p-4 }, /* i=97 70.209 */
  { 0x1.2d04029e0651p+0, -0x1.155a264ac80a8p+0, -0x3.ea139b794c7a2p-4, -0x5.28b10c156691cp-4, -0x4.71a4b49475188p-4, -0x5.88e7bda40cd74p-4, -0x6.880f0a1d71cbp-4, 0x6.27fffffffe63cp-4 }, /* i=98 70.142 */
  { 0x1.2bee6983af089p+0, -0x1.15d8617d9002dp+0, -0x3.f9a8903561fa2p-4, -0x5.3aaf7b9cb08acp-4, -0x4.8db44ddf18d04p-4, -0x5.b0cd20e952cc4p-4, -0x6.c4d388733e06p-4, 0x6.38000000025ap-4 }, /* i=99 70.078 */
  { 0x1.2ad85133b22p+0, -0x1.165892b5b5bc3p+0, -0x4.09742a304bffcp-4, -0x5.4d1fbd4960da4p-4, -0x4.aa8efd58bb97cp-4, -0x5.da26b26379dbcp-4, -0x7.042675fd71ba4p-4, 0x6.4800000022028p-4 }, /* i=100 70.012 */
  { 0x1.29c1b7b49be5cp+0, -0x1.16dac0dd6a3c1p+0, -0x4.1977c3b182358p-4, -0x5.60050c319a738p-4, -0x4.c83c2ef1f19d4p-4, -0x6.0504348668ae8p-4, -0x7.462aa050ceff8p-4, 0x6.580000002ec9p-4 }, /* i=101 69.945 */
  { 0x1.28aa9b05f6c7ap+0, -0x1.175ef30ac53c1p+0, -0x4.29b4c0deffbc4p-4, -0x5.7362c1bb8d688p-4, -0x4.e6c39f598d894p-4, -0x6.31763204b07bcp-4, -0x7.8afee166220a8p-4, 0x6.680000001554p-4 }, /* i=102 69.884 */
  { 0x1.2792f920207e8p+0, -0x1.17e5308105f7ep+0, -0x4.3a2c901a07354p-4, -0x5.873c56e841abp-4, -0x5.062d5ffafc4f4p-4, -0x6.5f8e0917ba26cp-4, -0x7.d2ce91d150604p-4, 0x6.77ffffffc4804p-4 }, /* i=103 69.813 */
  { 0x1.267acff40ffabp+0, -0x1.186d80b1e69ccp+0, -0x4.4ae0aa60b4c78p-4, -0x5.9b9565afec1a8p-4, -0x5.2681db49c47f8p-4, -0x6.8f5df784eef3cp-4, -0x8.1db6fb02a97b8p-4, 0x6.87ffffffe0a84p-4 }, /* i=104 69.749 */
  { 0x1.25621d6b45916p+0, -0x1.18f7eb3ee6efep+0, -0x4.5bd293b124e88p-4, -0x5.b071aa6c94f38p-4, -0x5.47c9d92ca328p-4, -0x6.c0f9275c8f4acp-4, -0x8.6beb9f79e598p-4, 0x6.97fffffff986cp-4 }, /* i=105 69.679 */
  { 0x1.2448df677b8b6p+0, -0x1.198477fac1edcp+0, -0x4.6d03db74ad3d4p-4, -0x5.c5d5055bad368p-4, -0x5.6a0e83d094dbcp-4, -0x6.f473bc99313acp-4, -0x8.bd8db0f22eea8p-4, 0x6.a80000004f038p-4 }, /* i=106 69.612 */
  { 0x1.232f13c29537bp+0, -0x1.1a132eead180dp+0, -0x4.7e761cebd8f0cp-4, -0x5.dbc37c2e74f5cp-4, -0x5.8d596c9ca22ccp-4, -0x7.29e2e3921598cp-4, -0x9.12cdb61ed961p-4, 0x6.b7fffffff040cp-4 }, /* i=107 69.545 */
  { 0x1.2214b84e34163p+0, -0x1.1aa41848ac9dcp+0, -0x4.902affa4a903cp-4, -0x5.f2413bb5cf1a4p-4, -0x5.b1b4918e79edp-4, -0x7.615ce0733b7b8p-4, -0x9.6bdd72694573p-4, 0x6.c7ffffffe1cd4p-4 }, /* i=108 69.475 */
  { 0x1.20f9cad3b7ef8p+0, -0x1.1b373c839d872p+0, -0x4.a22437ef412acp-4, -0x6.0952999b9e32p-4, -0x5.d72a62d570b84p-4, -0x7.9af91f9f81abcp-4, -0x9.c8eb58588884p-4, 0x6.d800000019828p-4 }, /* i=109 69.41 */
  { 0x1.1fde4914008dfp+0, -0x1.1bcca4424573ep+0, -0x4.b463875b9cec4p-4, -0x6.20fc1637d0d1p-4, -0x5.fdc5c8c2a14a4p-4, -0x7.d6d0473ae2dd8p-4, -0xa.2a2f48ab82e8p-4, 0x6.e7fffffff0d5cp-4 }, /* i=110 69.335 */
  { 0x1.1ec230c71572p+0, -0x1.1c6458645d9e9p+0, -0x4.c6eabd3e8937cp-4, -0x6.39425e8022178p-4, -0x6.25922a1c556p-4, -0x8.14fc49dcd1898p-4, -0xa.8fea73231751p-4, 0x6.f7fffffff4b6cp-4 }, /* i=111 69.268 */
  { 0x1.1da57f9c0822bp+0, -0x1.1cfe6204698ffp+0, -0x4.d9bbb738b3028p-4, -0x6.522a4e0ce151p-4, -0x6.4e9b72d41108p-4, -0x8.55987a82ab85p-4, -0xa.fa4e07b2bf0fp-4, 0x7.080000004c9p-4 }, /* i=112 69.196 */
  { 0x1.1c883338c5f3dp+0, -0x1.1d9aca7981d9fp+0, -0x4.ecd861c4cd1p-4, -0x6.6bb8f13a4215p-4, -0x6.78ee1b0d20e9cp-4, -0x8.98c1a1ce85e18p-4, -0xb.699daa840366p-4, 0x7.17fffffff9f04p-4 }, /* i=113 69.125 */
  { 0x1.1b6a4939a326cp+0, -0x1.1e399b5956faap+0, -0x5.0042b8d0e92d4p-4, -0x6.85f3876d140acp-4, -0x6.a4972eb0999fp-4, -0x8.de9614d5c0dep-4, -0xb.de1d70d38091p-4, 0x7.27fffffff3284p-4 }, /* i=114 69.054 */
  { 0x1.1a4bbf3156ea1p+0, -0x1.1edade7a08c53p+0, -0x5.13fcc85718824p-4, -0x6.a0df856c44fbp-4, -0x6.d1a4555cab698p-4, -0x9.2735cd5e5e5e8p-4, -0xc.5817223928e08p-4, 0x7.37ffffffcc574p-4 }, /* i=115 68.981 */
  { 0x1.192c92a8929fep+0, -0x1.1f7e9df447dbep+0, -0x5.2808ad03f39d8p-4, -0x6.bc8297e415f2p-4, -0x7.0023dae4d54acp-4, -0x9.72c283e002248p-4, -0xc.d7d3f16dc9cd8p-4, 0x7.47ffffffea114p-4 }, /* i=116 68.911 */
  { 0x1.180cc11ddad6ep+0, -0x1.2024e42566747p+0, -0x5.3c6894e03a12p-4, -0x6.d8e2a605fe7fp-4, -0x7.3024b846b0c68p-4, -0x9.c15fcb467047p-4, -0xd.5da879da75a8p-4, 0x7.57ffffffe905cp-4 }, /* i=117 68.838 */
  { 0x1.16ec48052a33bp+0, -0x1.20cdbbb19c57fp+0, -0x5.511ec00625f9cp-4, -0x6.f605d4528a5ap-4, -0x7.61b69d3f863bp-4, -0xa.13332eae52cb8p-4, -0xd.e9eab4fce00fp-4, 0x7.67ffffffeb144p-4 }, /* i=118 68.764 */
  { 0x1.15cb24c7b45dep+0, -0x1.21792f864e7c8p+0, -0x5.662d815cf16f8p-4, -0x7.13f287875d598p-4, -0x7.94e9fa70bacbp-4, -0xa.68645136c25e8p-4, -0xe.7cf59724f1918p-4, 0x7.77fffffffae7cp-4 }, /* i=119 68.688 */
  { 0x1.14a954c3a1c3bp+0, -0x1.22274adc711fcp+0, -0x5.7b973f5db1efcp-4, -0x7.32af67b775f6cp-4, -0x7.c9d00c28fe458p-4, -0xa.c11d1014fb86p-4, -0xf.172de9fd6423p-4, 0x7.87ffffffb59c4p-4 }, /* i=120 68.613 */
  { 0x1.1386d54ba021ap+0, -0x1.22d8193b1bd46p+0, -0x5.915e74e4cabbp-4, -0x7.5243639460a4cp-4, -0x8.007ae5ea9c32p-4, -0xb.1d89a72463278p-4, -0xf.b8fa751831068p-4, 0x7.98000000340f8p-4 }, /* i=121 68.538 */
  { 0x1.1263a3a6de6f8p+0, -0x1.238ba679f577bp+0, -0x5.a785b20543e94p-4, -0x7.72b5b3dafcc6cp-4, -0x8.38fd7e93a3cap-4, -0xb.7dd8d7f687fap-4, -0x1.062cd246e6acp+0, 0x7.a7ffffffe48dcp-4 }, /* i=122 68.464 */
  { 0x1.113fbd1058181p+0, -0x1.2441fec425bd9p+0, -0x5.be0f9cf3d31b4p-4, -0x7.940ddf059a4a8p-4, -0x8.736bbd7a859bp-4, -0xb.e23c13f8ccfdp-4, -0x1.11515efa000e9p+0, 0x7.b7fffffffb334p-4 }, /* i=123 68.387 */
  { 0x1.101b1eb6d081p+0, -0x1.24fb2e9af69e9p+0, -0x5.d4fef2efb359p-4, -0x7.b653bd2047d44p-4, -0x8.afda883473b48p-4, -0xc.4ae7a96d34208p-4, -0x1.1d0587127d8ep+0, 0x7.c80000000676p-4 }, /* i=124 68.308 */
  { 0x1.0ef5c5bc514efp+0, -0x1.25b742d8df698p+0, -0x5.ec568940a9c88p-4, -0x7.d98f7be3e575cp-4, -0x8.ee5fd16fd91d8p-4, -0xc.b812f3e30d58p-4, -0x1.295124f7c97ccp+0, 0x7.d80000003ef6p-4 }, /* i=125 68.232 */
  { 0x1.0dcfaf35ef7afp+0, -0x1.267648b481c97p+0, -0x6.04194e3bf22bcp-4, -0x7.fdc9a305d004cp-4, -0x9.2f12a8a1b41ep-4, -0xd.29f8902ec46ep-4, -0x1.363e049e0de42p+0, 0x7.e7ffffffcbaacp-4 }, /* i=126 68.154 */
  { 0x1.0ca8d82b329d9p+0, -0x1.27384dc402fc9p+0, -0x6.1c4a4a5dc8d6cp-4, -0x8.230b18d3ad99p-4, -0x9.720b4af4e481p-4, -0xd.a0d6946f08498p-4, -0x1.43d523412b4f6p+0, 0x7.f7fffffff6a04p-4 }, /* i=127 68.075 */
  { 0x1.0b813d95feb08p+0, -0x1.27fd600031d0cp+0, -0x6.34eca16651844p-4, -0x8.495d27044eb5p-4, -0x9.b763351ef59c8p-4, -0xe.1ceecc1215bep-4, -0x1.5220c8aed571fp+0, 0x8.080000001a72p-4 }, /* i=128 67.996 */
  { 0x1.0a58dc620307cp+0, -0x1.28c58dc81f30cp+0, -0x6.4e03938f1558cp-4, -0x8.70c97fe15c2e8p-4, -0x9.ff353696dc068p-4, -0xe.9e86f8a1db368p-4, -0x1.612c29c4f4832p+0, 0x8.180000000f7dp-4 }, /* i=129 67.915 */
  { 0x1.092fb16c4224cp+0, -0x1.2990e5e4ccab4p+0, -0x6.67927ece45de4p-4, -0x8.995a43bd3c68p-4, -0xa.499d861a679cp-4, -0xf.25e917811d738p-4, -0x1.71026c5d36295p+0, 0x8.280000010859p-4 }, /* i=130 67.834 */
  { 0x1.0805b9830e811p+0, -0x1.2a5f778cb0eccp+0, -0x6.819ce0205a5ccp-4, -0x8.c31a06ace94b8p-4, -0xa.96b9d76a05c48p-4, -0xf.b363acd154988p-4, -0x1.81b0dea82a6fbp+0, 0x8.37fffffff1f18p-4 }, /* i=131 67.753 */
  { 0x1.06daf164d8746p+0, -0x1.2b31526835582p+0, -0x6.9c2654fd6e394p-4, -0x8.ee13d6c598598p-4, -0xa.e6a972f701088p-4, -0x1.0474a14c9c252p+0, -0x1.93441d35ee076p+0, 0x8.47ffffffe8958p-4 }, /* i=132 67.67 */
  { 0x1.05af55c07f3edp+0, -0x1.2c06869553dc3p+0, -0x6.b7329cbfaa96cp-4, -0x9.1a53428528a38p-4, -0xb.398d4eb0276e8p-4, -0x1.0e1f4daaad957p+0, -0x1.a5ca92f3e11c3p+0, 0x8.57ffffffa4728p-4 }, /* i=133 67.589 */
  { 0x1.0482e3345ce53p+0, -0x1.2cdf24ac42f69p+0, -0x6.d2c59a39cde7p-4, -0x9.47e45fcf112bp-4, -0xb.8f8828faf2358p-4, -0x1.183c2173ed5aap+0, -0x1.b953da2657c35p+0, 0x8.6800000023ebp-4 }, /* i=134 67.505 */
  { 0x1.0355964e22ff3p+0, -0x1.2dbb3dc3be91fp+0, -0x6.eee355502382cp-4, -0x9.76d3d339abbc8p-4, -0xb.e8bea55c64718p-4, -0x1.22d15d681aac4p+0, -0x1.cdf0460a07ee5p+0, 0x8.77ffffffe8e58p-4 }, /* i=135 67.42 */
  { 0x1.02276b89f6e02p+0, -0x1.2e9ae3760d62ap+0, -0x7.0b8ffcb964c18p-4, -0x9.a72ed7eed282p-4, -0xc.45576b6377d38p-4, -0x1.2de5a85dbb599p+0, -0x1.e3b13364fa1bap+0, 0x8.880000003bafp-4 }, /* i=136 67.337 */
  { 0x1.00f85f524826dp+0, -0x1.2f7e27e5b41cp+0, -0x7.28cfe7c61ff24p-4, -0x9.d90347eed13ep-4, -0xc.a57b4795fa1e8p-4, -0x1.398016ae983e1p+0, -0x1.faa9b4611c082p+0, 0x8.980000004176p-4 }, /* i=137 67.25 */
  { 0xf.fc86dff0d9398p-4, -0x1.30651dc2d07fbp+0, -0x7.46a7984faf95p-4, -0xa.0c5fa4f34eafp-4, -0xd.09554ef90a3c8p-4, -0x1.45a83245eaf41p+0, -0x2.12ee79734040cp+0, 0x8.a7fffffff8df8p-4 }, /* i=138 67.165 */
  { 0xf.e9793d53d1a2p-4, -0x1.314fd85084a52p+0, -0x7.651bbcbca7b9p-4, -0xa.415321db625e8p-4, -0xd.7113052aa67bp-4, -0x1.5266035748ecap+0, -0x2.2c94dd35b78fp+0, 0x8.b7ffffffc7398p-4 }, /* i=139 67.078 */
  { 0xf.d65cd063f8fc8p-4, -0x1.323e6b6aa34b6p+0, -0x7.8431322241f98p-4, -0xa.77edacb528dcp-4, -0xd.dce4853fb458p-4, -0x1.5fc219cf7457ep+0, -0x2.47b47d75ae664p+0, 0x8.c7fffffff7d08p-4 }, /* i=140 66.99 */
  { 0xf.c3315af63356p-4, -0x1.3330eb8b9e115p+0, -0x7.a3ed06844f9f4p-4, -0xa.b03ff96fec59p-4, -0xe.4cfcadc29f138p-4, -0x1.6dc597904e504p+0, -0x2.64664ed2b4a66p+0, 0x8.d7fffffffefe8p-4 }, /* i=141 66.901 */
  { 0xf.aff69d91c949p-4, -0x1.34276dd2e0208p+0, -0x7.c4547b3ad035cp-4, -0xa.ea5b8d4cb7a2p-4, -0xe.c19150086d218p-4, -0x1.7c7a3b8ea772cp+0, -0x2.82c5f27e5affcp+0, 0x8.e800000003b7p-4 }, /* i=142 66.812 */
  { 0xf.9cac5766f703p-4, -0x1.3522080b5339ep+0, -0x7.e56d077639954p-4, -0xb.2652cb0b19e3p-4, -0xf.3adb6326bdd3p-4, -0x1.8bea6de0f203p+0, -0x2.a2f0206fe2138p+0, 0x8.f7fffffff9c28p-4 }, /* i=143 66.722 */
  { 0xf.89524643e38c8p-4, -0x1.3620d0b246484p+0, -0x8.073c5aeb67e7p-4, -0xb.6438ffecae0dp-4, -0xf.b9173ac4e1d18p-4, -0x1.9c214cdb2ebd2p+0, -0x2.c50503ebd7d18p+0, 0x9.07ffffffb7eb8p-4 }, /* i=144 66.632 */
  { 0xf.75e8268838b58p-4, -0x1.3723defeb8854p+0, -0x8.29c860a99967p-4, -0xb.a422719fb1728p-4, -0x1.03c84c256fd27p+0, -0x1.ad2abb5072b7bp+0, -0x2.e925f8a4b93a8p+0, 0x9.17ffffffd1d78p-4 }, /* i=145 66.54 */
  { 0xf.626db31af415p-4, -0x1.382b4ae8de57fp+0, -0x8.4d174216d6438p-4, -0xb.e6246d16b2fdp-4, -0x1.0c567bccd8972p+0, -0x1.bf13701381eb3p+0, -0x3.0f7834c946064p+0, 0x9.28000000398dp-4 }, /* i=146 66.448 */
  { 0xf.4ee2a55e8028p-4, -0x1.39372d321827bp+0, -0x8.712f6a1c6ea2p-4, -0xc.2a55566cc7d48p-4, -0x1.154080976fdaep+0, -0x1.d1e906db3965bp+0, -0x3.38230b8570a16p+0, 0x9.37ffffffe96e8p-4 }, /* i=147 66.354 */
  { 0xf.3b46b51fb5538p-4, -0x1.3a479f6d917e1p+0, -0x8.9617888d6d738p-4, -0xc.70ccb9f86783p-4, -0x1.1e8b1ee701bcap+0, -0x1.e5ba12ae6d533p+0, -0x3.6351be7edc792p+0, 0x9.480000004b9ap-4 }, /* i=148 66.259 */
  { 0xf.2799988d7c908p-4, -0x1.3b5cbc08d1f65p+0, -0x8.bbd695b105a08p-4, -0xc.b9a35e73c9b9p-4, -0x1.283b668ee128p+0, -0x1.fa9631f3b636p+0, -0x3.913153da5d25ap+0, 0x9.580000011df3p-4 }, /* i=149 66.164 */
  { 0xf.13db042b44e38p-4, -0x1.3c769e5500aefp+0, -0x8.e273d60d9e618p-4, -0xd.04f35882400dp-4, -0x1.3256b8403fecep+0, -0x2.108e245b4c126p+0, -0x3.c1f6039355d84p+0, 0x9.680000002654p-4 }, /* i=150 66.069 */
  { 0xf.000aaab9a6b58p-4, -0x1.3d9562912adadp+0, -0x9.09f6de8302188p-4, -0xd.52d81fc3c5648p-4, -0x1.3ce2cb72f8163p+0, -0x2.27b3e2db8982p+0, -0x3.f5d4653c4b8b4p+0, 0x9.780000005858p-4 }, /* i=151 65.972 */
  { 0xe.ec283d3357338p-4, -0x1.3eb925f3e4a1cp+0, -0x9.3267987cb9a58p-4, -0xd.a36ea52086f2p-4, -0x1.47e5b4bb39c98p+0, -0x2.401ab9d187876p+0, -0x4.2d089b30b6f1cp+0, 0x9.8800000002a2p-4 }, /* i=152 65.874 */
  { 0xe.d8336ab413d48p-4, -0x1.3fe206b6a9169p+0, -0x9.5bce46942a128p-4, -0xd.f6d56b0eb8c8p-4, -0x1.5365ecb80ae04p+0, -0x2.59d765b2c4498p+0, -0x4.67d17653f05f8p+0, 0x9.980000004bb8p-4 }, /* i=153 65.776 */
  { 0xe.c42be06f18a3p-4, -0x1.41102420e572bp+0, -0x9.8633895d4d228p-4, -0xe.4d2c9f64aeaf8p-4, -0x1.5f6a578b54c8bp+0, -0x2.7500326626eb6p+0, -0x4.a674da60b3354p+0, 0x9.a7ffffff93068p-4 }, /* i=154 65.676 */
  { 0xe.b0114996acd28p-4, -0x1.42439e94764e5p+0, -0x9.b1a064a141968p-4, -0xe.a69637590a1bp-4, -0x1.6bfa4d00a99fcp+0, -0x2.91ad1db410e8cp+0, -0x4.e93e4f4cd18e4p+0, 0x9.b7ffffff34ae8p-4 }, /* i=155 65.575 */
  { 0xe.9be34f4ef4378p-4, -0x1.437c979a1cb4ap+0, -0x9.de1e44d6b3c8p-4, -0xf.03360d6f5b8ap-4, -0x1.791da160c712bp+0, -0x2.aff7fcff6a6f6p+0, -0x5.307f01a17d97cp+0, 0x9.c7ffffffa4828p-4 }, /* i=156 65.473 */
  { 0xe.87a1989b64bb8p-4, -0x1.44bb31eeeb01ep+0, -0xa.0bb70500bbb38p-4, -0xf.633201c4659d8p-4, -0x1.86dcaf0bdc21ap+0, -0x2.cffca6b88303ap+0, -0x5.7c9171676f43cp+0, 0x9.d800000033b1p-4 }, /* i=157 65.371 */
  { 0xe.734bca491137p-4, -0x1.45ff91929871bp+0, -0xa.3a74f4f79799p-4, -0xf.c6b21ce7229dp-4, -0x1.954060f1b1f3ep+0, -0x2.f1d91ff1b9cc4p+0, -0x5.cdd5e0d91f19cp+0, 0x9.e8000000a6e1p-4 }, /* i=158 65.267 */
  { 0xe.5ee186dba926p-4, -0x1.4749dbd66b9a8p+0, -0xa.6a62e015da318p-4, -0x1.02de0b55bb9c4p+0, -0x1.a4523df3bb35dp+0, -0x3.15adce7fd34bcp+0, -0x6.24b5ab830fb0cp+0, 0x9.f7ffffffee3f8p-4 }, /* i=159 65.162 */
  { 0xe.4a626e726f9bp-4, -0x1.489a376d74d4ap+0, -0xa.9b8c146b44fa8p-4, -0x1.098ea983bc81ap+0, -0x1.b41c7553e44d5p+0, -0x3.3b9db03b619dcp+0, -0x6.81a5d03a674f8p+0, 0xa.08000000671dp-4 }, /* i=160 65.057 */
  { 0xe.35ce1ebacea18p-4, -0x1.49f0cc7cc78dep+0, -0xa.cdfc6a4b70528p-4, -0x1.107ff34ca96bfp+0, -0x1.c4a9ec3cf085dp+0, -0x3.63ce97d295ac4p+0, -0x6.e523598631e44p+0, 0xa.18000000a8dep-4 }, /* i=161 64.95 */
  { 0xe.212432d2deae8p-4, -0x1.4b4dc4ada2061p+0, -0xb.01c04c87f3c4p-4, -0x1.17b50cbdf1819p+0, -0x1.d6064c91d451ep+0, -0x3.8e697000e91a8p+0, -0x7.4fb94945b364p+0, 0xa.280000000edbp-4 }, /* i=162 64.841 */
  { 0xe.0c64432e94f9p-4, -0x1.4cb14b406a7c6p+0, -0xb.36e4c136317cp-4, -0x1.1f314a31c802cp+0, -0x1.e83e151fd238ep+0, -0x3.bb9a85cb578cap+0, -0x7.c1fd2577e27ccp+0, 0xa.380000003a55p-4 }, /* i=163 64.732 */
  { 0xd.f78de584226ep-4, -0x1.4e1b8d203e315p+0, -0xb.6d7772ff5effp-4, -0x1.26f833c35e679p+0, -0x1.fb5eab4b68097p+0, -0x3.eb91da867b524p+0, -0x8.3c95edce3e4cp+0, 0xa.480000001a1cp-4 }, /* i=164 64.622 */
  { 0xd.e2a0acab44df8p-4, -0x1.4f8cb8f87e796p+0, -0xb.a586bb3b3d1cp-4, -0x1.2f0d891164b6bp+0, -0x2.0f766e7b99fc4p+0, -0x4.1e837ebef22dcp+0, -0x8.c0380dbbc9bep+0, 0xa.580000000c97p-4 }, /* i=165 64.51 */
  { 0xd.cd9c28823f488p-4, -0x1.5104ff4b34aa7p+0, -0xb.df21acb4e29cp-4, -0x1.3775454f57ccbp+0, -0x2.2494cd5a18a52p+0, -0x4.54a7f6e8e1b78p+0, -0x9.4daa973b25618p+0, 0xa.68000000924dp-4 }, /* i=166 64.397 */
  { 0xd.b87fe5d38c68p-4, -0x1.52849288b1478p+0, -0xc.1a581f3453e18p-4, -0x1.4033a3aee1a75p+0, -0x2.3aca5d253aap+0, -0x4.8e3cab11bc898p+0, -0x9.e5c77c671132p+0, 0xa.77ffffff63548p-4 }, /* i=167 64.283 */
  { 0xd.a34b6e2b17f88p-4, -0x1.540ba729ba7f7p+0, -0xc.573abc0f06f48p-4, -0x1.494d2431409b6p+0, -0x2.5228f35cc9c5ap+0, -0x4.cb846319172d4p+0, -0xa.89806f34e2de8p+0, 0xa.87ffffffd7168p-4 }, /* i=168 64.168 */
  { 0xd.8dfe47c30bffp-4, -0x1.559a73c98e71cp+0, -0xc.95db0b623c0d8p-4, -0x1.52c690d9d9187p+0, -0x2.6ac3c1e5771d4p+0, -0x5.0cc7d0931da4cp+0, -0xb.39ddfebc099a8p+0, 0xa.980000002c9p-4 }, /* i=169 64.051 */
  { 0xd.7897f55837828p-4, -0x1.573131433bfe7p+0, -0xc.d64b8292038e8p-4, -0x1.5ca503606ec3bp+0, -0x2.84af761ae8532p+0, -0x5.5256288e6d18cp+0, -0xb.f80268aa35ce8p+0, 0xa.a800000003d7p-4 }, /* i=170 63.933 */
  { 0xd.6317f6055cf1p-4, -0x1.58d01ad03e283p+0, -0xd.189f93cd49e08p-4, -0x1.66edeb624a635p+0, -0x2.a0025b02ef3bap+0, -0x5.9c85cee8e8034p+0, -0xc.c531363909ccp+0, 0xa.b800000029d5p-4 }, /* i=171 63.814 */
  { 0xd.4d7dc51f17738p-4, -0x1.5a776e28e8334p+0, -0xd.5cebbec01b078p-4, -0x1.71a71521d762ep+0, -0x2.bcd47f0afec4p+0, -0x5.ebb515865fbap+0, -0xd.a2cc3c037107p+0, 0xa.c80000007fcfp-4 }, /* i=172 63.693 */
  { 0xd.37c8da0bc4ep-4, -0x1.5c276ba71e8f6p+0, -0xd.a345a29d90868p-4, -0x1.7cd6b0e5a90ecp+0, -0x2.db3fddb52347p+0, -0x6.404b120acc92p+0, -0xe.925e34b468c4p+0, 0xa.d7ffffff53288p-4 }, /* i=173 63.571 */
  { 0xd.21f8a80f32fa8p-4, -0x1.5de0566c31dap+0, -0xd.ebc411bcf2e4p-4, -0x1.88835b0c4255ep+0, -0x2.fb608dc3a826ep+0, -0x6.9ab88d33ec13p+0, -0xf.9599114c4ca5p+0, 0xa.e800000009a2p-4 }, /* i=174 63.447 */
  { 0xd.0c0c9e2c189bp-4, -0x1.5fa274875f554p+0, -0xe.367f267caa37p-4, -0x1.94b424d56db05p+0, -0x3.1d54f4301d1f6p+0, -0x6.fb790eafd73ep+0, -0x1.0ae60d2a21c22p+4, 0xa.f7fffffffd5c8p-4 }, /* i=175 63.321 */
  { 0xc.f60426e6da11p-4, -0x1.616e0f21419d9p+0, -0xe.83905a3c10cp-4, -0x1.a1709e11597a6p+0, -0x3.413dfcdb06e3cp+0, -0x7.63140a255ccacp+0, -0x1.1decc2bace9c4p+4, 0xb.080000000ffdp-4 }, /* i=176 63.195 */
  { 0xc.dfdea818b84p-4, -0x1.634372a8b9516p+0, -0xe.d3129dfc0be88p-4, -0x1.aec0dfb44f60cp+0, -0x3.673f5957bcbbap+0, -0x7.d21e312a92574p+0, -0x1.3292e501f59c9p+4, 0xb.1800000025e7p-4 }, /* i=177 63.066 */
  { 0xc.c99b82b89faep-4, -0x1.6522ef03a2634p+0, -0xf.2522752ccc32p-4, -0x1.bcad977769f87p+0, -0x3.8f7fc6c5a6786p+0, -0x8.493aefc1788e8p+0, -0x1.4901f5629cae3p+4, 0xb.280000001985p-4 }, /* i=178 62.936 */
  { 0xc.b33a12a166e18p-4, -0x1.670cd7c2fc64fp+0, -0xf.79de12c33186p-4, -0x1.cb40149fec3dbp+0, -0x3.ba295b8f59408p+0, -0x8.c91e192e2755p+0, -0x1.616805ebdf193p+4, 0xb.380000003d47p-4 }, /* i=179 62.804 */
  { 0xc.9cb9ae573ffb8p-4, -0x1.6901845aa491p+0, -0xf.d16578c801dbp-4, -0x1.da82560d50653p+0, -0x3.e769de07c8d22p+0, -0x9.528dcbf8b685p+0, -0x1.7bf8a595e7c32p+4, 0xb.480000000657p-4 }, /* i=180 62.671 */
  { 0xc.8619a6c4c482p-4, -0x1.6b01505d970f5p+0, -0x1.02bda9abdde68p+0, -0x1.ea7f19b8b70adp+0, -0x4.1773252228c74p+0, -0x9.e6649543d1948p+0, -0x1.98ed63f054b21p+4, 0xb.5800000020f8p-4 }, /* i=181 62.535 */
  { 0xc.6f5946fa3e238p-4, -0x1.6d0c9bbe240eap+0, -0x1.0896182fcc849p+0, -0x1.fb41edc94a51p+0, -0x4.4a7b847c3f81cp+0, -0xa.8593dc6fdef8p+0, -0x1.b886c8c4d5608p+4, 0xb.68000000278bp-4 }, /* i=182 62.398 */
  { 0xc.5877d3e308fdp-4, -0x1.6f23cb13866d8p+0, -0x1.0ea207b824cf7p+0, -0x2.0cd743760f186p+0, -0x4.80be455c78028p+0, -0xb.3126a21a5862p+0, -0x1.db0cf573fdcfep+4, 0xb.7800000069f5p-4 }, /* i=183 62.259 */
  { 0xc.41748bfbeafa8p-4, -0x1.714747e450ae1p+0, -0x1.14e403a64884bp+0, -0x2.1f4c83e43bdbp+0, -0x4.ba7c2e431d474p+0, -0xb.ea449eafcbe48p+0, -0x2.00d1167250cf2p+4, 0xb.87ffffffdd438p-4 }, /* i=184 62.118 */
  { 0xc.2a4ea6fbf6a38p-4, -0x1.737780f775926p+0, -0x1.1b5ec1289036ap+0, -0x2.32b0274b4de2cp+0, -0x4.f7fc1b3a2fb44p+0, -0xc.b235cf99e9d88p+0, -0x2.2a2e65dc98f1cp+4, 0xb.980000000c4bp-4 }, /* i=185 61.975 */
  { 0xc.1305558526a98p-4, -0x1.75b4eaaa78e78p+0, -0x1.2215228b4bd63p+0, -0x2.4711ce94e43fcp+0, -0x5.398ba8fadbcbp+0, -0xd.8a66837acbaf8p+0, -0x2.578b8a79cbbcp+4, 0xb.a800000005d4p-4 }, /* i=186 61.831 */
  { 0xb.fb97c0ced40dp-4, -0x1.77ffff4e90137p+0, -0x1.290a3adaa9d69p+0, -0x2.5c825fd298e6p+0, -0x5.7f7ff5b00d0ecp+0, -0xe.746bf9cb3b598p+0, -0x2.895c879b1b5b6p+4, 0xb.b7fffff7d0848p-4 }, /* i=187 61.684 */
  { 0xb.e4050a09f1048p-4, -0x1.7a593f93c2a1bp+0, -0x1.304151f251a6ep+0, -0x2.731426249575cp+0, -0x5.ca367b453ab08p+0, -0xf.7209b02a2582p+0, -0x2.c0245764ffebap+4, 0xb.c7ffffff70458p-4 }, /* i=188 61.534 */
  { 0xb.cc4c4a6c05538p-4, -0x1.7cc132eb4a884p+0, -0x1.37bde8bd19f81p+0, -0x2.8adaf49baf224p+0, -0x6.1a1603f5b78b8p+0, -0x1.085376f69fa61p+4, -0x2.fc776922c34dp+4, 0xb.d7ffffffe7128p-4 }, /* i=189 61.383 */
  { 0xb.b46c924e2e73p-4, -0x1.7f386809bc4fbp+0, -0x1.3f83be29922a7p+0, -0x2.a3ec4e08a85b4p+0, -0x6.6f8fc228d64ccp+0, -0x1.1b0284955276ap+4, -0x3.3efddaa411b94p+4, 0xb.e80000000e16p-4 }, /* i=190 61.229 */
  { 0xb.9c64e8eb3ce7p-4, -0x1.81bf7562d8fefp+0, -0x1.4796d46955bap+0, -0x2.be5f90d2c9c28p+0, -0x6.cb208bb6b409cp+0, -0x1.2f5529ded2291p+4, -0x3.887707a8ec552p+4, 0xb.f7ffffffc86d8p-4 }, /* i=191 61.073 */
  { 0x1.71f93c6d63c51p+0, -0x2.1f53124f173dp-4, 0x8.972b8c4a20448p-8, -0x2.e4f3dcbb95cb6p-8, 0x1.1e4f5193c9038p-8, -0x7.7e6f5037aee1p-12, 0xc.07ffffffe3c58p-4 }, /* i=192 69.804 */
  { 0x1.71d74fd079481p+0, -0x2.1e40b7836dd52p-4, 0x8.8e8361e012ap-8, -0x2.e07f4aa881dc4p-8, 0x1.1bfae59d98d05p-8, -0x7.6ab46d8822b6p-12, 0xc.180000062602p-4 }, /* i=193 69.821 */
  { 0x1.71b57450b185bp+0, -0x2.1d2f70e89a9d2p-4, 0x8.85e887448c588p-8, -0x2.dc13fdfffce02p-8, 0x1.19ac9a3a4d121p-8, -0x7.573852d0b665cp-12, 0xc.2800000040d7p-4 }, /* i=194 69.838 */
  { 0x1.7193a9dcb1183p+0, -0x2.1c1f3cd52faeap-4, 0x8.7d5ae0c25bb08p-8, -0x2.d7b1de61e151p-8, 0x1.17645bed0ee29p-8, -0x7.43fa15c498a34p-12, 0xc.37fffffffcbe8p-4 }, /* i=195 69.857 */
  { 0x1.7171f06367b2p+0, -0x2.1b1019a4b978cp-4, 0x8.74da52f953a1p-8, -0x2.d358d3c1bfedcp-8, 0x1.1522178791ebcp-8, -0x7.30f8d013572f4p-12, 0xc.480000002707p-4 }, /* i=196 69.872 */
  { 0x1.715047d3d6663p+0, -0x2.1a0205b5e4eecp-4, 0x8.6c66c2cebd1f8p-8, -0x2.cf08c65e35b18p-8, 0x1.12e5ba1c65e02p-8, -0x7.1e339f31034ecp-12, 0xc.57ffffffb4608p-4 }, /* i=197 69.891 */
  { 0x1.712eb01d12eabp+0, -0x2.18f4ff6a93158p-4, 0x8.6400156d5d68p-8, -0x2.cac19ec04a7a6p-8, 0x1.10af310c6e4f7p-8, -0x7.0ba9a47253cf4p-12, 0xc.6800000023c9p-4 }, /* i=198 69.908 */
  { 0x1.710d292e52438p+0, -0x2.17e90528277fp-4, 0x8.5ba6304751768p-8, -0x2.c68345bbbf5aap-8, 0x1.0e7e69f7eb7cfp-8, -0x6.f95a04c347eecp-12, 0xc.78000000529fp-4 }, /* i=199 69.925 */
  { 0x1.70ebb2f6e1063p+0, -0x2.16de1557432bp-4, 0x8.5358f9134d4c8p-8, -0x2.c24da46d0594ep-8, 0x1.0c5352c06ae07p-8, -0x6.e743e8c00ed5cp-12, 0xc.880000001087p-4 }, /* i=200 69.942 */
  { 0x1.70ca4d66233aep+0, -0x2.15d42e63bcc4ep-4, 0x8.4b1855cbca8cp-8, -0x2.be20a4383606p-8, 0x1.0a2dd99306652p-8, -0x6.d5667c8fbba3p-12, 0xc.980000001774p-4 }, /* i=201 69.961 */
  { 0x1.70a8f86b98b8p+0, -0x2.14cb4ebcbc85ep-4, 0x8.42e42caf539ep-8, -0x2.b9fc2ec89ae2p-8, 0x1.080decd700a23p-8, -0x6.c3c0efd7fdf5cp-12, 0xc.a7ffffffdcd38p-4 }, /* i=202 69.978 */
  { 0x1.7087b3f6d7d5p+0, -0x2.13c374d48b368p-4, 0x8.3abc643e6dd3p-8, -0x2.b5e02e0f07e3ep-8, 0x1.05f37b3589c31p-8, -0x6.b25275b70af1cp-12, 0xc.b80000000307p-4 }, /* i=203 69.994 */
  { 0x1.70667ff7922dcp+0, -0x2.12bc9f20b157ep-4, 0x8.32a0e33c02058p-8, -0x2.b1cc8c417804ep-8, 0x1.03de739bfc68dp-8, -0x6.a11a449bed2bcp-12, 0xc.c7fffffff7858p-4 }, /* i=204 70.011 */
  { 0x1.70455c5d9026ep+0, -0x2.11b6cc19ccf26p-4, 0x8.2a9190ab821e8p-8, -0x2.adc133d989a7ep-8, 0x1.01cec52ee141ep-8, -0x6.901796493c1dcp-12, 0xc.d8000000021fp-4 }, /* i=205 70.028 */
  { 0x1.70244918b34fcp+0, -0x2.10b1fa3b9df52p-4, 0x8.228e53d0c0b6p-8, -0x2.a9be0f93d84b2p-8, 0xf.fc45f545aabf8p-12, -0x6.7f49a7ba9da44p-12, 0xc.e80000000789p-4 }, /* i=206 70.046 */
  { 0x1.70034618f5379p+0, -0x2.0fae2804f6838p-4, 0x8.1a97142eecc8p-8, -0x2.a5c30a6eea076p-8, 0xf.dbf31b07f00fp-12, -0x6.6eafb92d00a74p-12, 0xc.f8000000017p-4 }, /* i=207 70.063 */
  { 0x1.6fe2534e65ee6p+0, -0x2.0eab53f7a8daep-4, 0x8.12abb9877d658p-8, -0x2.a1d00faa180d4p-8, 0xf.bbf2c20d19ec8p-12, -0x6.5e490dea020e8p-12, 0xd.08000000a1p-4 }, /* i=208 70.078 */
  { 0x1.6fc170a932b07p+0, -0x2.0da97c98b5602p-4, 0x8.0acc2bdb15678p-8, -0x2.9de50ac574af4p-8, 0xf.9c43eb9e32928p-12, -0x6.4e14ec4b84b24p-12, 0xd.17fffffff9918p-4 }, /* i=209 70.096 */
  { 0x1.6fa09e1997e16p+0, -0x2.0ca8a06fd665ap-4, 0x8.02f8536573568p-8, -0x2.9a01e77f3c5dap-8, 0xf.7ce59d21b8b8p-12, -0x6.3e129dbe1c1fcp-12, 0xd.2800000039efp-4 }, /* i=210 70.114 */
  { 0x1.6f7fdb8ff0e27p+0, -0x2.0ba8be07f63bap-4, 0x7.fb3018a08bdp-8, -0x2.962691d4da7c4p-8, 0xf.5dd6df0e03fb8p-12, -0x6.2e416e9a0d6b4p-12, 0xd.37ffffffe12d8p-4 }, /* i=211 70.131 */
  { 0x1.6f5f28fcab6a7p+0, -0x2.0aa9d3eec5fap-4, 0x7.f3736440d4e88p-8, -0x2.9252f6008d6f8p-8, 0xf.3f16bd835d5a8p-12, -0x6.1ea0ae1de35dcp-12, 0xd.47ffffffd1938p-4 }, /* i=212 70.148 */
  { 0x1.6f3e86504d97dp+0, -0x2.09abe0b4e70e2p-4, 0x7.ebc21f360d9ap-8, -0x2.8e87007947b4p-8, 0xf.20a448166881p-12, -0x6.0f2fae61bfd94p-12, 0xd.5800000044f3p-4 }, /* i=213 70.165 */
  { 0x1.6f1df37b7785fp+0, -0x2.08aee2edf17cep-4, 0x7.e41c32aaf19ap-8, -0x2.8ac29df209558p-8, 0xf.027e91b25ba68p-12, -0x5.ffedc452bd8f4p-12, 0xd.67ffffffed178p-4 }, /* i=214 70.181 */
  { 0x1.6efd706edaac7p+0, -0x2.07b2d9302afbp-4, 0x7.dc8188028afdp-8, -0x2.8705bb580f63ap-8, 0xe.e4a4b0799de8p-12, -0x5.f0da478ce6bbp-12, 0xd.78000000198bp-4 }, /* i=215 70.198 */
  { 0x1.6edcfd1b440abp+0, -0x2.06b7c214d0254p-4, 0x7.d4f208d9f0fa4p-8, -0x2.835045d333f76p-8, 0xe.c715be021b258p-12, -0x5.e1f4926b75b44p-12, 0xd.87ffffffa2bb8p-4 }, /* i=216 70.215 */
  { 0x1.6ebc9971913b1p+0, -0x2.05bd9c37ba2p-4, 0x7.cd6d9f0517c8cp-8, -0x2.7fa22ac3e308ap-8, 0xe.a9d0d6f824c78p-12, -0x5.d33c01f27e5acp-12, 0xd.9800000034bbp-4 }, /* i=217 70.231 */
  { 0x1.6e9c4562bce7ap+0, -0x2.04c46637b948cp-4, 0x7.c5f4349116ba8p-8, -0x2.7bfb57c3bfae8p-8, 0xe.8cd51b402fddp-12, -0x5.c4aff5b121748p-12, 0xd.a8000000147cp-4 }, /* i=218 70.247 */
  { 0x1.6e7c00dfd2646p+0, -0x2.03cc1eb62fd4ep-4, 0x7.be85b3c0a8d44p-8, -0x2.785bbaa377514p-8, 0xe.7021adb144a7p-12, -0x5.b64fcfc7c150cp-12, 0xd.b800000029p-4 }, /* i=219 70.265 */
  { 0x1.6e5bcbd9f4eep+0, -0x2.02d4c4574446cp-4, 0x7.b722070d3fbb8p-8, -0x2.74c3416ad4246p-8, 0xe.53b5b4686703p-12, -0x5.a81af4ece85p-12, 0xd.c800000051cep-4 }, /* i=220 70.281 */
  { 0x1.6e3ba6425edefp+0, -0x2.01de55c1d5c2p-4, 0x7.afc919263969p-8, -0x2.7131da57e7854p-8, 0xe.3790586f3532p-12, -0x5.9a10cc368766cp-12, 0xd.d7ffffffb41d8p-4 }, /* i=221 70.297 */
  { 0x1.6e1b900a5a905p+0, -0x2.00e8d19f401dap-4, 0x7.a87ad4eea8b38p-8, -0x2.6da773dd85c9ep-8, 0xe.1bb0c5c26405p-12, -0x5.8c30bf27f154p-12, 0xd.e800000031bcp-4 }, /* i=222 70.314 */
  { 0x1.6dfb89234e8ddp+0, -0x1.fff4369bb4117p-4, 0x7.a137257f892b4p-8, -0x2.6a23fca3e8ff8p-8, 0xe.00162b4ca90a8p-12, -0x5.7e7a39b94478p-12, 0xd.f8000000406ap-4 }, /* i=223 70.33 */
  { 0x1.6ddb917eb3508p+0, -0x1.ff008365e33e6p-4, 0x7.99fdf624d327cp-8, -0x2.66a76386d8cf6p-8, 0xd.e4bfbac12daep-12, -0x5.70ecaa1dee078p-12, 0xe.07fffffff1e08p-4 }, /* i=224 70.347 */
  { 0x1.6dbba90e13891p+0, -0x1.fe0db6aefd6a9p-4, 0x7.92cf325cfea88p-8, -0x2.63319795036fep-8, 0xd.c9aca8c809ebp-12, -0x5.638780ddb9bp-12, 0xe.18000000bdap-4 }, /* i=225 70.363 */
  { 0x1.6d9bcfc31739cp+0, -0x1.fd1bcf2affcdap-4, 0x7.8baac5daf31bp-8, -0x2.5fc2881082754p-8, 0xd.aedc2cdb77788p-12, -0x5.564a30cd04eb4p-12, 0xe.27ffffffda998p-4 }, /* i=226 70.379 */
  { 0x1.6d7c058f6fabfp+0, -0x1.fc2acb901799p-4, 0x7.84909c80f8bacp-8, -0x2.5c5a246c0386ap-8, 0xd.944d80ed6acp-12, -0x5.49342edeaa6dp-12, 0xe.38000000392ep-4 }, /* i=227 70.396 */
  { 0x1.6d5c4a64eb89bp+0, -0x1.fb3aaa9735b61p-4, 0x7.7d80a264b07ap-8, -0x2.58f85c4c4abe8p-8, 0xd.79ffe1ec9b718p-12, -0x5.3c44f24f347c8p-12, 0xe.48000000578p-4 }, /* i=228 70.412 */
  { 0x1.6d3c9e356adc1p+0, -0x1.fa4b6afbae969p-4, 0x7.767ac3cbd8aa8p-8, -0x2.559d1f864017ap-8, 0xd.5ff28f481efp-12, -0x5.2f7bf46169da4p-12, 0xe.5800000034eap-4 }, /* i=229 70.428 */
  { 0x1.6d1d00f2e3d2ap+0, -0x1.f95d0b7b5983ap-4, 0x7.6f7eed2cd6438p-8, -0x2.52485e1ece198p-8, 0xd.4624cb0eb8368p-12, -0x5.22d8b07eecd6p-12, 0xe.67fffffed5cc8p-4 }, /* i=230 70.445 */
  { 0x1.6cfd728f57671p+0, -0x1.f86f8ad63635fp-4, 0x7.688d0b2bae1b4p-8, -0x2.4efa08490ede4p-8, 0xd.2c95d9dde21bp-12, -0x5.165aa41e7955cp-12, 0xe.77fffffff1688p-4 }, /* i=231 70.461 */
  { 0x1.6cddf2fce8cddp+0, -0x1.f782e7cf189dap-4, 0x7.61a50a9eabd2cp-8, -0x2.4bb20e68234a6p-8, 0xd.13450303aeacp-12, -0x5.0a014ebf75ccp-12, 0xe.87fffffffe5b8p-4 }, /* i=232 70.477 */
  { 0x1.6cbe822dc67a3p+0, -0x1.f697212af6f8ep-4, 0x7.5ac6d888cdaf4p-8, -0x2.4870610c2bc72p-8, 0xc.fa31900078328p-12, -0x4.fdcc31d67acbp-12, 0xe.97ffffffc3098p-4 }, /* i=233 70.493 */
  { 0x1.6c9f201432a22p+0, -0x1.f5ac35b125683p-4, 0x7.53f2621b267ccp-8, -0x2.4534f0f2944f2p-8, 0xc.e15acd0fde2c8p-12, -0x4.f1bad0d9da49p-12, 0xe.a7fffffff4148p-4 }, /* i=234 70.51 */
  { 0x1.6c7fcca286aafp+0, -0x1.f4c2242b6adccp-4, 0x7.4d2794b51c5e8p-8, -0x2.41ffaf05d5f7p-8, 0xc.c8c008b0ade08p-12, -0x4.e5ccb11ad575cp-12, 0xe.b7ffffff60bd8p-4 }, /* i=235 70.525 */
  { 0x1.6c6087cb29418p+0, -0x1.f3d8eb65b2754p-4, 0x7.46665de1c5228p-8, -0x2.3ed08c5bdd24ap-8, 0xc.b06093bb5fbap-12, -0x4.da0159d5ea39cp-12, 0xe.c7ffffffee1b8p-4 }, /* i=236 70.541 */
  { 0x1.6c4151809b324p+0, -0x1.f2f08a2e66ea4p-4, 0x7.3faeab5a33628p-8, -0x2.3ba77a36c6c3cp-8, 0xc.983bc193ca8bp-12, -0x4.ce585430d738p-12, 0xe.d7fffffff4578p-4 }, /* i=237 70.556 */
  { 0x1.6c2229b56bc1ep+0, -0x1.f208ff561717cp-4, 0x7.39006b02772dcp-8, -0x2.38846a031eaf2p-8, 0xc.8050e79b894c8p-12, -0x4.c2d12afef9eep-12, 0xe.e80000002114p-4 }, /* i=238 70.574 */
  { 0x1.6c03105c3fb8bp+0, -0x1.f12249afa616cp-4, 0x7.325b8aeaab1a8p-8, -0x2.35674d5808508p-8, 0xc.689f5dc228bp-12, -0x4.b76b6b00c89f4p-12, 0xe.f7ffffffc9918p-4 }, /* i=239 70.589 */
  { 0x1.6be40567cbc8bp+0, -0x1.f03c68101d418p-4, 0x7.2bbff94d4a78p-8, -0x2.325015f621ce6p-8, 0xc.51267e1971bf8p-12, -0x4.ac26a2a62a35p-12, 0xf.07ffffffb8ca8p-4 }, /* i=240 70.606 */
  { 0x1.6bc508cad8c84p+0, -0x1.ef57594ec74e8p-4, 0x7.252da48fa4678p-8, -0x2.2f3eb5c7668cep-8, 0xc.39e5a4d02c578p-12, -0x4.a1026225c45c8p-12, 0xf.17fffffff0eb8p-4 }, /* i=241 70.621 */
  { 0x1.6ba61a784300fp+0, -0x1.ee731c4526d7ap-4, 0x7.1ea47b414191cp-8, -0x2.2c331ede94614p-8, 0xc.22dc30677c828p-12, -0x4.95fe3b5e0f83p-12, 0xf.27fffffff3028p-4 }, /* i=242 70.636 */
  { 0x1.6b873a62f7b52p+0, -0x1.ed8fafcedfd4fp-4, 0x7.18246c1aec09p-8, -0x2.292d437666b6ap-8, 0xc.0c098156e721p-12, -0x4.8b19c1de2dc7p-12, 0xf.37fffffff4c78p-4 }, /* i=243 70.652 */
  { 0x1.6b68687df638fp+0, -0x1.ecad12c9bb975p-4, 0x7.11ad65fe7b1b8p-8, -0x2.262d15f12e972p-8, 0xb.f56cfa4c1af58p-12, -0x4.80548ad59cabcp-12, 0xf.480000003e4fp-4 }, /* i=244 70.668 */
  { 0x1.6b49a4bc4fb6p+0, -0x1.ebcb4415a2f0cp-4, 0x7.0b3f57f65805p-8, -0x2.233288d84abe8p-8, 0xb.df06001451878p-12, -0x4.75ae2d1c3ae4cp-12, 0xf.5800000137bep-4 }, /* i=245 70.684 */
  { 0x1.6b2aef112d00ep+0, -0x1.eaea4294c4d8ap-4, 0x7.04da313647304p-8, -0x2.203d8edc36eep-8, 0xb.c8d3f979c1afp-12, -0x4.6b26411818c7cp-12, 0xf.680000004d18p-4 }, /* i=246 70.699 */
  { 0x1.6b0c476fbe603p+0, -0x1.ea0a0d2b1b796p-4, 0x6.fe7de1179a9dp-8, -0x2.1d4e1ad27d2c8p-8, 0xb.b2d64f1effe88p-12, -0x4.60bc60b0909cp-12, 0xf.780000004f59p-4 }, /* i=247 70.716 */
  { 0x1.6aedadcb4d608p+0, -0x1.e92aa2beeaea9p-4, 0x6.f82a571c7f1e8p-8, -0x2.1a641fb6f1cfap-8, 0xb.9d0c6be17ffcp-12, -0x4.56702775473b4p-12, 0xf.880000001af1p-4 }, /* i=248 70.731 */
  { 0x1.6acf221732c71p+0, -0x1.e84c023873904p-4, 0x6.f1df82ed7acf4p-8, -0x2.177f90aa40cecp-8, 0xb.8775bc756c568p-12, -0x4.4c41325ac273p-12, 0xf.97ffffffc27d8p-4 }, /* i=249 70.747 */
  { 0x1.6ab0a446d8578p+0, -0x1.e76e2a81fb3cfp-4, 0x6.eb9d545964758p-8, -0x2.14a060f1a14p-8, 0xb.7211af78d8788p-12, -0x4.422f1fdadf3a4p-12, 0xf.a7ffffff921e8p-4 }, /* i=250 70.761 */
  { 0x1.6a92344db8424p+0, -0x1.e6911a87c53d3p-4, 0x6.e563bb54dfc84p-8, -0x2.11c683f651252p-8, 0xb.5cdfb55752cp-12, -0x4.38398fed3149p-12, 0xf.b80000004b85p-4 }, /* i=251 70.778 */
  { 0x1.6a73d21f61d8p+0, -0x1.e5b4d13830a31p-4, 0x6.df32a7faedf1cp-8, -0x2.0ef1ed4590aeep-8, 0xb.47df409951498p-12, -0x4.2e6023e6d831p-12, 0xf.c80000009bc8p-4 }, /* i=252 70.792 */
  { 0x1.6a557daf74927p+0, -0x1.e4d94d8390583p-4, 0x6.d90a0a8b846ep-8, -0x2.0c22908fb6d7p-8, 0xb.330fc54ef2078p-12, -0x4.24a27e8c56c24p-12, 0xf.d7ffffffa51a8p-4 }, /* i=253 70.809 */
  { 0x1.6a3736f19a9c5p+0, -0x1.e3fe8e5bfff7fp-4, 0x6.d2e9d36a10a58p-8, -0x2.095861a73faeep-8, 0xb.1e70b9983501p-12, -0x4.1b00440c49fep-12, 0xf.e7ffffffd3a88p-4 }, /* i=254 70.824 */
  { 0x1.6a18fde8aede9p+0, -0x1.e32493228aa7ep-4, 0x6.ccd1f6299285cp-8, -0x2.069355e2ce8eep-8, 0xb.0a019f6843748p-12, -0x4.11791ea2793p-12, 0xf.f7fff7fffd878p-4 }, /* i=255 70.84 */
};

static const double pi_hi = 0x1.921fb54442d18p1;
static const double pi_lo = 0x1.1a62633145c07p-53;

/* For the slow path, use polynomials of degree DEGREE (thus with DEGREE+1
   coefficients), with coefficients of degree < LARGE represented as
   double-double, and coefficients of degree >= LARGE as double only.
   Thus each polynomial needs DEGREE+LARGE+1 'doubles'. */
#define DEGREE 11
#define LARGE 8

/* Each entry contains a degree 11-polynomial and the evaluation point xmid,
   so that:
   * for |x| < 0.5, i.e., 0 <= i < 64, acos(x)-pi/2 ~ p(x-xmid);
   * for 0.5 <= |x| < 1, acos(1-x) ~ sqrt(x)*p(x-xmid).
   The coefficients of degree < LARGE are stored as two double numbers (most
   significant first).
   Generated with:
   Gen_P_aux(11,128,threshold=64,large=8,verbose=true,slow=1,prec=128). */
static const double T2[128][DEGREE+LARGE+2] = {
  { 0x7.3f7b5c6f3abcp-120, 0x6.aa86b9baf9cab614bp-176, -0x1p+0, -0x4.045deep-104, 0x5.e6990ab0c2dd8p-92, -0x7.40b83a2b1b71454c98p-148, -0x2.aaaaaaaaaaaaap-4, -0xa.aaaae0ef10118e456p-60, 0x1.0160b49131a62p-68, -0x7.963a9045d1a6daedf2p-124, -0x1.3333333333333p-4, -0x6.015892f570086cb23ep-60, 0x4.f0883d2b50248p-52, -0xc.e10e4b3894d2b935c8p-108, -0xb.6db6db6e11cd8p-8, 0x2.bcd79079f2c4abcc19p-60, 0x4.4f35e4ec4bd2p-36, -0x7.c71c9380dc4fcp-8, 0x9.eba442f49964p-24, -0x5.bb98b11f1845p-8, 0x0p+0 }, /* i=0 117.1 */
  { -0x3.0004801239fb4p-8, 0x1.c7305887fe50696224p-64, -0x1.0004801e60e3dp+0, -0x7.023ec56b5ce7e77d04p-56, -0x1.801440e3d9586p-8, -0x6.9ef9beeedeb98a93f8p-64, -0x2.ab16b2941ef82p-4, 0xb.58b641169bbd0d67e8p-60, -0x1.202a33a27b15fp-8, -0x2.6b29211b91f914266p-68, -0x1.33ba46964a456p-4, 0x4.828c4ec205b6156d5ap-60, -0xf.044f1d0a0794p-12, 0x3.3528852e9ea44128f5p-64, -0xb.779131c48b988p-8, 0x2.06e1efb066f0b5a5fdp-60, -0xd.263c033a0cbbp-12, -0x7.d2335f752895p-8, -0xb.d88311b4bd9ep-12, -0x5.c672cfed3749p-8, 0x3p-8 }, /* i=1 115.6 */
  { -0x5.0014d63fc2f6p-8, 0x9.a9afbebc108d40843p-64, -0x1.000c80ea73147p+0, -0x1.603503f2bf69be97eep-56, -0x2.805dcb72fdedep-8, -0x7.ac245975dcd852336p-64, -0x2.abd6e7bd684fcp-4, 0x7.88b3d8f7ca980ee054p-60, -0x1.e0c37ec3b2803p-8, -0x5.22092fe7b71318476cp-64, -0x1.34aac8e3e936fp-4, -0x1.77055e8d9b0c603406p-60, -0x1.913f80fa00244p-8, -0x6.0cc691b2f88cf6cdcap-64, -0xb.8920e8cb460dp-8, 0x2.12732af294d8d7d6fp-60, -0x1.5fce7e060e262p-8, -0x7.e5fdcce61d648p-8, -0x1.3d71c9158346fp-8, -0x5.dc43c52ba8a18p-8, 0x5p-8 }, /* i=2 114.8 */
  { -0x7.00392f97c0c1p-8, -0x1.8d9822a9352bd84d68p-60, -0x1.00188384efb56p+0, -0x1.9e38be72f1d290f8fep-56, -0x3.81017d9a51848p-8, 0xf.fcc8c7cd6f1c44381cp-64, -0x2.acf7956d0c9c2p-4, 0x7.ddc8ba05bfed073b6cp-60, -0x2.a218ebb0302p-8, -0xe.00782f64664d278648p-64, -0x1.361472d1cef2cp-4, 0x5.0c23d206cb3877fa52p-60, -0x2.336e056048622p-8, -0x1.048359d7b3775c4604p-64, -0xb.a3944f235132p-8, -0x1.12ae0f735fdf980c9p-60, -0x1.eef7e59d50c2bp-8, -0x8.03dc4d4a96668p-8, -0x1.bfb5e8527c3b4p-8, -0x5.fd45609a9ba9p-8, 0x7p-8 }, /* i=3 114.3 */
  { -0x9.0079914fef96p-8, 0x7.9919e3720f5e815e2p-64, -0x1.0028899ee96f2p+0, 0x1.f02d504dc33579b604p-56, -0x4.8223988e5f804p-8, 0x1.24e896e2c72391bde8p-64, -0x2.ae792caf62116p-4, -0x5.23e2c1d0e9d717da98p-60, -0x3.6476851ae3a1p-8, -0xa.bf28a954b0a89afd8p-68, -0x1.37f85a3916f82p-4, -0x7.37fa8a155055dd7a5ep-60, -0x2.d74dd09a504ap-8, -0xf.065d667385196c9ee8p-64, -0xb.c70d046f05548p-8, -0x1.e9e639fb20fcf44c84p-60, -0x2.809742069052p-8, -0x8.2c07a71b58808p-8, -0x2.454e4a438e228p-8, -0x6.29cf5a1f56cp-8, 0x9p-8 }, /* i=4 113.9 */
  { -0xb.00de0491731ap-8, -0x2.33c42de1b15e0e5137p-60, -0x1.003c957ad60cbp+0, 0x2.e8ebcafcc589a16c38p-56, -0x5.83e88f0f3777p-8, 0x1.0ca00567fa6daae4dp-60, -0x2.b05c44c66c612p-4, -0x9.43d6bdadd054bcc944p-60, -0x4.282920e5796c4p-8, 0xc.a9d711debd8dcd7f08p-64, -0x1.3a57f3b4be11ap-4, -0x1.0b22af1783f22e176ap-60, -0x3.7d5e5bcc9102p-8, -0xe.b26e9187cbbc2f8c64p-64, -0xb.f3b840a1318ep-8, -0x1.9da6e2fa8461507dcap-60, -0x3.15688dc464eeap-8, -0x8.5ecc775a3bfb8p-8, -0x2.cf3d87d2daedcp-8, -0x6.6258878f46474p-8, 0xbp-8 }, /* i=5 113.6 */
  { -0xd.016e979c75fbp-8, 0x2.42eeb4bd53f53ac77bp-60, -0x1.0054a9ed716d4p+0, 0x6.a4f0e93b528c03909ep-56, -0x6.867513d1768d4p-8, 0x7.f61f66a00ca15692c8p-64, -0x2.b2a19b9e59876p-4, 0x4.9efcb191227740e90cp-60, -0x4.ed7e9be790fb4p-8, -0x1.5956c1dc457be11ep-64, -0x1.3d35147cb93b5p-4, 0xe.cbf0a41a4c5245662p-64, -0x4.2621ec3ede32cp-8, -0x1.f8d91354acd800b6068p-60, -0xc.29cf2a40e5b38p-8, -0x3.0c42f9cf57cd1f32fdp-60, -0x3.ae2ddec277414p-8, -0x8.9c8bfe0a25d4p-8, -0x3.5e91d54f523b2p-8, -0x6.a77883a9acc4cp-8, 0xdp-8 }, /* i=6 113.3 */
  { -0xf.02335eedc94ep-8, 0x2.1eba39d3cb0b0724d5p-60, -0x1.0070ca5eda78bp+0, 0x4.4698a04828924874e8p-56, -0x7.89ee2818a4314p-8, -0x1.74c41536e8ac2c8f83p-60, -0x2.b54a166032e88p-4, 0x9.7f98a8a1503dd5dab4p-60, -0x5.b4c61705d932p-8, -0x1.7a7d3b17e346b366d88p-60, -0x1.4091f4be69f0dp-4, 0x5.e5e6853ed886151de8p-60, -0x4.d21e3d4b4b37cp-8, -0x1.d321a1f6c23dd76398p-60, -0xc.699743e9e5b78p-8, 0x1.6d23f3e8c1b2b8e291p-60, -0x4.4bb0e2a02b81p-8, -0x8.e5bd23096b05p-8, -0x3.f467e9c697c0cp-8, -0x6.f9e9ccbf379cp-8, 0xfp-8 }, /* i=7 113.1 */
  { -0x1.103347666f892p-4, 0x5.ac543b551ba7a06cep-64, -0x1.0090facbeaf37p+0, 0x3.e9d3c9cdc30da5abccp-56, -0x8.8e792a8e57c6p-8, -0x2.d75c9c599f2963dac1p-60, -0x2.b856c2237e82ep-4, 0xc.bedd5ddd6414dcfd4p-64, -0x6.7e5035e8c7bcp-8, -0x5.838a49ddecc6b7b438p-64, -0x1.44713276eea77p-4, 0x5.ae8683b5279e546554p-60, -0x5.81dd30f34aa1p-8, -0x1.99d60dd06b9805dcd28p-60, -0xc.b362f22f93a8p-8, 0x2.858dff82832e431499p-60, -0x4.eec46f68d846p-8, -0x9.3aedb6f4a101p-8, -0x4.91ee1a7de9498p-8, -0x7.5a8c64dcf1858p-8, 0x1.1p-4 }, /* i=8 112.9 */
  { -0x1.3047a02794911p-4, 0x7.5ac7a969f75addd40ep-60, -0x1.00b53fc7cb3bdp+0, -0x7.5d094975c91de18538p-56, -0x9.943be661eacc8p-8, -0x1.876d7d06081e9e5ac5p-60, -0x2.bbc8d4bfb4916p-4, 0x3.16d54e90819c9634e4p-60, -0x7.4a6f5fa5bf064p-8, -0x1.8ad6654a718ae19b4a8p-60, -0x1.48d5d4d4ce7fdp-4, 0x6.3d5871f3fa3128f8e6p-60, -0x6.35ed8879fad14p-8, -0x1.31014e960e85e3b368p-60, -0xd.07921b4bbbd58p-8, 0x2.c730a8b0ffd5c4d042p-60, -0x5.98462d19700bp-8, -0x9.9cc3f496c07b8p-8, -0x5.3867b62c3bf9cp-8, -0x7.ca68ffa1a3448p-8, 0x1.3p-4 }, /* i=9 112.7 */
  { -0x1.5060c31541da5p-4, -0x6.764f6ab0927a4a3a52p-60, -0x1.00dd9e7dc32bcp+0, 0x3.059549c81809fe094p-56, -0xa.9b5ca2bcd89c8p-8, 0x1.a4f72628dbed900aedp-60, -0x2.bfa1adbe985p-4, 0x3.aa66eb723afdfa9d4cp-60, -0x8.197801b139b9p-8, -0x2.11f0447a80867a97a5p-60, -0x1.4dc350278e557p-4, -0x3.ba7f55aaac29deb9d2p-60, -0x6.eee3a6784cd04p-8, 0x1.af0cc75f9c608cc0ccp-60, -0xd.6692e24b6478p-8, 0x2.7700bfae94b563b806p-60, -0x6.49205cad4c6dp-8, -0xa.0c00485f61528p-8, -0x5.e930acc157e9p-8, -0x8.4ab4ccac41998p-8, 0x1.5p-4 }, /* i=10 112.5 */
  { -0x1.707f33c173a99p-4, 0x4.219874125e4d6fbc06p-60, -0x1.010a1cb349899p+0, -0x5.248e08ea83600a2dep-60, -0xb.a40232985a81p-8, -0x4.3d23c1089ff87ee2p-68, -0x2.c3e2d770aef52p-4, -0x2.d8d086aabe3c43f278p-60, -0x8.ebc0d57504af8p-8, 0x1.545af62ae17c7823ap-64, -0x1.533d8a64de43dp-4, -0x4.a0c9b8699b39a42e7ap-60, -0x7.ad5a5bf99080cp-8, -0x1.820d8993ddddab4bc68p-60, -0xd.d0e27fbaa24b8p-8, -0x3.683aa4518d38a3351bp-60, -0x7.024bc1d797d34p-8, -0xa.897f64655fe28p-8, -0x6.a5c191ee2fb8p-8, -0x8.dcd5f08f398f8p-8, 0x1.7p-4 }, /* i=11 112.3 */
  { -0x1.90a376809684p-4, -0x5.2435c5b68d2c3a5204p-60, -0x1.013ac0ca53a83p+0, 0x2.1f08fb5d1d627cc792p-56, -0xc.ae5405014de58p-8, 0x1.56ab9d57cd0e1038f6p-60, -0x2.c88e08253eb86p-4, 0xf.344d66e572456c1bacp-60, -0x9.c1a328eabd9p-8, 0x1.24ab251dba978ccee6p-60, -0x1.5948e04c44b89p-4, -0x5.d8f7f5babd05ee859p-60, -0x8.71f3c24d51c88p-8, 0x2.6188a395951b5c9dcfp-60, -0xe.470e3a26adb5p-8, -0x5.231dd02d10e247984p-64, -0x7.c4d1b50ea10e4p-8, -0xb.163ca8a8aaaap-8, -0x7.6fb40a52fedfcp-8, -0x9.8268c1b2fba68p-8, 0x1.9p-4 }, /* i=12 112.1 */
  { -0x1.b0ce107d3f69p-4, -0x7.c8e4c30e758d9ec1dap-60, -0x1.016f91c3e7096p+0, -0x5.f82ca03392f122b39p-56, -0xd.ba7a35d8032p-8, 0x3.9ee9e0d1178da88f4ep-60, -0x2.cda52387715p-4, -0x2.6022c3923504cec20cp-60, -0xa.9b7b2aa0c284p-8, 0xf.4238579cf985fd2b2p-64, -0x1.5fea2b337ac83p-4, -0xa.a66c4594b7aabdaap-68, -0x9.3d5a235ec4ca8p-8, -0x1.3d91ea217005d02f42p-60, -0xe.c9b4812c5a9ep-8, 0x4.f45716eb5784996d2p-64, -0x8.91ce642495688p-8, -0xb.b354e8974117p-8, -0x8.48c7b53893098p-8, -0xa.3d45e2882bd5p-8, 0x1.bp-4 }, /* i=13 112.0 */
  { -0x1.d0ff87cc3a7a5p-4, 0x2.ddab47f13b2d58bef2p-60, -0x1.01a89742ef0f4p+0, 0x5.2b6ac9787ebb47ac52p-56, -0xe.c89d9f1a377dp-8, -0x1.f159171ce24c99081ep-60, -0x2.d32a3c226a11p-4, 0x9.463d6bfe442e0ea61cp-60, -0xb.79a839968fc1p-8, -0xf.e837e6f2394f3a8f9p-64, -0x1.6726c78718976p-4, -0x5.7ae9796351e274d83ep-60, -0xa.1040f28ac7e8p-8, 0x8.0cee182f857e4799p-68, -0xf.59862e3224898p-8, -0xc.e5e7d1e4563487d4bp-64, -0x9.6a73485883dcp-8, -0xc.620998566504p-8, -0x9.32e7a71a1b99p-8, -0xb.0f8956f3b24cp-8, 0x1.dp-4 }, /* i=14 111.8 */
  { -0x1.f1386380f8b9ap-4, -0x4.9c5b12bdd5536e8c08p-60, -0x1.01e5d98f59023p+0, -0x4.45eff70ad666bfde72p-56, -0xf.d8e7eac635b18p-8, -0x3.a68cad6f53b93ead07p-60, -0x2.d91f950e583c8p-4, 0x6.932e39fb5cc7aa1a94p-60, -0xc.5c8d396544bf8p-8, 0xe.79e7b00e522176f7fp-64, -0x1.6f049c0cbc9c2p-4, 0x2.6d9229c6d4b46e8228p-60, -0xa.eb65d81c5a4p-8, 0x3.fad11ca5521eeaa6cdp-60, -0xf.f747ee61e58e8p-8, 0x2.df6aa762e23fe3e99ap-60, -0xa.5009d9967f258p-8, -0xd.23c46e0185468p-8, -0xa.30308c0f8f3e8p-8, -0xb.fb9ab830b4d58p-8, 0x1.fp-4 }, /* i=15 111.6 */
  { -0x2.11792bc265cdep-4, -0x6.328f56db8108a55208p-60, -0x1.0227619978f0bp+0, 0x3.0656c72a57a6533c08p-56, -0x1.0eb83a566e63dp-4, 0x7.251e63a30fa9cfc33ap-60, -0x2.df87a3c8ca8dp-4, 0x9.883cde197831344cd8p-60, -0xd.4490eb31cb188p-8, -0xe.168557f111584a942p-64, -0x1.778a21f58a06dp-4, -0x3.3e6e5784ce083ad74p-60, -0xb.cf91d1ccffee8p-8, -0x3.9f7367505278b48c68p-60, -0x1.0a3d3d9f8f18ep-4, 0x7.6ebf6d08dabbd16c54p-60, -0xb.43f6878a23f98p-8, -0xd.fa1b84fd953c8p-8, -0xb.42f78c71bdc8p-8, -0xd.0436b08dd4288p-8, 0x2.1p-4 }, /* i=16 111.5 */
  { -0x2.31c269e02f15ep-4, -0xf.41b37d98e3a044cdp-60, -0x1.026d38fdba2bep+0, -0x7.794dc2a8f64a0c491p-56, -0x1.2009c51597f42p-4, -0x4.6f04486e9add4e2774p-60, -0x2.e665123aca37cp-4, -0x4.6ef3df6bbc32896facp-60, -0xe.321e4beff543p-8, -0x3.350188a881b79f85efp-60, -0x1.80be6dd1b731p-4, -0x2.eb9e3c4f6c9c44d032p-60, -0xc.bd9a6ae3b949p-8, 0xb.5ed30daa57b2eb371p-64, -0x1.1601b3d97dadp-4, -0x5.8714b39c36f458825p-60, -0xc.47bc023855868p-8, -0xe.e6d613b3f3a78p-8, -0xc.6dd211ec9b468p-8, -0xe.2c79ef890a248p-8, 0x2.3p-4 }, /* i=17 111.3 */
  { -0x2.5214a86885ad4p-4, 0x1.d1c5881b1ed7659904p-60, -0x1.02b76a089e6c7p+0, 0x3.16a7a2c3fc121f5c98p-56, -0x1.3185e7ade82bfp-4, -0xe.91cc6b81d3e03829cp-64, -0x2.edbac0ef89fdcp-4, 0x7.56c7a4942c839bc0cp-64, -0xf.25a4f889006fp-8, -0x2.a50982f79cd284c945p-60, -0x1.8aa93977f279bp-4, 0x5.479cb46e445b88c8c4p-60, -0xd.b6630ed2dfa1p-8, 0x3.fd7d4834c1b6bc186cp-60, -0x1.22d289ae01667p-4, 0x7.7eb10a5064cc1cbf28p-60, -0xd.5cfedd111bab8p-8, -0xf.ebf1c6900fc18p-8, -0xd.b39e8fba3c12p-8, -0xf.77edcf4c490f8p-8, 0x2.5p-4 }, /* i=18 111.1 */
  { -0x2.7270733e656b8p-4, 0xb.af991a19247bcb8204p-60, -0x1.0305ffbb0ef02p+0, 0xc.dca961eab7a38a93cp-60, -0x1.432f7cd08d4c2p-4, -0x7.7bc78583db43102ed6p-60, -0x2.f58bc97ec2268p-4, 0xc.56636f9e4c8336df38p-60, -0x1.01f9998845356p-4, -0x7.4f8b620d04d1bd813ep-60, -0x1.9552eef5a94d6p-4, -0x3.10dccd6c8321f3cf9cp-60, -0xe.bade79847de38p-8, 0x3.102e8058f824009ee1p-60, -0x1.30c21e667e797p-4, -0x3.a7900ba2879eab74b4p-60, -0xe.858998d793b78p-8, -0x1.10ba8d60520b3p-4, -0xf.178e76104ce98p-8, -0x1.0ea96ee54d91cp-4, 0x2.7p-4 }, /* i=19 111.0 */
  { -0x2.92d657b06b572p-4, 0x1.5f053e118fc42bc28cp-60, -0x1.035905cf0331fp+0, -0x5.f6db1637cbc28a235ap-56, -0x1.55097279eec0cp-4, -0x2.e403f26c56fcd84e3ep-60, -0x2.fddb812e2519p-4, -0x8.7eeb9bac4d70bed068p-60, -0x1.120764fdf4128p-4, -0x6.e6c6a98a5af9c13ee2p-60, -0x1.a0c4b49fc1d9cp-4, 0x7.b242896efd644807acp-60, -0xf.cc1048ded555p-8, -0x2.bc0e71168febf440f9p-60, -0x1.3fe4b09dfb36fp-4, 0x6.32a4339d39f90a0fa4p-60, -0xf.c351226350db8p-8, -0x1.24878f0be7f28p-4, -0x1.09d317f0270d7p-4, -0x1.289060bb8114fp-4, 0x2.9p-4 }, /* i=20 110.8 */
  { -0x2.b346e490467d4p-4, -0x1.30e178922672fd2c68p-60, -0x1.03b088bc812e1p+0, 0x2.119e8ed372a11f815ap-56, -0x1.6716cb603786fp-4, -0x5.6317b8423125af174p-64, -0x3.06ad7bcdb27cap-4, -0x5.16a305d39a37fef34p-64, -0x1.228bb38d02e07p-4, 0x3.a14c45243915b3026p-60, -0x1.ad087a5e24c33p-4, 0x2.450edc78a8178d70dep-60, -0x1.0eb0eb371206cp-4, -0x4.fdf30a636a9e21f00ap-60, -0x1.50508b232eab5p-4, -0x6.5c264af8217838722ep-60, -0x1.11879d6252324p-4, -0x1.3a52b1701817cp-4, -0x1.248828a6eef18p-4, -0x1.4586b822402b8p-4, 0x2.bp-4 }, /* i=21 110.6 */
  { -0x2.d3c2aa4acea2ep-4, -0xc.315ad6689e72168f98p-60, -0x1.040c95befb7d1p+0, 0x3.b1a461ad7a28d28abap-56, -0x1.795aa073c56c8p-4, -0x9.33d886ca0fdd40dbap-64, -0x3.10058ed3113fap-4, -0x3.d1aba8e7be20147b64p-60, -0x1.338eee642e1f4p-4, 0x4.5714dad9d3e98eaa44p-60, -0x1.ba29084f92f49p-4, -0x5.5370ea96251be21f8p-64, -0x1.2190468acf9bp-4, -0x5.5127ded530c6bd86a6p-60, -0x1.621e3727e7c78p-4, -0x2.a2f5d1dc888a5e5bfcp-60, -0x1.2875d1a8b50d1p-4, -0x1.524dc88cab4e2p-4, -0x1.41df6471ee683p-4, -0x1.65eadcde5c842p-4, 0x2.dp-4 }, /* i=22 110.5 */
  { -0x2.f44a3b00ccd56p-4, 0x6.702830319a344089e4p-60, -0x1.046d3adb11f07p+0, 0x3.dc6ff72fe1cb76278p-60, -0x1.8bd8227304d5ep-4, 0x2.6dd1e7608ace4a54e8p-60, -0x3.19e7d4b884948p-4, 0xd.5afa3bae7a3c1de7f4p-60, -0x1.4519eefef5627p-4, 0x3.e458efeec61765868cp-60, -0x1.c8320ee6bf6efp-4, -0x5.d64f0cc41088973abap-60, -0x1.357329f94ff5fp-4, -0x4.d738dc6af07c93d69cp-60, -0x1.7568b45ec53aap-4, 0x2.ac4be3f96228aee4e8p-60, -0x1.4128fa5eaf9aep-4, -0x1.6cb08ed6dcd48p-4, -0x1.6228bf162570bp-4, -0x1.8a2839d638ff6p-4, 0x2.fp-4 }, /* i=23 110.3 */
  { -0x3.14de2aa082788p-4, -0x9.9347563cf5335bf33cp-60, -0x1.04d286e4b9beep+0, 0x6.5ff7c4d67d8fe41c8p-60, -0x1.9e929b935d8cbp-4, 0x5.4f4bed6482d1c36ed8p-60, -0x3.2458b0a485f12p-4, -0x6.54509f2931def72818p-60, -0x1.573608902158fp-4, 0x1.9b9b9cc35fe91832fcp-60, -0x1.d73038a58a9b9p-4, 0x6.ecafec3c958a32ca6p-60, -0x1.4a6f3595f97d7p-4, -0x1.614178c211dd9d09dap-60, -0x1.8a4db7be4a7a4p-4, 0x5.6abb5e2f9ce0d3bd7cp-60, -0x1.5bce88783dcbcp-4, -0x1.89b95e198a6efp-4, -0x1.85be07b51b4b3p-4, -0x1.b2b920234581dp-4, 0x3.1p-4 }, /* i=24 110.2 */
  { -0x3.357f0efffc1dp-4, 0xd.ad54442f02da6edb6p-60, -0x1.053c8985d2bc9p+0, -0x1.8e1719e46120dc2566p-56, -0x1.b18d7140f4d4p-4, -0x6.2d362cc209fbea3cf8p-60, -0x3.2f5cd25f8fa7ap-4, -0x5.01653de3685c3b8d88p-60, -0x1.69ed12265aa66p-4, 0x5.137c99e920d196265cp-60, -0x1.e7313d9f9a867p-4, 0x1.7fa6787634d28cb1cp-64, -0x1.609bbde37b139p-4, 0x6.82e7720d713b345f08p-60, -0x1.a0edf1b99c119p-4, -0x2.087ba1015e96dd9e4p-60, -0x1.789890e99d2eap-4, -0x1.a9ae01fb55853p-4, -0x1.ad04480ebc302p-4, -0x1.e028f73f1b549p-4, 0x3.3p-4 }, /* i=25 110.0 */
  { -0x3.562d7ff83e1bap-4, 0x8.b23961f0e00283c984p-60, -0x1.05ab53452f914p+0, 0x6.da82e0193d6eae4efcp-56, -0x1.c4cc25f72be96p-4, 0x5.dbd1c3c38119a7eebap-60, -0x3.3af93a9e3219p-4, 0xe.f01d760fffc06bbb88p-60, -0x1.7d49719d7d85cp-4, 0x7.fd080542b5ea5c2b3ap-60, -0x1.f843f8f272342p-4, 0x7.2725662d988253ef1p-60, -0x1.7811f7b07f32cp-4, -0x2.9dc3fbf176ccab11f8p-60, -0x1.b96d5cdc146f9p-4, 0x1.c0d69123a6057eab9p-60, -0x1.97be5ae8fcbf8p-4, -0x1.ccdca7ef41792p-4, -0x1.d86d57ad21b5fp-4, -0x2.1316c7608e0f8p-4, 0x3.5p-4 }, /* i=26 109.8 */
  { -0x3.76ea178159bdep-4, 0x4.9138ee4b1d0db1583cp-60, -0x1.061ef58e175b7p+0, -0x7.4d9691355bf43c913ap-56, -0x1.d8525b33e883bp-4, -0x1.2dfb1df665b1796fa8p-60, -0x3.47333fa62a632p-4, 0xf.8abeaaa0d6895c0d7cp-60, -0x1.915627713e873p-4, -0x1.3657b3a8bb5a07ad22p-60, -0x2.0a788066acfp-4, 0xf.9549cd6ffaf04921dp-60, -0x1.90ed2881a903ap-4, -0x4.2ade80042745f060c8p-60, -0x1.d3f395d608c31p-4, -0x6.d43866e9b9fea9f814p-60, -0x1.b97d003824c88p-4, -0x1.f39cf0efaf657p-4, -0x2.0879ab5fc629ap-4, -0x2.4c382d81876a4p-4, 0x3.7p-4 }, /* i=27 109.7 */
  { -0x3.97b571cf7a8ccp-4, -0xb.dc6cf43cf0f67cc048p-60, -0x1.069782b847ba6p+0, 0x7.7811e9e7adfc10abcap-56, -0x1.ec23d387df2b4p-4, -0x1.aebaa64c1536545cd2p-60, -0x3.54109255de6d8p-4, 0x6.ed5d5e4683a839e18p-60, -0x1.a61edb87d792ap-4, -0x4.9c0917f45d483e4ea8p-60, -0x2.1de03e7436cfap-4, 0x5.da6682f40fc716203p-60, -0x1.ab4adc154b6cp-4, 0xb.a6b8c5ed7613b6fd6p-64, -0x1.f0ac3e2f2a27dp-4, -0x3.4b1755eb1f8ef63f28p-60, -0x1.de182202eef3ep-4, -0x2.1e512a199359ep-4, -0x2.3dba6adeadc92p-4, -0x2.8c5ccb63213cap-4, 0x3.9p-4 }, /* i=28 109.5 */
  { -0x3.b8902d70fc4bcp-4, -0xf.c95f78e45444028cf4p-60, -0x1.07150e107ec92p+0, -0x6.2f7adccb6183b5ep-68, -0x2.004474c645cep-4, -0x1.2ecbc9cead273601dcp-60, -0x3.619743966155ap-4, -0xa.d1321efd05edf21884p-60, -0x1.bbafeb0dd91c5p-4, 0x5.47dd35303e8f864a62p-60, -0x2.328e0eea5ab92p-4, 0x6.840939050bbcaee1bcp-60, -0x1.c74b1f9f53817p-4, 0x1.03c112c3b1db775332p-60, -0x2.0fc76afebfa5p-4, 0x9.ad2bee3dcc7aff22f8p-60, -0x2.05dab5423ec4p-4, -0x2.4d67b21350e38p-4, -0x2.78d3d825cc89cp-4, -0x2.d47249519f738p-4, 0x3.bp-4 }, /* i=29 109.3 */
  { -0x3.d97aeb6d9b1ep-4, -0xd.7d17d85154db0b236p-60, -0x1.0797abe195398p+0, -0x4.d76f40edee54d158fcp-56, -0x2.14b84a568837ap-4, 0xa.deb2e2360b9aa7eec8p-60, -0x3.6fcdca40f6a5p-4, 0x5.ec556d3f78065cfa8cp-60, -0x1.d216777eea9ecp-4, 0x2.d93ffb2e307de20a8cp-60, -0x2.48965e8464d2ap-4, -0xb.c15c1c05fcbaf1982cp-60, -0x1.e510c3707c1b8p-4, -0x7.456d295d5ea2087ffep-60, -0x2.317a213c721ccp-4, 0x6.18082801d94df180a8p-60, -0x2.3117ea064eae4p-4, -0x2.815c9225baa02p-4, -0x2.ba8016235f49ep-4, -0x3.25890384d8b9cp-4, 0x3.dp-4 }, /* i=30 109.2 */
  { -0x3.fa764f66d066p-4, -0x9.8ecb00191242f74af4p-60, -0x1.081f717e3165ep+0, -0x2.5ac18e65bde23e013ep-56, -0x2.298387aad08bep-4, 0x1.de6c5d9cef5ee883dcp-60, -0x3.7ebb0981f4e28p-4, 0xd.ee8d56acb6c5a425dcp-60, -0x1.e96076fa48a37p-4, -0x6.42508937471f6098b2p-60, -0x2.600f4dbc8791ap-4, -0xa.dfc1e746dc170eb274p-60, -0x2.04c1a3d374744p-4, 0x9.ec97b77fd46723a24p-64, -0x2.55fee176a75f4p-4, -0x3.abb47e8629cc5a9a04p-60, -0x2.602c319176876p-4, -0x2.babb53154de38p-4, -0x3.03925ec1906e8p-4, -0x3.80d981e823b4p-4, 0x3.fp-4 }, /* i=31 109.0 */
  { -0x4.1b82ffb97015cp-4, -0xe.dadddd6e6658663c9p-60, -0x1.08ac754b1cf11p+0, 0x4.9d7d625213f2e205e6p-56, -0x2.3eaa8ade74246p-4, -0xc.3e8db0e4005fe743d8p-60, -0x3.8e6657c3ebcbep-4, -0xe.834dfe9fcb8027a8fcp-60, -0x2.019cc60530faap-4, -0x8.415bbab249e7663ccp-64, -0x2.7910d738cff18p-4, -0xa.1e0b8236fc3a3a157p-60, -0x2.2686fa05a59c8p-4, -0x4.1c318dea901094fedp-60, -0x2.7d9644fec27b4p-4, 0x7.1116f548ca6259b82p-60, -0x2.937e67c979056p-4, -0x2.fa2117266fc22p-4, -0x3.54fabb3851ca6p-4, -0x3.e7cadf03c2002p-4, 0x4.1p-4 }, /* i=32 108.8 */
  { -0x4.3ca1a5a09b654p-4, 0x1.6b73b247c58dd21408p-60, -0x1.093ececa4756p+0, 0x5.7324aeffa310273ddep-56, -0x2.5431df7f99ddp-4, -0x9.82afce910b0870511p-60, -0x3.9ed7862f12f2cp-4, 0x8.490c08239577879cap-60, -0x2.1adb3ae13beaep-4, -0xa.60e5cc0babd81f5abcp-60, -0x2.93b4fa396efcap-4, -0xb.f7dd93ba3f12dc301cp-60, -0x2.4a8db6503b9f8p-4, 0x1.0f00a444a01e35249cp-60, -0x2.a887aee996584p-4, 0x2.f0537bd7860435c334p-60, -0x2.cb81253cfc14ap-4, -0x3.403f045366b5ep-4, -0x3.afca5508b9da2p-4, -0x4.5bfa49e9b3f08p-4, 0x4.3p-4 }, /* i=33 108.6 */
  { -0x4.5dd2ed5a1f3e4p-4, 0x1.a2fb3e84f654a73995p-56, -0x1.09d696a680a83p+0, 0x2.47c0929f0472e3c1ep-56, -0x2.6a1e4187baa92p-4, 0xa.e0e7f5e23873ca784p-60, -0x3.b016e8ca46396p-4, -0x3.384ca8faf9401f12dcp-60, -0x2.352cbaa0e9a52p-4, 0x3.9d37c9105b583df7acp-60, -0x2.b017e97ca701ep-4, 0xa.34ebe7a665a7a0b2cp-64, -0x2.7106e469a87fp-4, 0x4.a35d6e7decb4363dp-60, -0x2.d722139cbd4p-4, -0x5.4fec35f6f6a2ae3f4p-64, -0x3.08b43fd1f14e2p-4, -0x3.8ddd0ba261114p-4, -0x4.15387a248137cp-4, -0x4.df43d7ad25f5p-4, 0x4.5p-4 }, /* i=34 108.5 */
  { -0x4.7f17864c5624cp-4, 0xa.588412338764d188p-68, -0x1.0a73e6bff8b8p+0, 0x3.73c3cc4f16550cb484p-56, -0x2.8074a086f1088p-4, -0xc.6f32d68571adb308cp-64, -0x3.c22d5f4c3ad6ap-4, 0x7.1a154b9017c680211cp-60, -0x2.50a350297ed8ep-4, 0x3.40fffb9ed531128348p-60, -0x2.ce583f1b7d038p-4, 0x9.16bee573a017148034p-60, -0x2.9a281b70bd686p-4, 0xa.93b887fbe3cb66f344p-60, -0x3.09bcda0a8226p-4, -0x2.50a6e72de63fdebd7cp-60, -0x3.4ba6812e3417ap-4, -0x3.e3dd1ccbf3774p-4, -0x4.86a873b7747f4p-4, -0x5.73cce361b686p-4, 0x4.7p-4 }, /* i=35 108.3 */
  { -0x4.a070232da6fap-4, -0x2.340dca0c21a412c9fp-60, -0x1.0b16da398fd6cp+0, -0x6.5297fa72c22cd62bap-60, -0x2.973a23066887ep-4, 0x5.15269d5ed7ff29f15p-60, -0x3.d5245ead25ae2p-4, 0x4.34b0d7e81e89913dccp-60, -0x2.6d524556bd1cep-4, -0x9.ea4c62cfa22823fca8p-60, -0x2.ee9735f33027ep-4, -0xc.f88672bf824ced282cp-60, -0x2.c62bfb018f0dcp-4, -0x1.b5b2819713a3ed375cp-60, -0x3.40b8da38234b4p-4, 0xd.caa15b6486c8ccabd8p-60, -0x3.94f79b147c6dep-4, -0x4.433ed8112596p-4, -0x5.05b054f975cd4p-4, -0x6.1c104851e653cp-4, 0x4.9p-4 }, /* i=36 108.1 */
  { -0x4.c1dd7a2dbbae4p-4, -0x1.0b08737cd3bb19be828p-56, -0x1.0bbf8d8707978p+0, -0x5.7ea6a165241f7047fep-56, -0x2.ae742a26a4e8ap-4, 0x4.154069c9b9d97a9cf8p-60, -0x3.e905fb8ac84fap-4, 0xf.6b559b50cd20a39d74p-60, -0x2.8b4e3e7b2a12cp-4, -0xe.60dd8715de8f9d38a8p-60, -0x3.10f8e942b983cp-4, -0x4.c8b36c7550bd5dfd94p-60, -0x2.f552b71aa94dep-4, -0xa.6e66975ea2eb72764cp-60, -0x3.7c817d3877d5cp-4, 0x9.6c4e2898c21005f8e8p-60, -0x3.e55a634d58bb2p-4, -0x4.ad23d35140274p-4, -0x5.9420ee94d4658p-4, -0x6.daecd187ed12p-4, 0x4.bp-4 }, /* i=37 107.9 */
  { -0x4.e36045208cda4p-4, -0x6.fd9112162b0f9265e8p-60, -0x1.0c6e1e7c2324dp+0, -0x5.410e2454c525f3b8cap-56, -0x2.c628557eca058p-4, 0x6.05694bf57525186e1cp-60, -0x3.fddcf572cabf8p-4, -0xc.af1e2211e4fc5d7474p-60, -0x2.aaad587e8c45cp-4, 0x2.269e36b3f502fe846p-60, -0x3.35a49b3994068p-4, 0xc.5cd7f882f95efd295p-60, -0x3.27e2b4c7cee08p-4, -0x7.cc1c62f46174401fa4p-60, -0x3.bd8e036b04dbcp-4, 0x7.4368b0d9779e2624f8p-60, -0x4.3d97625041fd4p-4, -0x5.22d48b3db64p-4, -0x6.340f1bc91b7dp-4, -0x7.b3b64c8e5a86p-4, 0x4.dp-4 }, /* i=38 107.7 */
  { -0x5.04f941ab613b8p-4, -0x1.71db7c5df95849dda4p-56, -0x1.0d22ac5cc81dep+0, -0x2.3aa9d121ca63604f54p-56, -0x2.de5c87428b502p-4, 0xe.8505b354394a92ea18p-60, -0x4.13b4c329743d8p-4, -0xb.c2f3d1c0dfe84df798p-60, -0x2.cb8749e424328p-4, -0xd.32a02d64bc99a758f4p-60, -0x3.5cc5034e33d68p-4, 0x3.d69498399ff39b02a8p-60, -0x3.5e2939d00d50cp-4, -0x5.df01d900b72ef3b0a8p-60, -0x4.0462f6919d25cp-4, -0x1.a6e26789701839dea88p-56, -0x4.9e8fc1c83075cp-4, -0x5.a5c61e161a24p-4, -0x6.e7dea8d7c3e28p-4, -0x8.aa49c43bc20f8p-4, 0x4.fp-4 }, /* i=39 107.5 */
  { -0x5.26a93173e2498p-4, -0x1.c1d2f6ce8524972b8fp-56, -0x1.0ddd57ee32684p+0, 0x5.60c4175a3bbca84a36p-56, -0x2.f716e8b4f55c2p-4, -0x2.b8f5fdc59ec82bc12cp-60, -0x4.2a99a00544604p-4, 0x1.857d286f49098ea79e8p-56, -0x2.edf5870b14bdap-4, -0x3.783022d722664bf724p-60, -0x3.8688a54f659cep-4, -0x6.9acd36e498acf809bcp-60, -0x3.987b31f8f79fcp-4, -0x5.b2fa1b5b86c057f528p-60, -0x4.5193ce31e24fp-4, 0xf.94babdf63251b30788p-60, -0x5.0940aa500605p-4, -0x6.37a0f2c69d964p-4, -0x7.b24f1edcadep-4, -0x9.c32584684d0b8p-4, 0x5.1p-4 }, /* i=40 107.4 */
  { -0x5.4870da51697c8p-4, 0xa.f3f00c5933244cd12p-60, -0x1.0e9e43894f02dp+0, 0x1.fe18b5e7a90293b0d6p-56, -0x3.105deeeedefep-4, -0x8.64d0e1aaa6b57ae928p-60, -0x4.42989a7aa9258p-4, -0x1.c89e9060121412dc588p-56, -0x3.12136a059013ep-4, 0xf.1d4109a53db39c4adp-60, -0x3.b3223244de9aap-4, -0xe.0534f4060134750b4p-64, -0x3.d7360cd25c95ep-4, -0xa.b0e65b0c17438d2e78p-60, -0x4.a5c4cdc2890acp-4, 0x4.d9e757ba9613229c5p-60, -0x5.7ec7226f6c508p-4, -0x6.da4875c98504p-4, -0x8.968ad074a0a28p-4, -0xb.0385a9ced20bp-4, 0x5.3p-4 }, /* i=41 107.2 */
  { -0x5.6a510680ac7f4p-4, 0x1.35e8dbcb0816dda8298p-56, -0x1.0f65932e53afap+0, -0x9.0802b635672891fbap-60, -0x3.2a3860007bf4ep-4, -0x5.15265b489ff58876f4p-60, -0x4.5bbfa3f62053p-4, 0x1.4ece06c9f4263144e38p-56, -0x3.37fe5e6ee258ep-4, 0xc.b5b4a24a8b914a7aacp-60, -0x3.e2c8f45945354p-4, 0x3.d437351c678e2cae38p-60, -0x4.1ac0b76f21618p-4, -0x1.9ef0c5ca7010cb0246p-56, -0x5.01ad255d7975p-4, -0x1.6682f3d265ed42cea38p-56, -0x6.0064840deb804p-4, -0x7.8fe41d0280c98p-4, -0x9.9838951c8f3b8p-4, -0xc.71862b9971a28p-4, 0x5.5p-4 }, /* i=42 107.0 */
  { -0x5.8c4a84d9f15fp-4, 0x1.e8ea4df76e5db565p-60, -0x1.10336c99ab3f9p+0, -0x5.7a7f42452c8142d712p-56, -0x3.44ad58764b17p-4, -0x1.9b92b9eb996185b54p-60, -0x4.761da2268694cp-4, -0xa.bb3c172a1e749df9c8p-60, -0x3.5fd611b58fd5p-4, 0x7.facdf930857397df0cp-60, -0x4.15b94733b4fbcp-4, 0xb.30bd717825f8349478p-60, -0x4.638cb5f3d70e8p-4, 0x1.4cd2106ddd95676dd8p-60, -0x5.66195f209f848p-4, -0x1.740d074d19a658bf91p-56, -0x6.8f83a15a279a8p-4, -0x8.5ae9ed3aed6a8p-4, -0xa.bb90b7c5675p-4, -0xe.144b8020c43ep-4, 0x5.7p-4 }, /* i=43 106.8 */
  { -0x5.ae5e2909f6e54p-4, 0xd.8ee857e2bcf3b2675p-60, -0x1.1107f75a50748p+0, 0x4.c3dbfa6999c755e966p-56, -0x3.5fc4514a6ff64p-4, 0x7.e271422cf24cfbaafcp-60, -0x4.91c281dd3afcp-4, 0x1.0a443e6aeb49f630adp-56, -0x3.89bca85dae836p-4, -0x9.245287cd529efd5464p-60, -0x4.4c351e592b874p-4, 0x1.ad06c394272f3dd72ap-56, -0x4.b217617d7b168p-4, 0x1.4124567e5e8e4591fb8p-56, -0x5.d3ee251a258bp-4, -0xe.d50b7ca2cc58bcd51p-60, -0x7.2dbeb6a634bdcp-4, -0x9.3e2ac749498bp-4, -0xc.0575b8a3c2d5p-4, -0xf.f4333e9eb369p-4, 0x5.9p-4 }, /* i=44 106.6 */
  { -0x5.d08ccbcdc0adp-4, -0xc.3f474f7bd97c2dd928p-60, -0x1.11e35ce9b3cedp+0, 0x8.e374d3b0ae4b072b2p-60, -0x3.7b85264c76c1ap-4, 0x5.7439e8004d4f7fb608p-60, -0x4.aebf4bae216a8p-4, -0x1.6ade5cf49a71524d2f8p-56, -0x3.b5d6f8d0b339ep-4, 0x5.a961a2590568f798c8p-60, -0x4.86849b688328cp-4, -0x1.c8a0e9272c2111d1a48p-56, -0x5.06eb4f968e428p-4, 0x1.44c144de47dd20514e8p-56, -0x6.4c2b6d9ebd288p-4, 0x1.7c14fa60dd0b80c1458p-56, -0x7.dce64c1eaa5b8p-4, -0xa.3ce0d0b38a358p-4, -0xd.7b91a3f455048p-4, -0x1.21b0e6f1060c5p+0, 0x5.bp-4 }, /* i=45 106.4 */
  { -0x5.f2d74b317a3b4p-4, 0x4.6d49b267f40f1e9a58p-60, -0x1.12c5c8c54b3d2p+0, 0x3.1eb5c819784a14692p-56, -0x3.97f81d0a8bfd2p-4, -0x7.d58da1523c69dfd26cp-60, -0x4.cd263a7e8496p-4, -0xb.0ca9207e2f8ce63e28p-60, -0x3.e44ccc636d8f6p-4, 0x5.5590632799fe7b691p-60, -0x4.c4f6b646f72a4p-4, -0x1.9ccd99fbbd8b7dabe2p-56, -0x5.62a1e94304dep-4, 0x1.a8df74ab4852f127cp-56, -0x6.cff01e491510cp-4, 0x4.0837df25046c4ce838p-60, -0x8.9f0930b815d78p-4, -0xb.5ac05c4c79bd8p-4, -0xf.2478e8be5beep-4, -0x1.494679c1b0585p+0, 0x5.dp-4 }, /* i=46 106.2 */
  { -0x6.153e8ad2a868cp-4, 0x6.8fbd13590e196fbd6p-60, -0x1.13af6889ed89ep+0, 0x6.bb52720a3f089b8f5p-56, -0x3.b525ec485c2a8p-4, 0xe.b8962707447960bcacp-60, -0x4.ed0ad43775714p-4, -0x1.3a07bc876f6e62fe8f8p-56, -0x4.1549275145fccp-4, -0xc.6312689175b8fd1bdp-60, -0x5.07e1f9a133d38p-4, 0x1.b43eaaac308d989ad8p-56, -0x5.c5e538a23356p-4, -0x1.1cddb1d4d55dedc28b8p-56, -0x7.607e38b37a024p-4, 0xf.4b42cf83d04f0faf7p-60, -0x9.767dbf6ff25bp-4, -0xc.9c0bcb92192e8p-4, -0x1.107d3cd657638p+0, -0x1.76dd735689583p+0, 0x5.fp-4 }, /* i=47 106.0 */
  { -0x6.37c37425e4e1p-4, -0x8.53ee398078eac375f8p-60, -0x1.14a06c111e9abp+0, 0x1.576036fe881f2758ap-60, -0x3.d317c4111ddfep-4, -0x7.9ee9fca5c9c33c9adp-60, -0x5.0e8204d6a29dp-4, -0x3.63895582e73339cb38p-60, -0x4.48fa9885734dp-4, 0x1.7480b66193a4dcb5a9p-56, -0x5.4fa5567f8571cp-4, 0x1.cbeab9ff4470960633p-56, -0x6.3171f5479150cp-4, 0x1.eb8eef8ae36724d7cc8p-56, -0x7.ff3fa756c991p-4, 0x1.4a380399c690a8f4dep-56, -0xa.65eca95307fc8p-4, -0xe.05aafa41a849p-4, -0x1.32e8fde6df549p+0, -0x1.ab76969d8a808p+0, 0x6.1p-4 }, /* i=48 105.8 */
  { -0x6.5a66f6c0644e4p-4, 0x1.b53a9536565c3cfe55p-56, -0x1.159905907522cp+0, -0x1.2ab039a4dd5bef6922p-56, -0x3.f1d75673bfe0ap-4, 0xd.9feedc27251a377f98p-60, -0x5.31a23c1fdb2fp-4, -0xe.cb8ddf0402274f97p-64, -0x4.7f9392284f734p-4, -0x2.8043021e6e18e8481p-60, -0x5.9ca9121509eb4p-4, -0xf.8e5b920f3d973c3e08p-60, -0x6.a619d8c02b9fp-4, 0x1.2d7c958df8249bddbap-56, -0x8.adcbc4f82b35p-4, -0x3.6b9d43ecf246fe573p-56, -0xb.705d89fb56818p-4, -0xf.9d46e3b5653bp-4, -0x1.5a31b012142dap+0, -0x1.e84196012190fp+0, 0x6.3p-4 }, /* i=49 105.6 */
  { -0x6.7d2a08a58d1ecp-4, 0x1.71ef7327550cd438838p-56, -0x1.169969bb465b6p+0, -0x6.2868ede5b9c003a64ep-56, -0x4.116ee0f9da134p-4, -0x4.261e15125a91db9338p-60, -0x5.56838e39a5c6cp-4, -0xe.a453c41af73cf0d92p-60, -0x4.b94acc08d8834p-4, 0x3.15f21a5ac38b4909f8p-60, -0x5.ef5fd16c78ba4p-4, -0x3.26243b89f9eb5bf278p-60, -0x7.24c6464cc360cp-4, 0xe.5c650dbe668ee2ff1p-60, -0x9.6dedbdd3797fp-4, 0x2.02a1dad8accc2d9eafp-56, -0xc.994599d08fe3p-4, -0x1.1696a567a37ccp+0, -0x1.871aa3378a081p+0, -0x2.2ea66a3b041b8p+0, 0x6.5p-4 }, /* i=50 105.3 */
  { -0x6.a00da698eaaap-4, 0x4.de9e53952b1db5546p-60, -0x1.17a1cfe6c8b8dp+0, -0x1.a5eccfa2a616d84b2p-56, -0x4.31e936ebda1ccp-4, -0x5.10bce18cd1d4e48cb8p-60, -0x5.7d3fd788addb4p-4, -0x9.993a8a508577eab8p-60, -0x4.f65bb1217dba4p-4, -0xa.bdfeb076b3be78dcb8p-60, -0x6.4847c72448b7cp-4, 0x3.0cbe45f1ed98df0abp-60, -0x7.ae7b52b3d5c24p-4, 0x8.6d5d81e1c242c65198p-60, -0xa.41abef647e8d8p-4, -0x3.0d3feb175ea5951829p-56, -0xd.e498e1cec387p-4, -0x1.371a8b7be0a7cp+0, -0x1.ba88e799161fp+0, -0x2.8050aca9f9794p+0, 0x6.7p-4 }, /* i=51 105.1 */
  { -0x6.c312d474cee18p-4, -0x1.beb0350ac9f2184874p-56, -0x1.18b27230e36f8p+0, 0x3.2bb1d164c3df158772p-56, -0x4.5351cc76fa33p-4, 0x5.41dee430f0ff3148p-68, -0x5.a5f2e4287257p-4, -0x1.c076f617343e58b19a8p-56, -0x5.3706d9a479b9p-4, 0x4.d6f8e3b1230355f488p-60, -0x6.a7ec081372f28p-4, -0x1.9af47f4bb1b6bd88bf8p-56, -0x8.445b3b4191c38p-4, -0x2.bae3e05f1482056b0bp-56, -0xb.2b5072882d748p-4, 0x3.5f2bbf548f49ccb2c9p-56, -0xf.56de68970264p-4, -0x1.5becc0fcabeecp+0, -0x1.f589b0cede34dp+0, -0x2.df3d6d22f185ap+0, 0x6.9p-4 }, /* i=52 104.9 */
  { -0x6.e63a9d85fbce8p-4, -0xf.28344ce2b74e3c76bp-60, -0x1.19cb8da9e60cbp+0, -0x6.58873022d31ebfc1e8p-56, -0x4.75b4c2cbeafa4p-4, -0x1.7b46723fb3e861e3e68p-56, -0x5.d0ba9b6bce64cp-4, 0x1.8f9bef9d13bb6dd332p-56, -0x5.7b929330329a4p-4, -0x9.8716e98891aab3e85p-60, -0x7.0ee60c744e6ap-4, -0x1.27fee63fc381cc508d8p-56, -0x8.e7aa5da7a7b7p-4, -0x1.a03a5f9d2c306857ffp-56, -0xc.2d7303e798da8p-4, -0x1.60a9cdb5865d7202c2p-56, -0x1.0f547e8e66b4ap+0, -0x1.85b0bdd4565bap+0, -0x2.395a09a9134eap+0, -0x3.4dcc193d4e186p+0, 0x6.bp-4 }, /* i=53 104.7 */
  { -0x7.098614ecb5e14p-4, 0x8.a709737c8d31240f5p-60, -0x1.1aed628169723p+0, 0x1.70ca5a4bc35244accap-56, -0x4.991ef54fc1b14p-4, 0x7.94c97bc7258dae0d3p-60, -0x5.fdb72fdde53p-4, 0x1.e4a37cdadfb5ae1c8dp-56, -0x5.c44b7908a0818p-4, -0x1.7a0ff5e5d31279f3b98p-56, -0x7.7ddf64196b524p-4, 0x1.20ef92ad67d6c953698p-56, -0x9.99d3c589412cp-4, -0x3.341b990b368d1e64b6p-56, -0xd.4b048885efcep-4, 0x1.5eca46891bf9ba64b6p-56, -0x1.2c5cdbf4623e8p+0, -0x1.b524e819b4463p+0, -0x2.8770244685c46p+0, -0x3.ced331311974p+0, 0x6.dp-4 }, /* i=54 104.5 */
  { -0x7.2cf65603a8aep-4, -0x1.d9730350eb04931fcfp-56, -0x1.1c183436a246fp+0, 0xb.85360fe51c0909004p-60, -0x4.bd9e07fbbdf8p-4, -0xd.82e8346b83375743b8p-60, -0x6.2d0b544bf3914p-4, -0x1.df628c9dfc437adaaep-56, -0x6.11851e7908d88p-4, 0xe.05d0d983e70c2c5d1p-60, -0x7.f593a53f72a68p-4, -0x1.fa48707713e2b88ed98p-56, -0xa.5c6e683266ea8p-4, -0x1.99b32d390d69bf3392p-56, -0xe.875c7840eeaa8p-4, 0x5.39fa8ae4931eba1d2p-60, -0x1.4cf4fe1210b81p+0, -0x1.eb27782029e0ep+0, -0x2.e186a3d0c99d2p+0, -0x4.65b9b6d430f1p+0, 0x6.fp-4 }, /* i=55 104.2 */
  { -0x7.508c84cd112b4p-4, 0x1.c74e53df13c5e4059fp-56, -0x1.1d4c49cc747fdp+0, -0x9.942b1b1ecc27ab3c4p-60, -0x4.e340770be35a8p-4, 0x2.e99063efc8b8684278p-60, -0x6.5edc7672e3f7p-4, -0x1.b33547bf36009734f38p-56, -0x6.639acdd32cfp-4, -0x6.e7dc420c8d0390a298p-60, -0x8.76d29ecec0458p-4, 0x1.a9c28f91da5b9f3ep-56, -0xb.31432b65923c8p-4, 0x3.60807d2760da05d572p-56, -0xf.e648868e8c8c8p-4, 0x2.f5de4259e217fb536dp-56, -0x1.719bcdc04db07p+0, -0x2.28bc55af74b54p+0, -0x3.49aa542f6228cp+0, -0x5.1696811eca154p+0, 0x7.1p-4 }, /* i=56 104.0 */
  { -0x7.7449ce66aaf9cp-4, -0xe.835d60abc43f4d49cp-60, -0x1.1e89ee01a019cp+0, -0x2.ad7149335d0a1fe688p-56, -0x5.0a15a82041c18p-4, -0xb.1a65db2e7abbf47a28p-60, -0x6.935300009d9dp-4, 0xb.5a107375bed94ed938p-60, -0x6.baf05ee071f6p-4, 0x9.4001a7c1df05fecd2p-60, -0x9.0282d8593b928p-4, 0x1.6f456e41470c04ea09p-56, -0xc.1a53d9754a838p-4, -0x1.cdc30a9ede741b17a1p-56, -0x1.16c1ef356afdap+0, 0x7.0f9137a13255a64962p-56, -0x1.9ae4012b0f4fap+0, -0x2.6f14246d06f46p+0, -0x3.c24ae0eba3b18p+0, -0x5.e656edd602c4cp+0, 0x7.3p-4 }, /* i=57 103.8 */
  { -0x7.982f6984e9e3cp-4, -0xe.153898d98a350f5bp-64, -0x1.1fd16f8d68a0ep+0, -0x7.400f5083a7eb66c82ap-56, -0x5.322dfcf92f124p-4, -0x1.9c448e0268e10d16a68p-56, -0x6.ca9a9eb15e9c8p-4, -0x1.a54ad8c37e4b1e7cf4p-56, -0x7.17f3280903424p-4, -0x1.c59d4873fd287141e08p-56, -0x9.99a46be06fe98p-4, 0x2.b1f952322f1770b0ccp-56, -0xd.19e32b3ab069p-4, -0x3.431be98b1f2b551dfdp-56, -0x1.31dd4032c36bfp+0, 0x7.7e0ef36adbec1b116p-56, -0x1.c9778b142bf9bp+0, -0x2.bf94c238270a8p+0, -0x4.4e4f41232d628p+0, -0x6.daeed0c241848p+0, 0x7.5p-4 }, /* i=58 103.5 */
  { -0x7.bc3e96f615894p-4, -0x1.7b9360c3b82f34f4ef8p-56, -0x1.2123216123cd1p+0, 0x7.a28968ed6428a2990ep-56, -0x5.5b9ae7fbb865cp-4, 0x1.a1592351c7ba0f6bf28p-56, -0x7.04e2946d8cc1p-4, -0x1.7c68277d9d76bc245p-60, -0x7.7b1b0df77befp-4, 0x1.928b66e9d6381debde8p-56, -0xa.3d54478329628p-4, 0x2.7a93d95c425006dcb6p-56, -0xe.327e182028f2p-4, -0x1.e2d3c7d97a496f8a57p-56, -0x1.5011338438e55p+0, -0x2.30dc7f5b2a8b6cc104p-56, -0x1.fe1bb258a5094p+0, -0x3.1be38a2d455ecp+0, -0x4.f12ebcffa822p+0, -0x7.fb9403108ebbcp+0, 0x7.7p-4 }, /* i=59 103.3 */
  { -0x7.e078a22deb88p-4, 0x1.0803b3e82eddbc8b01p-56, -0x1.227f5aef2780cp+0, 0x4.16f507a6d4ecd6f6a4p-56, -0x5.866f02b14b018p-4, 0xb.2c4fb9cead1e6bb1ap-60, -0x7.425e106cf1d28p-4, 0xa.98b003e850f9d277fp-60, -0x7.e4ebb618da388p-4, 0x1.846e9d5cdd461044178p-56, -0xa.eecfe7aa1511p-4, 0x3.35dfb4d35b68d24dc1p-56, -0xf.6706a5b28ca9p-4, -0x3.c92520700083d74b43p-56, -0x1.71c5d0576da1ap+0, 0x3.02d3cfc1e3fbd320ap-60, -0x2.39b5f5077ee7ep+0, -0x3.85f1c316dc9ap+0, -0x5.af0fac9f3fb6cp+0, -0x9.51089f486fea8p+0, 0x7.9p-4 }, /* i=60 103.0 */
  { -0x8.04dee1da7c3b8p-4, 0x1.3064f6236756b8a9dbp-56, -0x1.23e678779d056p+0, 0x5.63783e8d83e13c1a4ep-56, -0x5.b2be267c215cp-4, -0xd.5aec252040d7b824f8p-60, -0x7.8344928a135d4p-4, 0x9.a8d69ccb7514ea703p-60, -0x8.55f5e10cb5b38p-4, 0x2.e25f882e83e8b2d73ap-56, -0xb.af799d3d205dp-4, -0x3.3baed139eb1285d95dp-56, -0x1.0bac07b569c76p+0, -0x3.b6ab4c0fa87ee9918p-60, -0x1.97729dbdb45cbp+0, 0x4.04c556d53dbe858c74p-56, -0x2.7d51e26e3dee2p+0, -0x4.000bb9682749p+0, -0x6.8ced5c79a4c7cp+0, -0xa.e5f7d2be77a48p+0, 0x7.bp-4 }, /* i=61 102.8 */
  { -0x8.2972b88307eap-4, 0x3.92a14219ff6137e4bp-56, -0x1.2558db5bceab9p+0, -0x3.ebc3d57b5eaf241456p-56, -0x5.e09d87c180cc8p-4, 0x5.e6ba8cd662253b88d8p-60, -0x7.c7d25a1c30238p-4, -0xb.9599fc71f13dcbed38p-60, -0x8.ced8f2ee29998p-4, -0x1.368aacdcf4b7a212f4p-56, -0xc.80dd76f16448p-4, -0x1.27ae29341b45dc6c97p-56, -0x1.2315f8c4daab1p+0, -0x6.449a99a24001933212p-56, -0x1.c1a13a3aacac7p+0, 0x4.bb596392381b38327ep-56, -0x2.ca2821967b29ap+0, -0x4.8ceb24f521e3cp+0, -0x7.90c6e51e9e1ccp+0, -0xc.c76a4cd68eca8p+0, 0x7.dp-4 }, /* i=62 102.6 */
  { -0x8.4e359531b5db8p-4, -0xd.ea1e9a3200cf0d924p-60, -0x1.26d6ea789648cp+0, -0x1.16f6b7d6efb71c02aep-56, -0x6.1023d3d371198p-4, 0xd.5d66587f875804e3fp-60, -0x8.1048e1e18136p-4, -0x2.0b38cc3907496ecd7cp-56, -0x9.5044b058fe218p-4, 0x2.d62ee4e9b9c8acf165p-56, -0xd.64b6e8f918b28p-4, -0x6.53c844f55d6b43487p-60, -0x1.3cf193ab14b55p+0, 0x5.ca89f73a03669acbf2p-56, -0x1.f0f06ecd5ff36p+0, -0x4.5ee92fe852d688ebb8p-56, -0x3.21a6e653f36c8p+0, -0x5.2fcda36496468p+0, -0x8.c1d9425108cbp+0, -0xf.0558cc82b0408p+0, 0x7.fp-4 }, /* i=63 102.3 */
  { 0x1.7af961ef39d23p+0, 0x1.7a0a884eb5ed49d27ep-56, -0x2.6d50b66b084e4p-4, -0xb.90201654085ef46fe4p-60, 0xb.3e1c32d40fabp-8, 0x1.0579cca1ffa5f743d5p-60, -0x4.5ee7d4a59065p-8, -0x8.812f079071ecc24bc8p-64, 0x1.f4b0a2c3879e3p-8, -0x7.b0af306ba75a443b4ep-64, -0xf.34930b8c1f0d8p-12, -0xd.8cb19026383c5bccfp-68, 0x7.cc7c604badf28p-12, -0x1.6fcfd5ece8a00dafd9p-64, -0x4.27dce47e41678p-12, 0x1.4caa3ad84917ab084ap-64, 0x2.473d20538eaaep-12, -0x1.46d8ef5a3db71p-12, 0xb.a65f48ebbe6d8p-16, -0x6.bd31ee5e96e7p-16, 0x8.1p-4 }, /* i=64 121.6 */
  { 0x1.7aabe4ae052a9p+0, 0x2.d9b9e43acfb16f4ep-56, -0x2.6a8472a78c1cep-4, 0x6.f2d0503fb0bc3370bp-60, 0xb.241168d4df21p-8, -0x2.f9cf91cf5205cee702p-60, -0x4.4f680592aaf64p-8, 0xa.bd10f61500464779a8p-64, 0x1.eb4cbd935ae08p-8, 0x1.e4456a1ff5ec6f994ap-64, -0xe.d85658cfc0edp-12, -0x3.ddb50ff0b16517928bp-64, 0x7.934a2607cd244p-12, -0x1.00e44c72c91527c7148p-64, -0x4.041e3a1bbd0d4p-12, -0xa.3f1ca68af9311ae9p-72, 0x2.30c2ca3d644bap-12, -0x1.38a4049141523p-12, 0xb.1619da740868p-16, -0x6.614defb93a4f8p-16, 0x8.3p-4 }, /* i=65 121.7 */
  { 0x1.7a5ec08d99337p+0, -0x8.08d97c04118cc8ecp-60, -0x2.67bea6095f8b6p-4, -0xc.fd2ece58f95128a3d8p-60, 0xb.0a62be1c7ed28p-8, 0x3.dc79537703aae54b91p-60, -0x4.40327148a50ecp-8, -0x1.020d24b7cb23d11f88p-64, 0x1.e221aa41fa701p-8, 0x2.6cc9d4ac5bc4e0d664p-64, -0xe.7ebc6d25f867p-12, 0x3.ed5dc82ecac6578bb8p-64, 0x7.5c02a8657e28p-12, -0x2.d7dfd22014087aae78p-68, -0x3.e1bf58b19bb2ap-12, 0x8.c8df31c03299a9be98p-68, 0x2.1b41f2d0c3b16p-12, -0x1.2b1e990fc00dfp-12, 0xa.8d7ec4e4ae07p-16, -0x6.0abe1256cf12cp-16, 0x8.5p-4 }, /* i=66 121.8 */
  { 0x1.7a11f4c07f71ap+0, 0x4.458d4425cf0f9a136p-60, -0x2.64ff39bfec2a6p-4, 0xd.12369a22b145f51a2p-60, 0xa.f10e7a90d50bp-8, 0x3.97e99ff224eff115c4p-60, -0x4.314557bc7e488p-8, -0x1.f6f30f9b7a3e395841p-60, 0x1.d92dca2b755dbp-8, 0x1.ea6af080340e8fb742p-64, -0xe.27aeb975dd95p-12, -0x2.f59f06dcd868fdd781p-64, 0x7.269315c81e518p-12, 0x1.fa9b6b1ac0389908568p-64, -0x3.c0b10947ba90ap-12, 0xc.5060d669958b980528p-68, 0x2.06ae97a5fceaep-12, -0x1.1e3f5faa715b2p-12, 0xa.0c1c848a2201p-16, -0x5.b92cd8af16958p-16, 0x8.7p-4 }, /* i=67 121.8 */
  { 0x1.79c5807c14aa5p+0, -0x5.d41ec3abaf6020fc0ap-56, -0x2.62461767553e4p-4, -0xc.b4beb0aa360d8f81f8p-60, 0xa.d812f0719a6cp-8, 0x2.0877f73e7cf09386d1p-60, -0x4.229f05a0aa27cp-8, -0x9.1f225d0fd29d4aebcp-64, 0x1.d06f8c7fbf422p-8, 0x4.56701e09bcaa23f534p-64, -0xd.d3178b90929bp-12, -0x1.f0231a7598c838d5d4p-64, 0x6.f2e96c6b9e24p-12, 0x1.84205ebdd18672b42fp-64, -0x3.a0e4cfeb3781ep-12, -0x2.766c996af34a7268ap-68, 0x1.f2fd58e07f381p-12, -0x1.11fd949cf46adp-12, 0x9.9188ad87842b8p-16, -0x5.6c4a863dc25b4p-16, 0x8.9p-4 }, /* i=68 121.9 */
  { 0x1.797962f87b74bp+0, -0x5.f1b3abc49ba293cbe4p-56, -0x2.5f932905eabecp-4, 0xf.da270b58f6e6872ee4p-60, 0xa.bf6e7c0d30dp-8, 0x3.99eaf80561839acc2cp-60, -0x4.143dd3f8913ap-8, 0x1.5346632151165938018p-60, 0x1.c7e56dbb9cc71p-8, -0x2.cd53fa51407a678bd2p-64, -0xd.80e204af622e8p-12, 0x1.74639a0af3599edae7p-64, 0x6.c0f4706ff1b6p-12, -0x1.21953f343f8aa9cccdp-64, -0x3.824ce1d1667c8p-12, 0x9.304b5c9d2f066b019p-68, 0x1.e0236fd39f661p-12, -0x1.0650f4f8310ap-12, 0x9.1d5f757900b9p-16, -0x5.23ccb6122895p-16, 0x8.bp-4 }, /* i=69 122.0 */
  { 0x1.792d9b708f1c5p+0, 0x7.368c5a300dbb09ee88p-56, -0x2.5ce65909aecd8p-4, 0x7.85070e8eda187b16c4p-60, 0xa.a71f8377f8abp-8, 0x7.5cdc8f8b1302438fp-64, -0x4.062027b0336f4p-8, -0x1.1c55d229efa5e0e3da8p-60, 0x1.bf8df7275cb1cp-8, 0x6.d5a861bda4f6a026a4p-64, -0xd.30fa10662875p-12, 0x2.765a119c5cfe7a1064p-64, 0x6.90a3a26a95db8p-12, 0x1.19a6d296ebd71d7888p-64, -0x3.64dc1c0bc4c24p-12, -0x8.71d81e519084923484p-68, 0x1.ce16a63c62266p-12, -0xf.b31b6a3a631e8p-16, 0x8.af433cca830ep-16, -0x4.df6df81631d94p-16, 0x8.dp-4 }, /* i=70 122.1 */
  { 0x1.78e22921d6d1bp+0, 0x7.8e6707b366f34a06fap-56, -0x2.5a3f9245ed0ccp-4, 0x6.b5075ddc8fcc8995e4p-60, 0xa.8f2476460ca88p-8, 0x1.5ccf72638d360c734dp-60, -0x3.f8447137be66ap-8, 0x8.f95dc2f12b17ec09p-64, 0x1.b767be5b16aebp-8, 0x7.9bd44b2d9672d37cfp-64, -0xc.e34c5c03ebcb8p-12, 0xa.f64332967c51d7f96p-68, 0x6.61e736767bee4p-12, 0x4.9a6fd78a7ed4c5f96p-68, -0x3.4885fac3ea45p-12, -0xa.f7479cc774ddfd43a4p-68, 0x1.bccd4e065276p-12, -0xf.09880e722796p-16, 0x8.46dc292438568p-16, -0x4.9eed774574e6cp-16, 0x8.fp-4 }, /* i=71 122.2 */
  { 0x1.78970b4c79267p+0, 0x3.265a87154419ef8974p-56, -0x2.579ebff0e33acp-4, -0xf.1d26bc2099be836d2cp-60, 0xa.777bcd474c58p-8, -0x2.90c91c3b5066763c8p-64, -0x3.eaa92c22ec8bap-8, 0xe.57bf669839ea0046cp-68, 0x1.af7164c82cd01p-8, -0x3.341c734ee706b78d18p-64, -0xc.97c64e5be459p-12, -0x3.83ac59015abb25ae11p-64, 0x6.34b00bbb27bfp-12, 0x1.885a7789a8097565d5p-64, -0x3.2d3e90f713438p-12, 0xe.b0f7f9cfe24e50a23cp-68, 0x1.ac3e399227916p-12, -0xe.67e657330371p-16, 0x7.e3d7c267e9574p-16, -0x4.620ea5fd69d1p-16, 0x9.1p-4 }, /* i=72 122.3 */
  { 0x1.784c41332fd3cp+0, 0x7.5fc971d200f945be8ap-56, -0x2.5503cda17a846p-4, 0x8.da9d6e582337274758p-60, 0xa.60240a459fa8p-8, -0x1.ac4fe19f084caae082p-60, -0x3.dd4cdecc03f4ep-8, 0xe.1380019445efd3ca38p-64, 0x1.a7a99747d031p-8, 0x2.4772395aa1097fea6cp-64, -0xc.4e55fff18512p-12, 0x3.474b75352a1abd40bbp-64, 0x6.08efa4643ce6cp-12, 0x5.6c9455b86155b04c3p-68, -0x3.12fa80a96ac24p-12, 0x8.d6b52520ff788e6f64p-68, 0x1.9c60b475a23f5p-12, -0xd.cdcd9df3ad39p-16, 0x7.85e89a1754f74p-16, -0x4.2898f105399d4p-16, 0x9.3p-4 }, /* i=73 122.3 */
  { 0x1.7801ca1b3bc8bp+0, -0x4.38ddc40444b8f850eep-56, -0x2.526ea74d11006p-4, -0x7.9016643506f6a118bp-60, 0xa.491bb7c55dedp-8, -0xe.c6a534566df9cf124p-64, -0x3.d02e19fa4dd8ap-8, 0xa.7d133899199c30a59cp-64, 0x1.a00f0dae4c7a6p-8, 0x3.a6ac83b282ee4165ecp-64, -0xc.06ea3382676a8p-12, -0x2.29ded2716725d44382p-64, 0x5.de981e0322bf4p-12, 0x1.60e3dc65d1da98d7cp-64, -0x2.f9aef389c303cp-12, -0xe.fa31eb07c67e5a2ad4p-68, 0x1.8d2c7caca7a99p-12, -0xd.3adb19520591p-16, 0x7.2cc5f39341a9cp-16, -0x3.f25777adee846p-16, 0x9.5p-4 }, /* i=74 122.4 */
  { 0x1.77b7a54c597ddp+0, -0x4.37614ae12d5a8124f4p-56, -0x2.4fdf394552d0cp-4, 0x9.1b9e33d907eabef02p-60, 0xa.326168c7c40f8p-8, -0x1.e04b36084aea51bp-60, -0x3.c34b788be14aap-8, -0x2.e26c5b3ec96325e9ccp-64, 0x1.98a08a62e120ep-8, -0x4.b8b5332fc43c419ca8p-64, -0xb.c1724ee944d38p-12, -0x3.ee65c53a209a5f4fc3p-64, 0x5.b59c2a54c9744p-12, -0x1.15a783cbd91efd8f108p-64, -0x2.e15193fee44a2p-12, -0x5.d6f7abbbaba01948cp-68, 0x1.7e99bc342c7bep-12, -0xc.aeb1822814ba8p-16, 0x6.d82b767eac1ap-16, -0x3.bf18c99f41a16p-16, 0x9.7p-4 }, /* i=75 122.5 */
  { 0x1.776dd210b58d7p+0, -0x7.1bfc77b3b477370c14p-56, -0x2.4d55703622646p-4, 0x1.d9c8db80b76bd2740cp-60, 0xa.1bf3b88f6658p-8, -0x3.83d592fce66f4a78fep-60, -0x3.b6a39f229dae2p-8, -0xe.0ac94eee640f3b8d1cp-64, 0x1.915cd9fbf224dp-8, 0x5.df8c8875f3b780f148p-64, -0xb.7dde54556e568p-12, 0x3.84121393476fe3d862p-64, 0x5.8def0865f7358p-12, 0x3.19136064b09ad2587p-68, -0x2.c9d886980cb1cp-12, -0xf.615d6eaa8df51bb4e8p-68, 0x1.70a1030809ff6p-12, -0xc.28f8be8b3974p-16, 0x6.87d8e8fc790cp-16, -0x3.8eaeaa7c767eap-16, 0x9.9p-4 }, /* i=76 122.6 */
  { 0x1.77244fb4e18d5p+0, -0x7.1eb4ca51ef90e020bap-56, -0x2.4ad139238f60ap-4, -0x9.263bc285ff2c243a3cp-60, 0xa.05d14a668b23p-8, 0x1.f937ecb8cabbc73925p-60, -0x3.aa353bd432f06p-8, -0xb.db3b507becdfb405bcp-64, 0x1.8a42d2df5ce0fp-8, 0x4.4a0aa91a2255ceecaap-64, -0xb.3c1edbd273d9p-12, -0x2.e9cea79b820eba105p-64, 0x5.67847e10db328p-12, -0x1.e26cda7d342854a02ap-64, -0x2.b33a63daa67a4p-12, 0x2.b576ae53124f6acc04p-68, 0x1.633b417d535d3p-12, -0xb.a95d924c513b8p-16, 0x6.3b91e7985bd48p-16, -0x3.60edd84199d94p-16, 0x9.bp-4 }, /* i=77 122.7 */
  { 0x1.76db1d87c928ep+0, 0x2.cfe7d3eb3ad9b6749cp-56, -0x2.48528167dbb8p-4, 0x2.b41955bcc59288dc8p-60, 0x9.eff8c9675c8d8p-8, -0x2.69ea4b89092fb0c7bp-64, -0x3.9dff05dd173p-8, -0x8.68144c3a1c394677fp-64, 0x1.835154e6c011cp-8, 0x3.2815a9b9782f1892c8p-64, -0xa.fc250d1be9bdp-12, -0x3.9bbc400c092ab8c461p-64, 0x5.4250d1cef4f3cp-12, 0x1.08bf20dd1d279960fb8p-64, -0x2.9d6e32679d2b8p-12, 0xf.6de3f790e2357a6acp-72, 0x1.5661c2f32f3fep-12, -0xb.2f915497591a8p-16, 0x5.f31da7767c774p-16, -0x3.35add71cfda42p-16, 0x9.dp-4 }, /* i=78 122.7 */
  { 0x1.76923adaa77a1p+0, -0x3.ad86440ee90b38e5f8p-56, -0x2.45d936b18e7b4p-4, -0x3.7956442f87501a90e4p-60, 0x9.da68e845dfep-8, 0xc.d007fb4215a735bedp-64, -0x3.91ffbd564af2cp-8, -0xf.f71fd0eb07d0820ce4p-64, 0x1.7c87490788ba2p-8, -0x7.dcd73f5472d61189ep-64, -0xa.bde299b97507p-12, -0x1.c3008be90c8fcc5554p-64, 0x5.1e48c4daaa6fcp-12, 0x1.0bf482ceaef1d30a2ap-64, -0x2.886b6173164b8p-12, -0x1.e9db741df7c66f6008p-68, 0x1.4a0e28d4a7008p-12, -0xa.bb49aa583f0a8p-16, 0x5.ae46bbb98f774p-16, -0x3.0cc8c0d6ccb0ap-16, 0x9.fp-4 }, /* i=79 122.8 */
  { 0x1.7649a700fc9e7p+0, 0x6.e76c8a8217af57eb7p-56, -0x2.4365470193f68p-4, 0x8.5531344ece0f8e7f4cp-60, 0x9.c520611ba244p-8, 0x1.ef17ad8cbe5d76541ap-60, -0x3.86362aedce744p-8, -0x4.2642efa23067b8fb44p-64, 0x1.75e3a0fe97dd4p-8, -0xb.ec152cd6b4a924eacp-68, -0xa.8149b75f7c648p-12, 0x3.52a2e5d806aefce17ep-64, 0x4.fb618d9c2e09p-12, -0xb.dc9b52f8ef1d1cc3d8p-68, -0x2.7429c389a83c2p-12, 0xf.382b60eb012d1ed564p-68, 0x1.3e3a65e63cd52p-12, -0xa.4c4045113368p-16, 0x5.6cdadffbf577p-16, -0x2.e61b17ccbe7fcp-16, 0xa.1p-4 }, /* i=80 122.9 */
  { 0x1.7601615083886p+0, -0x5.76787f5f94abd2f7e4p-56, -0x2.40f6a0a96aadp-4, -0xa.c41d85522e445be758p-60, 0x9.b01df5350ae98p-8, 0x3.7dd9f401af7d3bf51cp-60, -0x3.7aa11fa1abef2p-8, -0x4.36d4e13bf40ec85e6cp-64, 0x1.6f6556ff573d2p-8, 0x3.a09c6a9eb3e1411bp-64, -0xa.464d1a910f368p-12, 0x1.7c9d0b1fb1bf219a6dp-64, 0x4.d990d25d9804p-12, 0xb.9486d3899eec790c3p-68, -0x2.60a1899e79dfp-12, 0x6.959febd39aed34e298p-68, 0x1.32e0b9da767e4p-12, -0x9.e232a5d4e571p-16, 0x5.2eaac1cc10f64p-16, -0x2.c1839cb42cebp-16, 0xa.3p-4 }, /* i=81 123.0 */
  { 0x1.75b969212808cp+0, -0x7.d077f3826218b69862p-56, -0x2.3e8d32495cceap-4, -0x2.9f4c46011788fc8124p-60, 0x9.9b606ce04375p-8, -0x2.652b4d21d7c895c2f2p-60, -0x3.6f3f747d7c0a2p-8, 0xe.7690591e67c0ef4554p-64, 0x1.690b6d661560cp-8, 0x7.b441a4b288ca802af4p-64, -0xa.0cdff17fc3aep-12, 0x8.ed2421c1e8f52b4f4p-68, 0x4.b8cca4425518cp-12, 0x1.a1e81167a6aac95e368p-64, -0x2.4dcb3e5e002cap-12, 0x2.383c994bcd8f422824p-68, 0x1.27fbad28db261p-12, -0x9.7ce1e41cb211p-16, 0x4.f389d3e38729cp-16, -0x2.9ee3281e83faep-16, 0xa.5p-4 }, /* i=82 123.1 */
  { 0x1.7571bdccfd11bp+0, -0x4.9b83cc71f549f0b6d8p-56, -0x2.3c28eacec5c08p-4, 0x4.348988e2f26540b5ep-60, 0x9.86e6973da822p-8, 0x3.b5d85adbeb9c0a1f72p-60, -0x3.64100a5a4ac56p-8, 0x9.9d48271e3a2b3cb8p-68, 0x1.62d4ee6d83208p-8, 0x3.bdc61dc186e314e36p-64, -0x9.d4f5df267a088p-12, 0x8.f233cbbc529be0463p-68, 0x4.990b7a7e45eep-12, 0x1.d378189d119b1dcc958p-64, -0x2.3b9fc1c154682p-12, -0xe.622373479646382144p-68, 0x1.1d860d2334bf4p-12, -0x9.1c12783804868p-16, 0x4.bb4e20c6386dp-16, -0x2.7e1c866cc141p-16, 0xa.7p-4 }, /* i=83 123.1 */
  { 0x1.752a5eb0332f5p+0, -0x1.8f4de18b34f96a732ep-56, -0x2.39c9b972635b4p-4, -0x1.5121cde8f98c73d0b4p-60, 0x9.72af4a11b1aa8p-8, -0x1.b18007eeaf3e7c3ccap-60, -0x3.5911c9a0c46eap-8, -0x8.889cc9c65fe6db9828p-64, 0x1.5cc0ebe72ed35p-8, 0x2.344ec92b00756b6bbcp-64, -0x9.9e82f69c1d6ep-12, 0x3.eb23decc10506761a6p-64, 0x4.7a442dc921d1p-12, 0xf.0707c32d5956556d8p-72, -0x2.2a1844de5e696p-12, -0xd.e34774df2482cf10cp-68, 0x1.137ae8452b69ap-12, -0x8.bf8c0916f27cp-16, 0x4.85d01efdff1acp-16, -0x2.5f1455d80c258p-16, 0xa.9p-4 }, /* i=84 123.2 */
  { 0x1.74e34b290f35cp+0, 0x5.964e8d013a9900312p-60, -0x2.376f8db6b280ap-4, 0x5.02db14e778b012d1f8p-60, 0x9.5eb961984b7ap-8, -0x3.f39dd31bf0cb83fe9ep-60, -0x3.4e43a20d9336cp-8, -0x9.a5b6f020460fce318cp-64, 0x1.56ce7ef6dafc9p-8, -0x6.0fbb3527c41ab3feb4p-64, -0x9.697bb69ba2d7p-12, 0x8.ea89f5ad709927dbdp-68, 0x4.5c6df40ae20b8p-12, 0xe.c6a4611206310d2568p-68, -0x2.192e45f14ae4ap-12, 0x7.6fdbc8861dd1a87074p-68, 0x1.09d58aba93558p-12, -0x8.67193d3685d1p-16, 0x4.52ea906ca862cp-16, -0x2.41b0e8666c57ap-16, 0xa.bp-4 }, /* i=85 123.3 */
  { 0x1.749c8297e122ep+0, 0x7.6051f418afdb27ce28p-56, -0x2.351a576656afp-4, 0x4.37fae376c3762123dp-60, 0x9.4b03c059895d8p-8, -0x1.e5cd45d51b94de18f4p-60, -0x3.43a48a77d7112p-8, -0x5.fcdb542be076989218p-64, 0x1.50fcc7d090054p-8, 0x4.874da5515c947c3252p-64, -0x9.35d5053eac068p-12, 0x3.7a81b4dc38a0a7083dp-64, 0x4.3f805c3e2b1f8p-12, 0x1.ae9c9d9647d067baf68p-64, -0x2.08db8c9bfdfdcp-12, 0xe.330d31b7775b705514p-68, 0x1.00917b190956p-12, -0x8.12878e794f1ap-16, 0x4.227a571f5a118p-16, -0x2.25da25f468b8ep-16, 0xa.dp-4 }, /* i=86 123.4 */
  { 0x1.7456045efb329p+0, 0x5.a85f01650a0c59a7cap-56, -0x2.32ca06928c3d6p-4, -0x1.7d4ad3799e31897e5cp-60, 0x9.378d4effb0478p-8, -0x2.c7a0e3668cbbb6127ep-60, -0x3.39338099a2806p-8, 0x1.e2f4f32d8601327e58p-64, 0x1.4b4aed794a21cp-8, -0x4.7b5b345f330d61c364p-64, -0x9.03842be858078p-12, 0x3.57ee03d038450a2df5p-64, 0x4.23734a85d612cp-12, 0x1.791b1236329b1ecf4a8p-64, -0x1.f91a26586716ep-12, 0x3.4e6fb89613f92ffaf2p-68, 0xf.7aa7749ac75c8p-16, -0x7.c1a720ba90fap-16, 0x3.f45e579a5ff4cp-16, -0x2.0b7971f710cfcp-16, 0xa.fp-4 }, /* i=87 123.4 */
  { 0x1.740fcfe2a923cp+0, 0x3.d8938b59feccc93fp-60, -0x2.307e8b91a4ebcp-4, 0x1.b76b28ca17bdd78ddp-60, 0x9.2454fc2e864a8p-8, 0x1.df741e981adecf6606p-60, -0x3.2eef88da67e5ap-8, 0xd.63f456c462bb9f7854p-64, 0x1.45b81d8a26027p-8, 0x1.fc6938e6773411ced4p-64, -0x8.d27ed35dec16p-12, -0x1.63a83da6bad31a7e77p-64, 0x4.083ef472e617cp-12, -0x1.755b5240041d9816dbp-64, -0x1.e9e4631acbedbp-12, -0x7.0e0bbda869c355f60ap-68, 0xe.f1c719fff9848p-16, -0x7.744a9aede6588p-16, 0x3.c87758aee23a2p-16, -0x1.f27992741d949p-16, 0xb.1p-4 }, /* i=88 123.5 */
  { 0x1.73c9e48927ad8p+0, 0x2.f862a7522d2c91ea66p-56, -0x2.2e37d6fd8e77ap-4, -0xd.aeac0f3f434764f224p-60, 0x9.1159bc5bde67p-8, -0x1.77330b5e73a8fbf971p-60, -0x3.24d7ae1b43ccap-8, 0xc.62c30cc4e02869c7d4p-64, 0x1.40438bf5f05c4p-8, 0x5.62798e96844114fdfp-64, -0x8.a2bb000b1f48p-12, -0x3.12d39fdb9dedbd298dp-64, 0x3.eddbdd785e646p-12, -0xd.ccc573cdf57ed0d7c4p-68, -0x1.db34d2214f011p-12, 0x2.6726941687b6a0a1acp-68, 0xe.6e38e1b2c97ep-16, -0x7.2a4702a07ac38p-16, 0x3.9ea7e636fa49ap-16, -0x1.dac698d2c915p-16, 0xb.3p-4 }, /* i=89 123.6 */
  { 0x1.738441ba9c22bp+0, 0x5.942d2e2f8ce10954c4p-56, -0x2.2bf5d9b262e7ep-4, -0x2.0f834eed89d14e33dcp-60, 0x8.fe9a89a9552e8p-8, -0x3.f13712c833b9a547cfp-60, -0x3.1aeb01851189ap-8, 0x7.d28c9f9316df050758p-64, 0x1.3aec72d0fdaabp-8, -0x6.e8a19ec9fe56c0e522p-64, -0x8.742f0e6ff1038p-12, -0x8.65c7adbac45edf91ap-68, 0x3.d442d38a9122ep-12, 0xd.09e3b2c96e65f8e9bp-68, -0x1.cd063eee20191p-12, 0x4.41afa65775782f5768p-68, 0xd.efc1fcf124c3p-16, -0x6.e37399b3ef76cp-16, 0x3.76d433db595cp-16, -0x1.c44dcc0d95615p-16, 0xb.5p-4 }, /* i=90 123.7 */
  { 0x1.733ee6e10c42ep+0, -0x2.a441245403ac9589cp-60, -0x2.29b884cd02436p-4, -0x7.6251a2331425418cd8p-60, 0x8.ec1663bf33bdp-8, 0x1.15ca1b4065582ae689p-60, -0x3.11289a58384cep-8, 0x5.d59f11763ee7f3aa9cp-64, 0x1.35b2121b30d57p-8, 0x2.1b280dc35515366362p-64, -0x8.46d1afb60e518p-12, -0x1.7793b2692f029d119ep-64, 0x3.bb6cebe7a339p-12, -0x6.9987efd3fc3d0d23ccp-68, -0x1.bf53ae69edde5p-12, 0x3.e3dfdaa639e860e9b6p-68, 0xd.762a672b0a06p-16, -0x6.9fa9be2dcb898p-16, 0x3.50e206193302p-16, -0x1.aefd94f8a1336p-16, 0xb.7p-4 }, /* i=91 123.7 */
  { 0x1.72f9d36856376p+0, 0x3.d41b9e85b2e8b21d5cp-56, -0x2.277fc9a9b56c6p-4, -0xc.a6864301a590ee6578p-60, 0x8.d9cc4fa86ee5p-8, 0x2.6734f2f22de5b1d94ep-60, -0x3.078f95be1f964p-8, -0x8.bb1761fab9f35455p-64, 0x1.3093af8c1898cp-8, -0x6.e82933430d6a6c3a16p-64, -0x8.1a99e66be4898p-12, -0x1.6a0085b64cc6eae686p-64, 0x3.a3538007215b4p-12, 0xd.65f7ae4f80d2f9e9b8p-68, -0x1.b2185c2c54e1fp-12, 0x8.6fd1260fb789ec57p-72, 0xd.013cc11b12ebp-16, -0x6.5ec4cbf8ef0b4p-16, 0x3.2cb8964d984d2p-16, -0x1.9ac56aef5f783p-16, 0xb.9p-4 }, /* i=92 123.8 */
  { 0x1.72b506be28bbcp+0, -0x1.aa226731ca9f314574p-56, -0x2.254b99e2d9d8cp-4, -0xa.c0413d33f90214adp-64, 0x8.c7bb57afb8d1p-8, -0x1.48d1cd2af52db7706ap-60, -0x2.fe1f169c3abp-8, 0xf.8a0fb836ca202b0c58p-64, 0x1.2b9096610cb27p-8, -0x2.1ab98df2bd4ce2513ep-64, -0x7.ef7f036398bccp-12, -0x6.7348980d74ca667b48p-68, 0x3.8bf02aaea19dep-12, -0xc.cd3f40d94a5283eb54p-68, -0x1.a54fb7e82b808p-12, 0x7.ca9b99d888f26c95dp-68, 0xc.90c62deebcafp-16, -0x6.20a20078d807cp-16, 0x3.0a40810a4c554p-16, -0x1.8795c34d7c071p-16, 0xb.bp-4 }, /* i=93 123.9 */
  { 0x1.72708051fb6fcp+0, 0x5.c9e2fc1548c73359dep-56, -0x2.231be74f95e4ap-4, 0xe.8a6d89342cbab44284p-60, 0x8.b5e28b3d9bc4p-8, -0x1.b9deb19773bd138c4p-64, -0x2.f4d645689b894p-8, 0x6.1694779a19a2ec374p-64, 0x1.26a8172d34ef2p-8, 0x7.545d1ceb6979ab9192p-64, -0x7.c578a2b432dacp-12, 0x1.5d4e04f92d5a29e6048p-64, 0x3.753cc529898eep-12, -0xd.eb29cc915cb9daae38p-68, -0x1.98f562f9aacd8p-12, -0x5.3dfdae7a9d1e7f891ap-68, 0xc.2496325eec4cp-16, -0x5.e5205fd0cd948p-16, 0x2.e963ae4c84f4p-16, -0x1.7560008c51bccp-16, 0xb.dp-4 }, /* i=94 123.9 */
  { 0x1.722c3f950752ep+0, -0x7.54529dfdb2775a694p-56, -0x2.20f0a402957e4p-4, -0x8.9f13b6d73d32501534p-60, 0x8.a440feb794fap-8, 0x1.a0bad18fb6cd12db59p-60, -0x2.ebb44ffffdf82p-8, -0xb.17f5d03dbf2469fb8p-64, 0x1.21d987ab65376p-8, -0x2.e43e3722384971a714p-64, -0x7.9c7ea8db61bdcp-12, 0x3.ccbb07051eb29fe9ep-68, 0x3.5f3364a22c12p-12, -0x1.766571489c6087518p-68, -0x1.8d052e14935f8p-12, 0x6.703d498e6b1ac8c0dap-68, 0xb.bc7e959af0cap-16, -0x5.ac209bc3035acp-16, 0x2.ca0d3b52237a4p-16, -0x1.641662f3828fp-16, 0xb.fp-4 }, /* i=95 124.0 */
  { 0x1.71e843fa3f66fp+0, 0x4.2467dd5d3ac3a32262p-56, -0x2.1ec9c248ceecap-4, 0xc.218a5a3d779ac10aap-60, 0x8.92d5cb6026fd8p-8, -0x3.4e4efd83f413aab30ep-60, -0x2.e2b8697d3d08ep-8, 0x2.3ff213250d3f1066ap-64, 0x1.1d244291bab86p-8, -0x4.2ee7f107e142883cap-64, -0x7.74893ffe525p-12, -0x5.87c7923b22913984p-68, 0x3.49ce579a8bb4cp-12, -0xb.d777eec7b02e8ec784p-68, -0x1.817b171089849p-12, 0x7.de70287b9c4cdd8c8ap-68, 0xb.585343ea2aeap-16, -0x5.7584fc0db28a4p-16, 0x2.ac296df9957d2p-16, -0x1.53abfb767ad52p-16, 0xc.1p-4 }, /* i=96 124.1 */
  { 0x1.71a48cf6497abp+0, 0x3.75cd762c6a2e5813c4p-56, -0x2.1ca734a84f71cp-4, -0x7.df0c9d8f40dfdfaeb4p-60, 0x8.81a00f37db208p-8, 0x1.85647c0834444b29e1p-60, -0x2.d9e1ca12249ecp-8, -0x4.0513a3f978db2e8c34p-64, 0x1.1887a766e731ap-8, 0x3.65df6f536100ae6c22p-64, -0x7.4d90d54827cccp-12, -0x1.585650e0496f35e5ef8p-64, 0x3.3508238326352p-12, -0xf.2dd0882cd123cb9088p-68, -0x1.765346d1ff48cp-12, 0x5.7b918ccb120c48eb9ap-68, 0xa.f7ea32ea8e81p-16, -0x5.4131482df315p-16, 0x2.8fa59f03c97cep-16, -0x1.44149d884721fp-16, 0xc.3p-4 }, /* i=97 124.1 */
  { 0x1.716119ff77194p+0, 0x3.56bc3e4590c80c446p-60, -0x2.1a88eddf0f91ep-4, 0x9.bbbc459bf5f06ed5b4p-60, 0x8.709eecdf2a1dp-8, -0x8.0c1987fe03753924ap-64, -0x2.d12faee1923f8p-8, -0x2.6daea2545967dc452p-64, 0x1.14031a591850bp-8, 0x3.b8af5210674150e914p-64, -0x7.278e1664b5e34p-12, -0x9.2f15b47761d6578f5p-68, 0x3.20db826e448ap-12, -0xa.a1b4781d4d571b4ebcp-68, -0x1.6b8a0f4e1bde2p-12, 0x1.2ee1b7a1858a47e8d4p-68, 0xa.9b1b47537c7p-16, -0x5.0f0ab271b0b1cp-16, 0x2.74702998edd8cp-16, -0x1.3544d33dde181p-16, 0xc.5p-4 }, /* i=98 124.3 */
  { 0x1.711dea8dbe9f4p+0, 0x2.9806250e4940f52f5ep-56, -0x2.186ee0e1cec2p-4, -0x6.7096addbcb93200728p-60, 0x8.5fd18b79441b8p-8, 0x2.bfa16b048be72b38fep-60, -0x2.c8a159dad86c8p-8, 0x9.0db6c2f864340f351p-64, 0x1.0f96041669db4p-8, -0x6.26e4e992308504b6bap-64, -0x7.0279ef162f35p-12, 0x1.cbbd207a83f9062b338p-64, 0x3.0d4360de5ff94p-12, -0xe.fee3d0bdd7eb9ccd34p-68, -0x1.611be9a8278adp-12, 0x2.add039dd1839656ef4p-68, 0xa.41c03c26243b8p-16, -0x4.def7c443aa674p-16, 0x2.5a78603dda726p-16, -0x1.2731d22287146p-16, 0xc.7p-4 }, /* i=99 124.3 */
  { 0x1.70dafe1ab473fp+0, -0x7.f78554098a0424efp-56, -0x2.165900daf64b8p-4, -0x4.94a6376e3a8b674834p-60, 0x8.4f37168fb0c3p-8, 0x2.c8177e5e3311d2faa1p-60, -0x2.c03611965876ep-8, 0xe.63210cfe8aecf5e7f8p-64, 0x1.0b3fd1a6d232cp-8, -0x4.fa039ca8b53d71862p-64, -0x6.de4d86e48b6bcp-12, 0x1.c7a363d1950a3b0612p-64, 0x2.fa3adbae40c8ep-12, -0x4.dd1db766732851955p-68, -0x1.5705746719575p-12, 0x7.010c05377d2e4dbbp-68, 0x9.ebb48b363018p-16, -0x4.b0e04b9dd167p-16, 0x2.41ae79eb948fep-16, -0x1.19d1701d5d446p-16, 0xc.9p-4 }, /* i=100 124.4 */
  { 0x1.7098542184652p+0, -0x2.fa326dd1132bfe77a2p-56, -0x2.14474129832bap-4, -0xc.5225bad10fc42e494p-60, 0x8.3ecebdf6c0228p-8, 0x3.9a236b62a384e6fcf7p-60, -0x2.b7ed2133473acp-8, 0x3.5ba9261d3d47a2e78p-68, 0x1.06fff447797ddp-8, 0x6.3de479fcc9b9f00412p-64, -0x6.bb023ee577fa4p-12, -0x1.394f308b79378dcf25p-64, 0x2.e7bd3e118ded4p-12, -0xd.a1335e1cd9eef1d998p-68, -0x1.4d4371c1f86c5p-12, 0x7.0f9ca7458c10eec676p-68, 0x9.98d556f682598p-16, -0x4.84ad498fa7034p-16, 0x2.2a0387cef1338p-16, -0x1.0d1a19ac3a2a7p-16, 0xc.bp-4 }, /* i=101 124.5 */
  { 0x1.7055ec1eeb262p+0, -0x2.34dd6883b16beecbp-56, -0x2.1239955ff6becp-4, 0xb.57aca64f6a5459bf98p-60, 0x8.2e97b5b2c59c8p-8, 0x3.0d74a421f4e265f812p-60, -0x2.afc5d83691a8ep-8, -0x6.ae31fae8712c51b2b8p-64, 0x1.02d5e1476c7a9p-8, -0x4.bfd7c08a4e512d39p-68, -0x6.9891af9bb595cp-12, 0x5.55477ec7b08c1fa8e8p-68, 0x2.d5c5ffac96cfep-12, 0x5.c760cad9232f4d2a58p-68, -0x1.43d2c601d678ep-12, 0x6.b6e988ccd170c9c65cp-68, 0x9.490155775b9a8p-16, -0x4.5a48e1c7770c8p-16, 0x2.136966ceecee8p-16, -0x1.0102c88bd54fap-16, 0xc.dp-4 }, /* i=102 124.5 */
  { 0x1.7013c5912fef8p+0, 0x7.948d424df5a55f3704p-56, -0x2.102ff1434e032p-4, 0xd.5c38e69e07d72fc2p-60, 0x8.1e9135de1631p-8, -0xb.f183eaf8f7790f08dp-64, -0x2.a7bf8a6ad67a8p-8, -0xf.3fc27e800f20bdafd4p-64, 0xf.ec111e59cab58p-12, -0x5.55e60914a1851de66p-68, -0x6.76f5a6ecd2a2cp-12, -0x1.5bec5e4a195838c819p-64, 0x2.c450c2c23102p-12, 0x4.15d6e5cde25727dbap-68, -0x1.3ab075f83aaefp-12, 0x5.0f169e32bfee65ff54p-68, 0x8.fc18bc8430c28p-16, -0x4.319e4b0e6b3d4p-16, 0x1.fdd2b52c1089cp-16, -0xf.582fb08f7f1p-20, 0xc.fp-4 }, /* i=103 124.6 */
  { 0x1.6fd1dff81e402p+0, -0x5.9134ceae9d9839d7aap-56, -0x2.0e2a48c9ff458p-4, -0x8.196a1aa4d72712b648p-60, 0x8.0eba7a8fc3c7p-8, 0x3.43567fe2d44d815373p-60, -0x2.9fd98fc16adccp-8, 0x3.eba80699219863d164p-64, 0xf.ac1033020362p-12, -0x2.95e1e5d3c7e12937c6p-64, -0x6.5628262c4041p-12, -0x3.80575e19b6f0f65d98p-68, 0x2.b359527693672p-12, 0xd.c28b1bd03b6a061eecp-68, -0x1.31d9a588e5f94p-12, -0x1.aba66f9083de43aa7p-68, 0x8.b1fd2ed0b6388p-16, -0x4.0a99c0a882d74p-16, 0x1.e932c63b954d2p-16, -0xe.a92abcfe31d58p-20, 0xd.1p-4 }, /* i=104 124.7 */
  { 0x1.6f903ad4ffbd6p+0, 0x1.6822c081fed6975d4cp-56, -0x2.0c28901afe076p-4, -0x4.341eb435abad97c2f4p-60, 0x7.ff12c3c30f55p-8, -0x1.b145651c160e3fdc37p-60, -0x2.9813443460492p-8, -0xc.04def8aee284c50ef4p-64, 0xf.6d535e4a46a38p-12, -0x2.a237b81c8ce3625945p-64, -0x6.3623603acc4a8p-12, -0x1.8803944c4563a0ebbp-68, 0x2.a2dba126273e6p-12, -0x4.6696161d09de36cb98p-68, -0x1.294b9645e9153p-12, -0x1.1a28cec99e389ad0b4p-68, 0x8.6a91aa258bc18p-16, -0x3.e528748a73f36p-16, 0x1.d57d9922e564ap-16, -0xe.02a4a5bfb8548p-20, 0xd.3p-4 }, /* i=105 124.8 */
  { 0x1.6f4ed5aa96338p+0, -0x3.505c1d46315ef6dfb4p-56, -0x2.0a2abb8cc4f2p-4, 0x2.91bb355644d0d4bd0cp-60, 0x7.ef99553f8bf44p-8, -0x1.f5f1d08149dc76c2e58p-60, -0x2.906c07a9822a4p-8, -0xa.4bbfbad856d06c346cp-64, 0xf.2fd2e52067798p-12, -0x1.e4c0c8e3b16e60121dp-64, -0x6.16e1b7b996ac4p-12, 0x1.7bd7529dbb36717a0b8p-64, 0x2.92d3c6cf5568cp-12, 0x9.6e5b45ce5e015e409p-68, -0x1.2103a61d14467p-12, 0x6.2e7b309fa8e53de0b2p-68, 0x8.25ba767dfc558p-16, -0x3.c13882586581ap-16, 0x1.c2a7d013f3dd2p-16, -0xd.642b418d4c2d8p-20, 0xd.5p-4 }, /* i=106 124.7 */
  { 0x1.6f0daffd15b35p+0, -0x1.484aa516d8e6e46c6p-56, -0x2.0830bfa465abep-4, -0x4.cd7484a9799c1b89cp-64, 0x7.e04d7681ed38p-8, -0x1.020d1fe76efafc9fea8p-60, -0x2.88e33dd642488p-8, 0x5.c375728bb0f0eb1a7p-64, 0xe.f38743b06706p-12, -0x3.d775cc3d5ea66a0b88p-64, -0x5.f85dbd4fa3d9p-12, 0x1.a1399633b9817b7c55p-64, 0x2.833dff8e637fap-12, 0xf.45db62273952e054f4p-68, -0x1.18ff4e15d67bap-12, -0x4.8cf8408d638ede608p-68, 0x7.e35d1609319fp-16, -0x3.9eb8e3214dbdcp-16, 0x1.b0a6a23a0317ap-16, -0xc.cd52d16de84a8p-20, 0xd.7p-4 }, /* i=107 124.9 */
  { 0x1.6eccc9521ece5p+0, 0x2.cae6e35841f902e62cp-56, -0x2.063a91149e65ap-4, 0xb.947496e03731759f6p-60, 0x7.d12e72a57b3f4p-8, 0x1.68db92522924907712p-60, -0x2.81784e248b5c4p-8, -0x2.bfc8085413e18daf9cp-64, 0xe.b8692ba0935p-12, 0x2.aa303ebd116cf90cbp-64, -0x5.da922e0138e74p-12, -0x1.4bdcd9b43e6d3fe089p-64, 0x2.7416aa2a80b9ep-12, -0x3.5f655e1bca9dbdcp-80, -0x1.113c211eaeff9p-12, 0x2.38f84e1a9c68c6b66ep-68, 0x7.a3603601cd21p-16, -0x3.7d9961cb81aa2p-16, 0x1.9f6fdab62bceep-16, -0xc.3db5aafcad9bp-20, 0xd.9p-4 }, /* i=108 124.9 */
  { 0x1.6e8c2130b8ef7p+0, -0x7.4fe72109176175708ap-56, -0x2.044824bcf4f54p-4, 0xd.543789a246f9e71428p-60, 0x7.c23b984e273b4p-8, 0x1.7fcdbd44ef2b12046ap-60, -0x2.7a2aa3987179ep-8, 0xd.adf3d55c30f80be908p-64, 0xe.7e71825de6cbp-12, -0x2.f0e9902e713c13995cp-64, -0x5.bd79f1983820cp-12, 0x1.d188beec6253a67226p-64, 0x2.655a46b32e5ecp-12, -0xa.92c17e4af8b5b3202p-68, -0x1.09b7cae961483p-12, 0x6.079dc26706af93eed6p-68, 0x7.65aba0501d20cp-16, -0x3.5dca9027cd994p-16, 0x1.8ef9c983ad9e9p-16, -0xb.b4f3cf5dee878p-20, 0xd.bp-4 }, /* i=109 125.0 */
  { 0x1.6e4bb7214ccfep+0, 0x5.fec2690f3b325f8b3p-56, -0x2.02596fa8d7596p-4, 0xb.d70c342f67068897dp-60, 0x7.b37439933b4dp-8, -0x1.94516f3b25040baf268p-60, -0x2.72f9acb6b85f6p-8, 0x5.7acbe6a9f32d81be9cp-64, 0xe.45995f780259p-12, 0x3.5f6f01e7d67890fb0ep-64, -0x5.a110191cbe78p-12, -0xa.f014b293e88d2f264p-68, 0x2.5705752d3a96ep-12, 0x6.1432b6ca7a400c7d14p-68, -0x1.02700ed5262d5p-12, 0x7.a8bee6652cd332e128p-68, 0x7.2a282deb4b25cp-16, -0x3.3f3dbca10e6aep-16, 0x1.7f3b405c94ef3p-16, -0xb.32b29bc541fp-20, 0xd.dp-4 }, /* i=110 125.0 */
  { 0x1.6e0b8aad9f084p+0, 0x2.9efd2b776c82e73738p-56, -0x2.006e670ec17d2p-4, -0x2.72f9076bb7e832c4a4p-60, 0x7.a4d7abeaa0ce8p-8, 0x6.8ab3eebe66dddeeb6p-64, -0x2.6be4db6c2801p-8, -0xa.60ebfe3c4e1e471a48p-64, 0xe.0dda0b0c22e7p-12, 0x3.1c9c98a503b362f03ap-64, -0x5.854fdd5d4b91p-12, -0xb.b7c82dc145fbc1ddp-72, 0x2.4914f44e8eep-12, 0x6.2163f31748597249d4p-68, -0xf.b62c6e620d0d8p-16, 0x1.a7352f03725a5b73e3p-68, 0x6.f0bfb9ef0aba4p-16, -0x3.21e4e87ef3bfp-16, 0x1.702b892264734p-16, -0xa.b69c73549a9d8p-20, 0xd.fp-4 }, /* i=111 125.2 */
  { 0x1.6dcb9b60cabbcp+0, -0x5.5b8613fcf2e4159972p-56, -0x1.fe87004f68188p-4, -0x6.f61b59b07aa9a3428ap-60, 0x7.96654814b84a4p-8, 0x1.b426d2f674cf7cd868p-60, -0x2.64eba4f5a7fb2p-8, -0x7.69bc9f0a01d21455a4p-64, 0xd.d72cfc3e7a92p-12, 0xe.e4d8b4cc600dd92dp-72, -0x5.6a349d85c6af4p-12, 0xa.e147c52a1efb7ff9f8p-68, 0x2.3b85a0481d464p-12, -0xe.d8cfb7be78188085ccp-68, -0xf.48de2c967dd88p-16, -0x1.b4551fb89ddc11c3bp-72, 0x6.b95d155b9cep-16, -0x3.05b2beb322b84p-16, 0x1.61c26257ffebfp-16, -0xa.4060794bdecbp-20, 0xe.1p-4 }, /* i=112 125.2 */
  { 0x1.6d8be8c73c5e4p+0, 0x4.4100f0cc3c0cc83f5p-56, -0x1.fca330f4e888ap-4, -0x1.d995259992e599f1e4p-60, 0x7.881c6a08be888p-8, 0xb.b72128aac182d5577p-64, -0x2.5e0d81c91add2p-8, 0x6.e848ab17dead60313p-64, 0xd.a18bd7c15c04p-12, -0x3.00ee4eb32dd083d668p-64, -0x5.4fb9ddc4bb2d8p-12, 0xd.35c55cd74533fb703p-68, 0x2.2e54719d40d98p-12, 0x2.ef7b95368f5f47f068p-68, -0xe.def66e4fd7c08p-16, -0x2.c2a86d1eb75857833fp-68, 0x6.83ebfb76ddb7cp-16, -0x2.ea9a8b2887feep-16, 0x1.53f7f3ee5cabp-16, -0x9.cfb23ea1f3a2p-20, 0xe.3p-4 }, /* i=113 125.2 */
  { 0x1.6d4c726eac944p+0, -0x3.0c005fdd4e5d5ed16p-56, -0x1.fac2eeb1fd7ddp-4, -0x6.cd9a4308f7c0c1df2p-64, 0x7.79fc70e1ba4b4p-8, 0x1.8042474f55f24cf0918p-60, -0x2.5749ed7ef2876p-8, 0xc.948eff7789900b7a7p-64, 0xd.6cf06e69acf08p-12, 0x2.1b3cc4976780e768fep-64, -0x5.35db45fe2eed8p-12, 0xb.a96403914814a6d4p-72, 0x2.217e7c07eea3p-12, 0x9.d534dae411c7b5064cp-68, -0xe.7856b731986d8p-16, -0x2.9211255b30e7daa5cdp-68, 0x6.505906c54cd6cp-16, -0x2.d090328d213bap-16, 0x1.46c4cd3dd31e2p-16, -0x9.644987d4c7fc8p-20, 0xe.5p-4 }, /* i=114 125.2 */
  { 0x1.6d0d37e61b2b3p+0, 0x5.622cea774714dc2984p-56, -0x1.f8e62f61385dap-4, -0xe.01b50ad49cb3d7ed4p-64, 0x7.6c04becbee76p-8, 0x1.96244448ff9deb8a418p-60, -0x2.50a066bc771a2p-8, 0xd.ba5dd91a799e4a4e9cp-64, 0xd.3954bbd01aa58p-12, -0x4.b2f185a15955c1c32p-68, -0x5.1c94a08b7ca18p-12, -0x1.c674faa437d5fe44f2p-64, 0x2.1500ed691cae8p-12, -0x3.9ba2314dd9fd708c5cp-68, -0xe.14e1ba832223p-16, -0x2.d86d7fef87a3d00db5p-68, 0x6.1e91a692dcc6p-16, -0x2.b7882a8f0df38p-16, 0x1.3a21ddc3c840cp-16, -0x8.fde20aad9148p-20, 0xe.7p-4 }, /* i=115 125.4 */
  { 0x1.6cce38bdca2b4p+0, 0x7.63bedd9ef8f1de529ep-56, -0x1.f70ce9043f382p-4, 0x3.19ea6c6a61389f622cp-60, 0x7.5e34b8f2cc82cp-8, -0x1.69841536fde98309db8p-60, -0x2.4a106f1eba33ep-8, -0xd.56071d5b8c95cc715p-64, 0xd.06b2e4fe91548p-12, -0x1.74c58d693d9e767ap-64, -0x5.03e1d907a2ddcp-12, 0x1.61ace928d4f4cb6b66p-64, 0x2.08d90cc4cc9aap-12, 0xf.9590b435d3d8104p-80, -0xd.b47b4e336f918p-16, -0x2.af7ed1802c0d781739p-68, 0x5.ee841503954ep-16, -0x2.9f77728622e0ap-16, 0x1.2e0870fae4c64p-16, -0x8.9c3b34d5c9e88p-20, 0xe.9p-4 }, /* i=116 125.5 */
  { 0x1.6c8f748739006p+0, 0x5.da3781e570bdd788acp-56, -0x1.f53711c30f30dp-4, -0x7.5ca626a9c385c7046ap-60, 0x7.508bc76f634d8p-8, -0x7.e23d3cd6b51b481578p-64, -0x2.43998b26307dap-8, 0x8.d9d97300eca00a3fe8p-64, 0xc.d50537297c3dp-12, -0x2.123b2409ac957a72fap-64, -0x4.ebbefb277f8f8p-12, -0xd.bf7caa869eedc1b3ap-68, 0x1.fd0439492e90cp-12, 0x7.797cdcb60c711abffap-68, -0xd.57085e77b027p-16, 0x1.412de9c41dceb983bfp-68, 0x5.c01f4da4ad91cp-16, -0x2.88538c83aaf58p-16, 0x1.227228b21b496p-16, -0x8.3f17f3a7ba18p-20, 0xe.bp-4 }, /* i=117 125.5 */
  { 0x1.6c50ead51fbb2p+0, -0x6.3e63d5db2040214862p-56, -0x1.f3649feb433d2p-4, -0x2.94f737b2a19a36dd88p-60, 0x7.430955374664p-8, -0xb.c5cc2ea76057df299p-64, -0x2.3d3b4222ebbbcp-8, 0x5.912f204ec5dd36a8bp-64, 0xc.a44626745a038p-12, 0x3.6866c5703cafa77a1dp-64, -0x4.d428319d75d98p-12, 0x7.d496378f6016448a8p-68, 0x1.f17fe960588aap-12, -0x6.0f7c4384cc8a89cb08p-68, -0xc.fc6ee1fcfe09p-16, 0x1.b14dfe63db3f55fe5p-68, 0x5.9353047720428p-16, -0x2.721276c26497ep-16, 0x1.1758f787b308ep-16, -0x7.e63e7e57267bcp-20, 0xe.dp-4 }, /* i=118 125.5 */
  { 0x1.6c129b3b6a682p+0, 0x3.b51399eee26774b35ap-56, -0x1.f19589ef5f19cp-4, -0x5.a44b5db43db2ac66fcp-60, 0x7.35acd00bea23p-8, -0x1.967634028b474853b6p-60, -0x2.36f51e216fd9p-8, 0xa.beb543f532ebe1aap-68, 0xc.74704cc13584p-12, -0x3.de9c90904a2a411495p-64, -0x4.bd19c50802518p-12, -0x3.d10b2eeb3dafa2dedp-68, 0x1.e649a9cc14727p-12, 0x3.4d433fa89cdd54f5ecp-68, -0xc.a495cea5cb318p-16, 0x6.48942fb7c66adcc08p-72, 0x5.680f9d6d21a68p-16, -0x2.5caaa57126ee8p-16, 0x1.0cb71df9bc2afp-16, -0x7.91782847aa228p-20, 0xe.fp-4 }, /* i=119 125.7 */
  { 0x1.6bd4854f347e8p+0, -0x3.3a8337f26c9bfedf48p-56, -0x1.efc9c6661e5b7p-4, 0x5.e7294602d56a7a16ep-64, 0x7.2875a86a7113p-8, 0x1.0a5c5f3ffaaffe5286p-60, -0x2.30c6abd81d90ep-8, -0xc.8181e928677c6b398p-64, 0xc.457e688a981bp-12, -0x9.d4f076fb7ba433766p-68, -0x4.a6901aeac7658p-12, -0xf.50a7a2b6c0b2c8dbep-68, 0x1.db5f1ccb4df3dp-12, 0x2.346e4b0436d095ff52p-68, -0xc.4f650ecc0928p-16, -0x3.cd65404e295f92fae3p-68, 0x5.3e462454524a4p-16, -0x2.4812fcd2f58eap-16, 0x1.028725fa13a56p-16, -0x7.4091330bb3ca4p-20, 0xf.1p-4 }, /* i=120 125.8 */
  { 0x1.6b96a8a6c4633p+0, 0x5.0a970e7f2cd9ed9aaep-56, -0x1.ee014c09c77e1p-4, 0x7.d464d03b871440292cp-60, 0x7.1b63517be711cp-8, 0xb.07d499fba45384ca48p-64, -0x2.2aaf7a952d98cp-8, 0xc.e58d68a00a4968bc44p-64, 0xc.176b5bc78efbp-12, 0x3.29de1408c4a433783cp-64, -0x4.9087b4b1a1144p-12, 0x1.7682e4495ab33a6ff7p-64, 0x1.d0bdf948aec3cp-12, 0x1.38dcd162dd945a8008p-68, -0xb.fcc577016efa8p-16, 0x1.0eecb072f72eb815dp-68, 0x5.15e84526d5918p-16, -0x2.3442cbaf99448p-16, 0xf.8c3de5be7607p-20, -0x6.f358a21549f88p-20, 0xf.3p-4 }, /* i=121 125.8 */
  { 0x1.6b5904d98702bp+0, 0x1.c5c8c14396e8dd3a82p-56, -0x1.ec3c11b782d7dp-4, -0x5.57671e0774fad60926p-60, 0x7.0e754105e6f18p-8, 0x1.beb33416803cdb5a548p-60, -0x2.24af1c2d37552p-8, 0x4.bcc5f03c718719c83cp-64, 0xb.ea322ad961538p-12, 0x3.a1681b4c8092bc376ep-64, -0x4.7afd2ebd5458p-12, 0x1.09e5222f9f04fd7515p-64, 0x1.c6640a11fdc9fp-12, 0x4.37d9d39795ce05ccbp-68, -0xb.aca0bc479bae8p-16, 0x3.e21afb0a9b0ef2d2fap-68, 0x4.eee844bdc93b8p-16, -0x2.2131c61021a2ep-16, 0xe.f68575d780188p-20, -0x6.a9a011e3f68dp-20, 0xf.5p-4 }, /* i=122 125.7 */
  { 0x1.6b1b99800b7edp+0, 0x4.164447102409afcd2p-60, -0x1.ea7a0e6eb5564p-4, -0x7.ebf320a56064bfb1cp-64, 0x7.01aaef5ba75bcp-8, 0xf.ed6b1f4610752b63d8p-64, -0x2.1ec524ea3e5cep-8, 0xc.6f1ed7ccef62ba1fdp-64, 0xb.bdcdfb829944p-12, -0xb.00a642dc338f5c8abp-68, -0x4.65ed3f79738dcp-12, 0x1.5f84ad02da05e8a896p-64, 0x1.bc4f2d17dae72p-12, -0x5.a8c62428723945p-72, -0xb.5ee16aba28dbp-16, -0x2.1874803ca4519ceff7p-68, 0x4.c938f9dff1eacp-16, -0x2.0ed80042faeb4p-16, 0xe.66fdf2d360d18p-20, -0x6.633b93f3caef8p-20, 0xf.7p-4 }, /* i=123 125.9 */
  { 0x1.6ade6633fef17p+0, 0x5.d16674c9af0b94fd42p-56, -0x1.e8bb39505eeadp-4, 0x2.5f231a2b6f186bd012p-60, 0x6.f503d74f6bd8p-8, 0x1.9f5d6daa73b0630315p-60, -0x2.18f12b7b42348p-8, 0xf.7ae8c36cd99d65c018p-64, 0xb.923a13e705738p-12, 0x3.8bfca2b60d5b8068f2p-64, -0x4.5154b67b15cap-12, 0x1.83a44bb07f7e10930ep-64, 0x1.b27d52b576204p-12, 0x5.4bf61b3d7052a4fb14p-68, -0xb.1372dca5119e8p-16, 0x5.24547e2ffc5ab358p-72, 0x4.a4cdc6a7b610cp-16, -0x1.fd2dea23882d9p-16, 0xd.dd5ff8838a598p-20, -0x6.20018b8113e8cp-20, 0xf.9p-4 }, /* i=124 125.9 */
  { 0x1.6aa16a9028432p+0, -0x7.8414a0f94ed0af18fp-56, -0x1.e6ff899e7c8fap-4, 0x3.9e11e6af39b1b69c4cp-60, 0x6.e87f762456fa4p-8, -0x1.6d3e0c710666c64dc28p-60, -0x2.1332c8e44bdd8p-8, -0x9.c3bf9e9af144375d7cp-64, 0xb.6771d9945ceb8p-12, -0x2.cc5bc54a4c4735ede4p-64, -0x4.3d307ba80384cp-12, -0x3.d6fca78bdba7f15bcp-68, 0x1.a8ec7cffe5bd7p-12, 0x5.e000fcde0c04a6d8ap-72, -0xa.ca4132021e1e8p-16, -0x3.ee15f9c5ad67cb6b49p-68, 0x4.819a923bd50b8p-16, -0x1.ec2c4aa17005dp-16, 0xd.59678aeefe5bp-20, -0x5.dfca887f5932p-20, 0xf.bp-4 }, /* i=125 126.1 */
  { 0x1.6a64a6306414dp+0, -0x1.2728c849ad57b51148p-56, -0x1.e546f6bb6dcb8p-4, -0x6.78bcc9b64aceab1d7cp-60, 0x6.dc1d4b809acfp-8, 0xb.e063d6aca19e3be5f8p-64, -0x2.0d89986ef4fcap-8, -0xf.20e58a3606b2c2e1a4p-64, 0xb.3d70d09332588p-12, -0x3.61cb62faa1391338e1p-64, -0x4.297d8e65ff3c4p-12, -0xb.5cad6b584356f728dp-68, 0x1.9f9abf1cc3708p-12, 0x7.5b7c7241a81002a1e6p-68, -0xa.833948584b7cp-16, -0x2.6e6ff17b906b3b0fd1p-68, 0x4.5f93c2d66b9fcp-16, -0x1.dbcc3b84025fdp-16, 0xc.dad414bd0d0c8p-20, -0x5.a2712f93001cp-20, 0xf.dp-4 }, /* i=126 126.2 */
  { 0x1.6a2818c0bd494p+0, -0x4.56db480b4368ceeaf6p-56, -0x1.e39178965b6c5p-4, -0x4.20bb3f4eeda06f112ap-60, 0x6.cfdcdc6bf4918p-8, 0x1.854e6a9a87f2c589998p-60, -0x2.07f538fdead68p-8, 0xf.6663706c5be74c573p-64, 0xb.1432a4b775e7p-12, 0x3.37638caafde049b51p-64, -0x4.16390995676c4p-12, -0x1.951a08db4b9fb024f88p-64, 0x1.96863edd5ccb1p-12, 0x4.c97a59d57ef626d41ep-68, -0xa.3e48c3f324a28p-16, -0x2.f34dd5d39f24ddaa2ep-68, 0x4.3eae402c67548p-16, -0x1.cc072944ca63dp-16, 0xc.61682c3d319p-20, -0x5.67d221e8a35b4p-20, 0xf.effff8p-4 }, /* i=127 126.2 */
};

typedef union { double x; uint32_t i[2]; } union_t;

// begin_acospi
// exceptional cases for |x| < 0.5, should be sorted by increasing values
#define EXCEPTIONS 114
static const double exceptions[EXCEPTIONS][2] = {
    {-0x1.ac26cc49b0264p-3, 0x1.22530cb92c8f1p-1},
    {-0x1.f8374b05d52b2p-5, 0x1.0a09933a4589ap-1},
    {-0x1.7b9482807fe3cp-5, 0x1.078de0c9ad10cp-1},
    {-0x1.df5c9cfec76e3p-6, 0x1.04c4dbe40b566p-1},
    {-0x1.d8a6b018a513dp-7, 0x1.0259d1a82f985p-1},
    {-0x1.f0a606701a3e5p-9, 0x1.009e16a90fcbep-1},
    {-0x1.f7b97997721efp-10, 0x1.00502b98793edp-1},
    {-0x1.1c68c1ab0ea57p-10, 0x1.002d43e032ad6p-1},
    {-0x1.59ec5d1d2e86p-11, 0x1.001b87178f037p-1},
    {-0x1.2878e662f31bcp-12, 0x1.000bcbd77c281p-1},
    {-0x1.83d790f27f3e5p-14, 0x1.0003dba1b6271p-1},
    {-0x1.dd61dedc265f5p-15, 0x1.00025fd270a53p-1},
    {-0x1.29320aa46f06ep-15, 0x1.00017a66a9fb7p-1},
    {-0x1.99bc3ca2a27p-17, 0x1.0000826c3c68fp-1},
    {-0x1.8ec2e3d383d64p-19, 0x1.00001fbb7f8eap-1},
    {-0x1.7679a41f39c67p-19, 0x1.00001dccbecdcp-1},
    {-0x1.427b32daed80dp-20, 0x1.00000cd4c4593p-1},
    {-0x1.2a4da90d1ca6p-21, 0x1.000005ef3f1d5p-1},
    {-0x1.c0149992fc291p-22, 0x1.00000475070d1p-1},
    {-0x1.783beb73b1516p-23, 0x1.000001df093cep-1},
    {-0x1.5b9604ff5f529p-23, 0x1.000001ba8f6d1p-1},
    {-0x1.3a21ab6af61a5p-23, 0x1.0000018ff6f5p-1},
    {-0x1.71732d4eb3deap-25, 0x1.00000075997cbp-1},
    {-0x1.c68b45cd67cb4p-28, 0x1.0000001215f29p-1},
    {-0x1.dbbf2a5ec6647p-29, 0x1.0000000976f55p-1},
    {-0x1.d0e597b6ba3b6p-30, 0x1.000000049fd9ap-1},
    {-0x1.10c87354e03c8p-31, 0x1.000000015b517p-1},
    {-0x1.c277920f3fbadp-32, 0x1.000000011ec6bp-1},
    {-0x1.74c46decf3374p-32, 0x1.00000000ed4fap-1},
    {-0x1.65a1dd290660fp-36, 0x1.000000000e3adp-1},
    {-0x1.635e3d74befcap-38, 0x1.000000000388fp-1},
    {-0x1.635e3d74befcap-39, 0x1.0000000001c47p-1},
    {-0x1.6c6cbc45dc8dep-40, 0x1.0000000000e8p-1},
    {-0x1.5a4fbea3a16b6p-40, 0x1.0000000000dc7p-1},
    {-0x1.6c6cbc45dc8dep-41, 0x1.000000000074p-1},
    {-0x1.6c6cbc45dc8dep-42, 0x1.00000000003ap-1},
    {-0x1.6c6cbc45dc8dep-43, 0x1.00000000001dp-1},
    {-0x1.6c6cbc45dc8dep-44, 0x1.00000000000e8p-1},
    {-0x1.6c6cbc45dc8dep-45, 0x1.0000000000074p-1},
    {-0x1.6c6cbc45dc8dep-46, 0x1.000000000003ap-1},
    {-0x1.6c6cbc45dc8dep-47, 0x1.000000000001dp-1},
    {-0x1.921fb54442d18p-51, 0x1.0000000000002p-1},
    {-0x1.2d97c7f3321d2p-51, 0x1.0000000000001p-1},
    {-0x1.921fb54442d18p-52, 0x1.0000000000001p-1},
    {-0x1.921fb54442d18p-53, 0x1p-1},
    {0x0p+0, 0x1p-1},
    {0x1.921fb54442d19p-54, 0x1.fffffffffffffp-2},
    {0x1.3a28c59d5433bp-49, 0x1.ffffffffffff3p-2},
    {0x1.6c6cbc45dc8dep-49, 0x1.ffffffffffff1p-2},
    {0x1.3a28c59d5433bp-48, 0x1.fffffffffffe7p-2},
    {0x1.6c6cbc45dc8dep-48, 0x1.fffffffffffe3p-2},
    {0x1.6c6cbc45dc8dep-47, 0x1.fffffffffffc6p-2},
    {0x1.6c6cbc45dc8dep-46, 0x1.fffffffffff8cp-2},
    {0x1.21cfda23b228p-45, 0x1.fffffffffff47p-2},
    {0x1.6c6cbc45dc8dep-45, 0x1.fffffffffff18p-2},
    {0x1.6c6cbc45dc8dep-44, 0x1.ffffffffffe3p-2},
    {0x1.6c6cbc45dc8dep-43, 0x1.ffffffffffc6p-2},
    {0x1.b4e0b2cec917ep-43, 0x1.ffffffffffba7p-2},
    {0x1.6c6cbc45dc8dep-42, 0x1.ffffffffff8cp-2},
    {0x1.90a6b78a52d2ep-42, 0x1.ffffffffff807p-2},
    {0x1.6c6cbc45dc8dep-41, 0x1.ffffffffff18p-2},
    {0x1.6c6cbc45dc8dep-40, 0x1.fffffffffe3p-2},
    {0x1.67e57cdd4dc54p-39, 0x1.fffffffffc6b8p-2},
    {0x1.67e57cdd4dc54p-38, 0x1.fffffffff8d71p-2},
    {0x1.0dec1da5fa53fp-37, 0x1.fffffffff5429p-2},
    {0x1.93c05c9ed3cbcp-36, 0x1.ffffffffdfdedp-2},
    {0x1.39c6fd67805a7p-35, 0x1.ffffffffce0f9p-2},
    {0x1.a9adcc7f96cfp-35, 0x1.ffffffffbc405p-2},
    {0x1.39c6fd67805a7p-34, 0x1.ffffffff9c1f2p-2},
    {0x1.44bdb557e1dc1p-34, 0x1.ffffffff98a1bp-2},
    {0x1.bf9b3c6059d24p-34, 0x1.ffffffff7185cp-2},
    {0x1.ca91f450bb53ep-34, 0x1.ffffffff6e085p-2},
    {0x1.57bfc3f20ac4fp-32, 0x1.fffffffe4a533p-2},
    {0x1.74c46decf3374p-32, 0x1.fffffffe2560cp-2},
    {0x1.10c87354e03c8p-31, 0x1.fffffffd495d2p-2},
    {0x1.7ff87542efc93p-31, 0x1.fffffffc2e3a2p-2},
    {0x1.b09687e430521p-29, 0x1.ffffffeec9b1dp-2},
    {0x1.dbbf2a5ec6647p-29, 0x1.ffffffed12156p-2},
    {0x1.cc4070a3023a7p-28, 0x1.ffffffdb5fd43p-2},
    {0x1.8d69f36b1cce6p-27, 0x1.ffffffc0bfe6fp-2},
    {0x1.bdbf63f170fep-27, 0x1.ffffffb90e9bbp-2},
    {0x1.ac57910d7e867p-26, 0x1.ffffff77a78c9p-2},
    {0x1.71732d4eb3deap-25, 0x1.ffffff14cd06ap-2},
    {0x1.24e3786d7459ap-24, 0x1.fffffe8b152c5p-2},
    {0x1.94ee8417af86ep-23, 0x1.fffffbf8d9921p-2},
    {0x1.2a4da90d1ca6p-21, 0x1.fffff42181c56p-2},
    {0x1.769cbdfa9e1d8p-21, 0x1.fffff1183b0ffp-2},
    {0x1.427b32daed80dp-20, 0x1.ffffe656774dap-2},
    {0x1.73e9e2254c98fp-20, 0x1.ffffe267702edp-2},
    {0x1.8ec2e3d383d64p-19, 0x1.ffffc08900e2cp-2},
    {0x1.2226500110692p-18, 0x1.ffffa3a47a927p-2},
    {0x1.4c2c402a28e67p-18, 0x1.ffff96441ee89p-2},
    {0x1.59abdcf048a9bp-17, 0x1.ffff23f0548e1p-2},
    {0x1.5ef5605b74f7p-16, 0x1.fffe41255058bp-2},
    {0x1.969df9fda1b9p-16, 0x1.fffdfa476ec5cp-2},
    {0x1.29320aa46f06ep-15, 0x1.fffd0b32ac091p-2},
    {0x1.6638f8a2185cbp-15, 0x1.fffc6fcb3553bp-2},
    {0x1.dd61dedc265f5p-15, 0x1.fffb405b1eb59p-2},
    {0x1.41c4a3005149fp-14, 0x1.fff9993f9b0cfp-2},
    {0x1.17485bab7607bp-13, 0x1.fff4e34013127p-2},
    {0x1.a0ea629bbf42p-13, 0x1.ffef69553769ep-2},
    {0x1.2878e662f31bcp-12, 0x1.ffe8685107afep-2},
    {0x1.cf88e9268bdcfp-12, 0x1.ffdb1cf235d97p-2},
    {0x1.0d0c56d741955p-9, 0x1.ff54b7ee08148p-2},
    {0x1.f0a606701a3e5p-9, 0x1.fec3d2ade0684p-2},
    {0x1.730786d8eea35p-8, 0x1.fe27969253c05p-2},
    {0x1.bfa7fc1b1e67cp-8, 0x1.fdc605a201b7bp-2},
    {0x1.271adcddb2534p-7, 0x1.fd1082e4a2bcdp-2},
    {0x1.1f6bbb1fac712p-6, 0x1.fa481a48acacep-2},
    {0x1.20588554d6e24p-6, 0x1.fa43642302402p-2},
    {0x1.5fa453f6bb40ep-6, 0x1.f900f7075f3f5p-2},
    {0x1.69a68ac73eaa2p-4, 0x1.e32ee8bcc5f06p-2},
    {0x1.da83e335e379ap-3, 0x1.b3c936985f665p-2},
    {0x1.f608aa4e62781p-2, 0x1.58fc14707d797p-2},
  };
static const int8_t exceptions_rnd[EXCEPTIONS] = {
    1, /* -0x1.ac26cc49b0264p-3 */
    -1, /* -0x1.f8374b05d52b2p-5 */
    -1, /* -0x1.7b9482807fe3cp-5 */
    -1, /* -0x1.df5c9cfec76e3p-6 */
    -1, /* -0x1.d8a6b018a513dp-7 */
    -1, /* -0x1.f0a606701a3e5p-9 */
    1, /* -0x1.f7b97997721efp-10 */
    -1, /* -0x1.1c68c1ab0ea57p-10 */
    1, /* -0x1.59ec5d1d2e86p-11 */
    1, /* -0x1.2878e662f31bcp-12 */
    1, /* -0x1.83d790f27f3e5p-14 */
    1, /* -0x1.dd61dedc265f5p-15 */
    1, /* -0x1.29320aa46f06ep-15 */
    -1, /* -0x1.99bc3ca2a27p-17 */
    -1, /* -0x1.8ec2e3d383d64p-19 */
    -1, /* -0x1.7679a41f39c67p-19 */
    1, /* -0x1.427b32daed80dp-20 */
    1, /* -0x1.2a4da90d1ca6p-21 */
    -1, /* -0x1.c0149992fc291p-22 */
    -1, /* -0x1.783beb73b1516p-23 */
    -1, /* -0x1.5b9604ff5f529p-23 */
    -1, /* -0x1.3a21ab6af61a5p-23 */
    1, /* -0x1.71732d4eb3deap-25 */
    -1, /* -0x1.c68b45cd67cb4p-28 */
    1, /* -0x1.dbbf2a5ec6647p-29 */
    -1, /* -0x1.d0e597b6ba3b6p-30 */
    1, /* -0x1.10c87354e03c8p-31 */
    1, /* -0x1.c277920f3fbadp-32 */
    1, /* -0x1.74c46decf3374p-32 */
    -1, /* -0x1.65a1dd290660fp-36 */
    -1, /* -0x1.635e3d74befcap-38 */
    1, /* -0x1.635e3d74befcap-39 */
    1, /* -0x1.6c6cbc45dc8dep-40 */
    1, /* -0x1.5a4fbea3a16b6p-40 */
    1, /* -0x1.6c6cbc45dc8dep-41 */
    1, /* -0x1.6c6cbc45dc8dep-42 */
    1, /* -0x1.6c6cbc45dc8dep-43 */
    1, /* -0x1.6c6cbc45dc8dep-44 */
    1, /* -0x1.6c6cbc45dc8dep-45 */
    1, /* -0x1.6c6cbc45dc8dep-46 */
    1, /* -0x1.6c6cbc45dc8dep-47 */
    -1, /* -0x1.921fb54442d18p-51 */
    1, /* -0x1.2d97c7f3321d2p-51 */
    -1, /* -0x1.921fb54442d18p-52 */
    -1, /* -0x1.921fb54442d18p-53 */
    0, /* 0x0p+0 */
    1, /* 0x1.921fb54442d19p-54 */
    1, /* 0x1.3a28c59d5433bp-49 */
    1, /* 0x1.6c6cbc45dc8dep-49 */
    -1, /* 0x1.3a28c59d5433bp-48 */
    -1, /* 0x1.6c6cbc45dc8dep-48 */
    -1, /* 0x1.6c6cbc45dc8dep-47 */
    -1, /* 0x1.6c6cbc45dc8dep-46 */
    1, /* 0x1.21cfda23b228p-45 */
    -1, /* 0x1.6c6cbc45dc8dep-45 */
    -1, /* 0x1.6c6cbc45dc8dep-44 */
    -1, /* 0x1.6c6cbc45dc8dep-43 */
    1, /* 0x1.b4e0b2cec917ep-43 */
    -1, /* 0x1.6c6cbc45dc8dep-42 */
    1, /* 0x1.90a6b78a52d2ep-42 */
    -1, /* 0x1.6c6cbc45dc8dep-41 */
    -1, /* 0x1.6c6cbc45dc8dep-40 */
    1, /* 0x1.67e57cdd4dc54p-39 */
    -1, /* 0x1.67e57cdd4dc54p-38 */
    1, /* 0x1.0dec1da5fa53fp-37 */
    1, /* 0x1.93c05c9ed3cbcp-36 */
    -1, /* 0x1.39c6fd67805a7p-35 */
    -1, /* 0x1.a9adcc7f96cfp-35 */
    -1, /* 0x1.39c6fd67805a7p-34 */
    1, /* 0x1.44bdb557e1dc1p-34 */
    -1, /* 0x1.bf9b3c6059d24p-34 */
    1, /* 0x1.ca91f450bb53ep-34 */
    1, /* 0x1.57bfc3f20ac4fp-32 */
    -1, /* 0x1.74c46decf3374p-32 */
    -1, /* 0x1.10c87354e03c8p-31 */
    -1, /* 0x1.7ff87542efc93p-31 */
    1, /* 0x1.b09687e430521p-29 */
    -1, /* 0x1.dbbf2a5ec6647p-29 */
    1, /* 0x1.cc4070a3023a7p-28 */
    1, /* 0x1.8d69f36b1cce6p-27 */
    -1, /* 0x1.bdbf63f170fep-27 */
    -1, /* 0x1.ac57910d7e867p-26 */
    -1, /* 0x1.71732d4eb3deap-25 */
    1, /* 0x1.24e3786d7459ap-24 */
    1, /* 0x1.94ee8417af86ep-23 */
    -1, /* 0x1.2a4da90d1ca6p-21 */
    -1, /* 0x1.769cbdfa9e1d8p-21 */
    -1, /* 0x1.427b32daed80dp-20 */
    -1, /* 0x1.73e9e2254c98fp-20 */
    1, /* 0x1.8ec2e3d383d64p-19 */
    -1, /* 0x1.2226500110692p-18 */
    -1, /* 0x1.4c2c402a28e67p-18 */
    1, /* 0x1.59abdcf048a9bp-17 */
    1, /* 0x1.5ef5605b74f7p-16 */
    -1, /* 0x1.969df9fda1b9p-16 */
    1, /* 0x1.29320aa46f06ep-15 */
    -1, /* 0x1.6638f8a2185cbp-15 */
    1, /* 0x1.dd61dedc265f5p-15 */
    1, /* 0x1.41c4a3005149fp-14 */
    -1, /* 0x1.17485bab7607bp-13 */
    -1, /* 0x1.a0ea629bbf42p-13 */
    -1, /* 0x1.2878e662f31bcp-12 */
    -1, /* 0x1.cf88e9268bdcfp-12 */
    -1, /* 0x1.0d0c56d741955p-9 */
    1, /* 0x1.f0a606701a3e5p-9 */
    1, /* 0x1.730786d8eea35p-8 */
    1, /* 0x1.bfa7fc1b1e67cp-8 */
    1, /* 0x1.271adcddb2534p-7 */
    -1, /* 0x1.1f6bbb1fac712p-6 */
    -1, /* 0x1.20588554d6e24p-6 */
    -1, /* 0x1.5fa453f6bb40ep-6 */
    -1, /* 0x1.69a68ac73eaa2p-4 */
    -1, /* 0x1.da83e335e379ap-3 */
    1, /* 0x1.f608aa4e62781p-2 */
  };
// end_acospi

/* slow path, assumes |x| < 1 */
static double
slow_path (double x)
{
  double absx, y, h, l, u, v;
  union_t w;
  const double *p;
  absx = x > 0 ? x : -x;
  w.x = 1.0 + absx; /* 1 <= w.x <= 2 */
  /* Warning: w.x might be 2 for rounding up or nearest. */
  int i = (w.x == 2.0) ? 127 : (w.i[1] >> 13) & 127;
  if (i < 64) /* |x| < 0.5 */
  {
    // begin_acospi
    // deal with exceptional cases
    int a, b, c;
    for (a = 0, b = EXCEPTIONS; a + 1 != b;)
    {
      c = (a + b) / 2;
      if (exceptions[c][0] <= x)
        a = c;
      else
        b = c;
    }
    if (x == exceptions[a][0])
    {
      double h = exceptions[a][1];
      int8_t l = (h > 0) ? exceptions_rnd[a] : -exceptions_rnd[a];
      return h + h * 0x1p-54 * (double) l;
    }
    // end_acospi

    p = T2[i];
    y = absx - p[DEGREE+LARGE+1]; /* exact */
    h = p[DEGREE+LARGE];
    l = 0;
    for (int j = DEGREE - 1; j >= 0; j--)
    {
      if (j >= LARGE) /* use double only */
        h = p[LARGE+j] + h * y, l = 0;
      else /* use double-double */
      {
        /* multiply by y */
        mul2_1 (&u, &v, h, l, y);
        /* add coefficient of degree j */
        /* Warning: for i=0, the coefficients p[2*j] are tiny for j even,
           and fast_two_sum(h,l,a,b) requires |a| >= |b|, or at least
           exp(a) >= exp(b).  */
        if (__builtin_fabs (p[2*j]) >= __builtin_fabs (u))
          fast_two_sum (&h, &l, p[2*j], u);
        else
          fast_two_sum (&h, &l, u, p[2*j]);
        l += p[2*j+1] + v;
      }
    }
    /* acos(x) ~ pi/2 + (h + l) for x > 0, pi/2 - (h + l) for x < 0 */
    if (x < 0)
      h = -h, l = -l;
    fast_two_sum (&u, &v, pi_hi / 2, h);
    v += pi_lo / 2 + l;

    // acospi_begin
    // multiply by 1/pi
    d_mul (&u, &v, u, v, ONE_OVER_PIH, ONE_OVER_PIL);
    // acospi_end

    return u + v;
  }
  else /* 0.5 <= |x| < 1 */
  {
    double h1, l1;

    // exceptional cases
    if (x == 0x1.b32b7ac93ddefp-1)
      return 0x1.69c0e1dfbf177p-3 - 0x1.16beadd718bafp-108;
    if (x == 0x1.e55a7fa9a24c4p-1)
      return 0x1.a67c4d04a9236p-4 - 0x1.edc5fcb35e5e3p-110;

    h1 = 1.0 - absx; /* exact since |x| >= 0.5 */
    h1 = sqrt_dbl_dbl (h1, &l1);
    p = T2[i];
    y = absx - p[DEGREE+LARGE+1]; /* exact */
    h = p[DEGREE+LARGE];
    l = 0;
    for (int j = DEGREE - 1; j >= 0; j--)
    {
      if (j >= LARGE) /* use double only */
        h = p[LARGE+j] + h * y, l = 0;
      else /* use double-double */
      {
        /* multiply by y */
        mul2_1 (&u, &v, h, l, y);
        /* add coefficient of degree j */
        fast_two_sum (&h, &l, p[2*j], u);
        l += p[2*j+1] + v;
      }
    }
    /* acos(x) ~ (h1 + l1) * (h + l) */
    dekker (&u, &v, h1, h);
    v += l1 * h + h1 * l;
    if (x < 0)
    {
      fast_two_sum (&u, &l, pi_hi, -u);
      v = l + pi_lo - v;
    }

    // acospi_begin
    // multiply by 1/pi
    d_mul (&u, &v, u, v, ONE_OVER_PIH, ONE_OVER_PIL);
    // acospi_end

    return u + v;
  }
}

double
cr_acospi (double x)
{
  union_t u, v;
  int32_t k;
  u.x = x;
  u.i[1] &= 0x7fffffff; /* set sign bit to 0 */
  double absx = u.x;
  k = u.i[1];
  if (k < 0x3fe80000) { /* |x| < 0.75 */
    if (__builtin_expect (k == 0 && u.i[0] == 0, 0)) return 0.5;
    /* approximate acos(x) by pi/2 +/- p(x-xmid), where [0,0.75) is split
       into 192 sub-intervals */
    v.x = 1.0 + absx; /* 1 <= v.x < 2 */
    /* v.i[1] contains 20 significant bits in its low bits, we shift by 12
       to get the upper 8 (ignoring the implicit leading bit) */
    int i = (v.i[1] >> 12) & 255;
    if (__builtin_expect (i == 192, 0))
      i = 191;
    const double *p = T[i];
    double y = absx - p[7]; /* p[7] = xmid */
    double zh, zl;
    double yy = y * y;
    /* evaluate in parallel p[1] + p[2] * y and p[3] + p[4] * y, and
       p[5] + p[6] * y using Estrin's scheme */
    double p56 = p[5] + p[6] * y;
    double p34 = p[3] + p[4] * y;
    zh = p34 + yy * p56;
    zh = p[2] + y * zh;
    fast_two_sum (&zh, &zl, p[1], y * zh);
    double u, v;
    fast_two_sum (&u, &v, p[0], zh * y);
    v += zl * y;
    /* Special case for i=0, since we are obliged to use xmid=0 (so that
       x-xmid is exact) thus we can't use Gal's trick.  This costs about
       0.5 cycle in the average time (for both branches).  */
    if (i == 0)
      v += 0x4.6989e4b05fa3p-56;
    /* acos(x) ~ u + v for x > 0, pi - (u + v) for x < 0 */
    if (x < 0) /* acos(-x) = pi-acos(x) */
    {
      fast_two_sum (&u, &zl, pi_hi, -u);
      v = pi_lo + zl - v;
    }

    // acospi_begin
    /* We multiply the approximation u+v, with maximal error 0x1.5fp-61,
       by 1/pi. The maximal value of |u+v| is less than 2.42 (for x=-0.75).
       The maximal error is the sum of several terms:
       * 0x1.5fp-61 * (ONE_OVER_PIH + ONE_OVER_PIL) < 2^-62.196
       * (u+v)*|ONE_OVER_PIH+ONE_OVER_PIL-1/pi| < 2.42*2^-109.523 < 2^-108
       * the ignored term v*ONE_OVER_PIL in d_mul. The maximal observed value
         of v is 0x1.06d413839cafcp-51 for x=-0x1.6a01f2fb71p-1 (rndd),
         we conjecture |v| < 2^-50. Then |v*ONE_OVER_PIL| < 2^-105
       * the rounding error in d_mul. The d_mul call decomposes into:
         a_mul (u, lo1, u_in, ONE_OVER_PIH)
         lo2 = __builtin_fma (u_in, ONE_OVER_PIL, lo1)
         v = __builtin_fma (v_in, ONE_OVER_PIH, lo2)
         Since |u| <= acos(-0.75)/pi < 0.8 we have |lo1| <= ulp(0.8) <= 2^-53.
         Then since |u_in| <= 2.42, |lo2| <= |2.42*ONE_OVER_PIL|+2^-53
                                          < 2^-52.485
         Then |v| <= 2^-50+ONE_OVER_PIH*2^-52.485 < 2^-49.920.
         The rounding error is bounded by ulp(lo2)+ulp(v) <= 2^-105+2^-102
         < 2^-101.83.
       The total error is thus bounded by:
       2^-62.196 + 2^-108 + 2^-105 + 2^-101.83 < 2^-62.195.
       */
    d_mul (&u, &v, u, v, ONE_OVER_PIH, ONE_OVER_PIL);
    // acospi_end
    
    static const double err = 0x1.c0p-63; // acospi_specific, 2^-62.195 < 0x1.c0p-63
    double left  = u + (v - err), right = u + (v + err);
    if (left != right)
      return slow_path (x); /* hard to round case */
    return left;
  }
  /*--------------------------- 0.75 <= |x| < 1 ---------------------*/
  else
  if (k < 0x3ff00000) { /* |x| < 1 */
    /* approximate acos(x) by sqrt(1-x)*p(x-xmid) where p is a polynomial,
       and [0.75,1) is split into 64 sub-intervals */
    v.x = 1.0 + absx; /* 1 <= v.x <= 2 */
    /* The low 20 bits of v.i[1] are the upper bits (except the
       implicit leading bit) of the significand of 1+|x|.
       Warning: v.x might be 2 for rounding up or nearest. */
    int i = (v.x == 2.0) ? 255 : (v.i[1] & 0xff000) >> 12;
    const double *p = T[i];
    double y = absx - p[6]; /* exact (p[6] = xmid) */
    double h1, l1;
    h1 = 1.0 - absx; /* exact since |x| >= 0.5 */
    h1 = sqrt_dbl_dbl (h1, &l1);
    double zh, zl;
    /* use Estrin's scheme to evaluate p2 + p3*y + p4*y^2 + p5*y^3 */
    double yy = y * y;
    double p45 = p[4] + p[5] * y;
    double p23 = p[2] + p[3] * y;
    zh = p23 + p45 * yy;
    zh = p[1] + zh * y;
    fast_two_sum (&zh, &zl, p[0], zh * y);
    double l1zh = l1 * zh; /* compute earlier */
    double h1zl = h1 * zl;
    /* acos(x) ~ (h1 + l1) * (zh + zl) */
    double u, v;
    dekker (&u, &v, h1, zh);
    v += l1zh + h1zl;
    if (x < 0) /* acos(x) = pi - (u+v) */
    {
      fast_two_sum (&u, &zl, pi_hi, -u);
      /* acos(x) = u + zl + pi_lo - v */
      v = zl + pi_lo - v;
    }

    // acospi_begin
    /* Similar analysis as above.
       We multiply the approximation u+v, with maximal error 0x1.afp-65
       by 1/pi. The maximal value of |u+v| is pi (for x=-1).
       The maximal error is the sum of several terms:
       * 0x1.afp-65 * (ONE_OVER_PIH + ONE_OVER_PIL) < 2^-65.899
       * (u+v)*|ONE_OVER_PIH+ONE_OVER_PIL-1/pi| < pi*2^-109.523 < 2^-107
       * the ignored term v*ONE_OVER_PIL in d_mul. The maximal observed value
         of v is 0x1.4586d502c6913p-51 for x=-0x1.fffcc87dece8p-1 (rndd),
         we conjecture |v| < 2^-50. Then |v*ONE_OVER_PIL| < 2^-105
       * the rounding error in d_mul. The d_mul call decomposes into:
         a_mul (u, lo1, u_in, ONE_OVER_PIH)
         lo2 = __builtin_fma (u_in, ONE_OVER_PIL, lo1)
         v = __builtin_fma (v_in, ONE_OVER_PIH, lo2)
         Since |u| <= acospi(-1) < 1 we have |lo1| <= ulp(1-) <= 2^-53.
         Then since |u_in| <= pi, |lo2| <= |pi*ONE_OVER_PIL|+2^-53
                                          < 2^-52.361.
         Then |v| <= 2^-50+ONE_OVER_PIH*2^-52.361 < 2^-49.913.
         The rounding error is bounded by ulp(lo2)+ulp(v) <= 2^-105+2^-102
         < 2^-101.83.
       The total error is thus bounded by:
       2^-65.899 + 2^-107 + 2^-105 + 2^-101.83 < 2^-65.898
       */
    d_mul (&u, &v, u, v, ONE_OVER_PIH, ONE_OVER_PIL);
    // acospi_end

    static const double err = 0x1.13p-66; // acospi_specific, 2^-65.898 < 0x1.13p-66
    double left  = u + (v - err), right = u + (v + err);
    if (left != right)
      return slow_path (x); /* hard to round case */
    return left;
  }    /*   else  if (k < 0x3ff00000)    */

  /*---------------------------- |x|>=1 -----------------------*/
  else
    if (k==0x3ff00000 && u.i[0]==0) return (x>0) ? 0 : 1; // acospi_specific
  else
  if (k>0x7ff00000 || (k == 0x7ff00000 && u.i[0] != 0)) return x + x;
  else {
    u.i[1]=0x7ff00000;
    v.i[1]=0x7ff00000;
    u.i[0]=0;
    v.i[0]=0;
    return u.x/v.x;
  }
}

#ifndef SKIP_C_FUNC_REDEF
/* just to compile since glibc does not provide this function */
double acospi (double x)
{
  return acos (x) / M_PI;
}
#endif
