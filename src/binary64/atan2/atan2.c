/* Correctly-rounded atan2 function for two binary64 values.

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

#define TRACEY -0x1.49343d4e26bb6p+53
#define TRACEX 0x1.d786165a1b544p+51

#include <stdio.h>
#include <stdint.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union { double f; uint64_t u; } d64u64;
#include "tint.h"

// return non-zero if u (with sign bit cleared) encodes Inf or NaN
static inline int inf_or_nan (uint64_t u)
{
  return (u >> 52) == 0x7ff;
}

#define MASK 0x7ffffffffffffffful

static inline int is_nan (uint64_t u)
{
  u = u & MASK;
  uint64_t e = u >> 52;
  return e == 0x7ff && u != (e << 52);
}

static inline int is_inf (uint64_t u)
{
  return (u & MASK) == (0x7fful << 52);
}

// PI_H+PI_L approximates pi with error bounded by 2^-108.041
#define PI_H 0x1.921fb54442d18p+1
#define	PI_L 0x1.1a62633145c07p-53
// PI_OVER2_H+PI_OVER2_L approximates pi/2 with error bounded by 2^-109.041
#define PI_OVER2_H 0x1.921fb54442d18p+0
#define	PI_OVER2_L 0x1.1a62633145c07p-54
// PI_OVER4_H+PI_OVER4_L approximates pi/4 with error bounded by 2^-110.041
#define PI_OVER4_H 0x1.921fb54442d18p-1
#define	PI_OVER4_L 0x1.1a62633145c07p-55

/* The following is a degree-15 polynomial with odd coefficients
   approximating atan(x) on [0,2^-11.2] with maximal relative error
   2^-192.031 (cf atan2_small.sollya). */
static const tint_t Psmall[] = {
  // degree 1: 1
  {.h=0x8000000000000000, .m=0x0, .l=0x0, .ex=1, .sgn=0},
  // degree 3: -0x1.5555555555555555555555555555555555555555569de77ep-2
  {.h=0xaaaaaaaaaaaaaaaa, .m=0xaaaaaaaaaaaaaaaa, .l=0xaaaaaaaaab4ef3bf, .ex=-1, .sgn=1},
  // degree 5: 0x1.9999999999999999999999999999999999be8f1d48d2355cp-3
  {.h=0xcccccccccccccccc, .m=0xcccccccccccccccc, .l=0xccdf478ea4691aae, .ex=-2, .sgn=0},
  // degree 7: -0x1.2492492492492492492492492492c75ep-3
  {.h=0x9249249249249249, .m=0x24924924924963af, .l=0x0, .ex=-2, .sgn=1},
  // degree 9: 0x1.c71c71c71c71c71c71c71dbb1af83beap-4
  {.h=0xe38e38e38e38e38e, .m=0x38e38edd8d7c1df5, .l=0x0, .ex=-3, .sgn=0},
  // degree 11: -0x1.745d1745d1745d18p-4
  {.h=0xba2e8ba2e8ba2e8c, .m=0x0, .l=0x0, .ex=-3, .sgn=1},
  // degree 13: 0x1.3b13b13b13d919bap-4
  {.h=0x9d89d89d89ec8cdd, .m=0x0, .l=0x0, .ex=-3, .sgn=0},
  // degree 15: -0x1.1111103a0ee2178ep-4
  {.h=0x8888881d07710bc7, .m=0x0, .l=0x0, .ex=-3, .sgn=1},
};

// case |y/x| < 2^-11.2 or 2^11.2 < |y/x|
static double
atan2_accurate_small_or_large (double y, double x)
{
  if (y == TRACEY && x == TRACEX) printf ("atan2_accurate_small_or_large: y=%la x=%la\n", y, x);
  /* first check when t=y/x is small and exact and x > 0, since for
     |t| <= 0x1.d12ed0af1a27fp-27, atan(t) rounds to t (to nearest) */
  double t = y / x;
  double corr = __builtin_fma (t, x, -y);
  if (corr == 0 && x > 0) // t is exact
    if (__builtin_fabs (t) <= 0x1.d12ed0af1a27fp-27)
      return __builtin_fma (t, -0x1p-54, t);

  int inv = __builtin_fabs (t) > 1;

  tint_t z[1], z2[1], p[1];
  if (inv == 0)
    div_tint_d (z, y, x);
  else
    div_tint_d (z, x, y);
  mul_tint (z2, z, z);
  cp_tint (p, Psmall+7); // degree 15
  for (int i = 6; i >= 0; i--)
  {
    mul_tint (p, p, z2);
    add_tint (p, p, Psmall+i);
  }
  // multiply by z
  mul_tint (p, p, z);
  if (inv)
  {
    if (z->sgn == 0) { // z > 0
      p->sgn = 1 - p->sgn;
      add_tint (p, &PI2, p);
    }
    else {
      add_tint (p, &PI2, p);
      p->sgn = 1 - p->sgn;
    }
  }
  // if x is negative we go to the opposite quadrant
  if (x < 0) {
    if (p->sgn == 0) { // 1st quadrant -> 3rd quadrant (subtract pi)
      p->sgn = 1;
      add_tint (p, &PI, p);
      p->sgn = 1;
    }
    else // 4th quadrant -> 2nd quadrant (add pi)
      add_tint (p, &PI, p);
  }
  return tint_tod (p);
}

/* The following polynomials P[] and Q[] are a rational approximation
   of atan(z) on [0.000425,1.0] generated by rminimax using the following
   command:
   ./ratapprox --function="atan(x)" --dom=[0.000425,1.0] --type=[29,29] --output=x1p.sollya --log --numF=[192] --denF=[192] --prec=1024
   The relative approximation error output by rminimax is:
   fpminimax approximation error                    = 2.1867e-58
   The content of the output x1p.sollya file is the following
   (Numerator should match P[], and Denominator should match Q[]):
Numerator = [|
-0x7.3d7a1dfb49d278ccb715e4114e96c6666b4446f6ed3474b8p-212,
0x1.94e9ffa36f41c07101c28621e2d24a6f1fb143956d7e7bbp-12,
0x9.5436fa854afbd456b8bbf4c71be08631ef581a18ac2797bp-12,
0x2.587f99d02aec2e24d2691aef5342878cac931658f167cf94p-8,
0x6.eef4426fdf55888b29f68ac4318c52f36ac8fe5d0828d24p-8,
0x1.0ce608c67ec76ca2552201613daef4d663bff51218de5f3ap-4,
0x2.26636bd1256c7d5d712cd447f6dcf522250a7bc4beb505ap-4,
0x3.d466647655c9dd8fd3cd2d886f9384528447bd2322f534b8p-4,
0x6.05fce37ce6622fa0980922ee2dcb797ba4bcd33351775718p-4,
0x8.785ab1f453705b902d155327664299e0fbaa9738f080e8p-4,
0xa.bb0b512681e73baeef6cf2783869cff93af2bba88207f6ep-4,
0xc.522372318d71f95282e6a19c16d337a00b88fdef9a34681p-4,
0xc.dee76309cb8a742eeed33eb46fba9b76184fef9199d1e9cp-4,
0xc.43a9b950d684e1ee645ec6e00a64e45b6a6cf4bfa91eba2p-4,
0xa.ac9b9a7a079078fefe3566c2f54b42d48bd2827334e6992p-4,
0x8.7cff89e9b24a17da7629d3b930113b85b3920b492a91ec9p-4,
0x6.29a8825dfa98afc0d6bcb914fd73dc7d65d791e364ea8cp-4,
0x4.13f4ce2f501b02cfc16940beacc255fae26acb10fe717408p-4,
0x2.73d0b89743d45bc7bca743cd22120307a840c2cafea455ep-4,
0x1.55bb9ea2ad2516ee842ce08950f2bb79d5881a388b2a6826p-4,
0xa.76d64812b753dbc71326ff30c5599fd5130e86af59d52d2p-8,
0x4.949539ab7bb896f8b3182c2e7bc1734ab8855685eceaf0ap-8,
0x1.c60fb3ead4fb41d4a84d4579865584a965ebf60cab52b158p-8,
0x9.9814bf383257da821a8fe242024748639bbb8ce304df07ap-12,
0x2.c7f21840b0d3c1a7f129cb96384cb4207e8c5691a00f973p-12,
0xa.c7cc3bf4146a522a4fd4e4b0930343a776a27cbc0516c9fp-16,
0x2.19fa5637351f3c0e014ecb3de333302bb27329721db6ee88p-16,
0x4.f5690187e14cbcafce2dbb0c62fcb6ddbe718943f354e838p-20,
0x7.e4c8db239c5a930539f8a20cfae9596525de3c1e6825ee2p-24,
0x6.61a1007cfa0b49905761b86a46f2626e4cfc098d3942713p-28|];
Denominator = [|
0x1.94e9ffa36f41c07101c28621e2d24a6f1fb143956d7e7b34p-12,
0x9.5436fa854afbd456b8bbf4c71be08631ef581a18ac2a7a4p-12,
0x2.60ef2478e7e8377c82727dbab2a6e919a7e6c7c1632ec43p-8,
0x7.20b567a7fb901cf7a3251fddaccbab14751ad43e52fac888p-8,
0x1.1944b344d1db581b3cfb815c6ff622d5e76c4b10fdb0b7b8p-4,
0x2.4a898eb6cef07ecd777f0e497e5fe453771800d544b4fb2p-4,
0x4.2ac536fb4141a5c3739836d52e5c908c34f96baea5a746b8p-4,
0x6.b4062e3066b07151635c322532a937cb7392bdeeb78e13c8p-4,
0x9.a8f3af8a77de0d14007f1f31ef89fd5a01248168e49f6e7p-4,
0xc.9100122c50ab76541f598d816d2c08c62e7a68a8a299eeap-4,
0xe.d932ad9ecaf82250674a27d7895765843db41d03bee0276p-4,
0x1p+0,
0xf.c49ed94a7267647763d93f0eaaa6a5adafaa4b63a84cdbfp-4,
0xe.3c2aa93d24dbbbf29c439b479305c86a6191280d431d255p-4,
0xb.c7cb990cce68d0692b52316098fa98e441344660792ac1p-4,
0x8.ef40e19a39b3eb5c76b49475360dad2a1709b63693dd9adp-4,
0x6.33f0be18ca51ef2556d1528c962b8e314b5d93fb1b12f0ap-4,
0x3.ef094bbd8728d9bab1cdb7369fbbf41ab27032886ea775cp-4,
0x2.458ad62b206610f532a7eb0669645297bf031a192a3ba9fp-4,
0x1.307a2c00697d7c2603e7acc134dbd2a55c91d5c441fac19ap-4,
0x8.fb28045098521ebc56c497d8ee5866b732abc4da0841659p-8,
0x3.caaa9d4e5c4e992d01834ae1b38543034c65c731463c5e9cp-8,
0x1.6aed4c0b712ef7381b0839dfec431450dca29356524e3b38p-8,
0x7.6945e7b98183547fd9f356fd30e43789b5ed5afc92fff97p-12,
0x2.13ee70fec48e5880d4dd182a6548a6b6e6b9c686f1485ed8p-12,
0x7.cba649b60294f0567ef6503a3a315ebd3e725a869e8f8d2p-16,
0x1.78b0b03c8604729d199e8beeb116acd983833c47417d2d58p-16,
0x3.5cf714cc1eadd7a74fbee173119e521697f9d8ed5b28d71p-20,
0x5.2fda0b1c5cc46dc139583528da5c10b9c572ba82134d3dfp-24,
0x4.10001e9d0cd93e9de87d4466d548a0363ccd56c5e00d1e6p-28|];
*/
static const tint_t P[30] = {
  {.h=0xe7af43bf693a4f19, .m=0x96e2bc8229d2d8cc, .l=0xcd6888dedda68e97, .ex=-209, .sgn=1},
   {.h=0xca74ffd1b7a0e038, .m=0x80e14310f1692537, .l=0x8fd8a1cab6bf3dd8, .ex=-11, .sgn=0},
   {.h=0x95436fa854afbd45, .m=0x6b8bbf4c71be0863, .l=0x1ef581a18ac2797b, .ex=-8, .sgn=0},
   {.h=0x961fe6740abb0b89, .m=0x349a46bbd4d0a1e3, .l=0x2b24c5963c59f3e5, .ex=-6, .sgn=0},
   {.h=0xddde884dfbeab111, .m=0x653ed15886318a5e, .l=0x6d591fcba1051a48, .ex=-5, .sgn=0},
   {.h=0x867304633f63b651, .m=0x2a9100b09ed77a6b, .l=0x31dffa890c6f2f9d, .ex=-3, .sgn=0},
   {.h=0x8998daf4495b1f57, .m=0x5c4b3511fdb73d48, .l=0x89429ef12fad4168, .ex=-2, .sgn=0},
   {.h=0xf519991d95727763, .m=0xf4f34b621be4e114, .l=0xa111ef48c8bd4d2e, .ex=-2, .sgn=0},
   {.h=0xc0bf9c6f9ccc45f4, .m=0x1301245dc5b96f2f, .l=0x74979a666a2eeae3, .ex=-1, .sgn=0},
   {.h=0x8785ab1f453705b9, .m=0x2d155327664299e, .l=0xfbaa9738f080e80, .ex=0, .sgn=0},
   {.h=0xabb0b512681e73ba, .m=0xeef6cf2783869cff, .l=0x93af2bba88207f6e, .ex=0, .sgn=0},
   {.h=0xc522372318d71f95, .m=0x282e6a19c16d337a, .l=0xb88fdef9a34681, .ex=0, .sgn=0},
   {.h=0xcdee76309cb8a742, .m=0xeeed33eb46fba9b7, .l=0x6184fef9199d1e9c, .ex=0, .sgn=0},
   {.h=0xc43a9b950d684e1e, .m=0xe645ec6e00a64e45, .l=0xb6a6cf4bfa91eba2, .ex=0, .sgn=0},
   {.h=0xaac9b9a7a079078f, .m=0xefe3566c2f54b42d, .l=0x48bd2827334e6992, .ex=0, .sgn=0},
   {.h=0x87cff89e9b24a17d, .m=0xa7629d3b930113b8, .l=0x5b3920b492a91ec9, .ex=0, .sgn=0},
   {.h=0xc535104bbf5315f8, .m=0x1ad797229fae7b8f, .l=0xacbaf23c6c9d5180, .ex=-1, .sgn=0},
   {.h=0x827e99c5ea036059, .m=0xf82d2817d5984abf, .l=0x5c4d59621fce2e81, .ex=-1, .sgn=0},
   {.h=0x9cf42e25d0f516f1, .m=0xef29d0f3488480c1, .l=0xea1030b2bfa91578, .ex=-2, .sgn=0},
   {.h=0xaaddcf5156928b77, .m=0x42167044a8795dbc, .l=0xeac40d1c45953413, .ex=-3, .sgn=0},
   {.h=0xa76d64812b753dbc, .m=0x71326ff30c5599fd, .l=0x5130e86af59d52d2, .ex=-4, .sgn=0},
   {.h=0x9292a7356f7712df, .m=0x16630585cf782e69, .l=0x5710aad0bd9d5e14, .ex=-5, .sgn=0},
   {.h=0xe307d9f56a7da0ea, .m=0x5426a2bcc32ac254, .l=0xb2f5fb0655a958ac, .ex=-7, .sgn=0},
   {.h=0x99814bf383257da8, .m=0x21a8fe2420247486, .l=0x39bbb8ce304df07a, .ex=-8, .sgn=0},
   {.h=0xb1fc86102c34f069, .m=0xfc4a72e58e132d08, .l=0x1fa315a46803e5cc, .ex=-10, .sgn=0},
   {.h=0xac7cc3bf4146a522, .m=0xa4fd4e4b0930343a, .l=0x776a27cbc0516c9f, .ex=-12, .sgn=0},
   {.h=0x867e958dcd47cf03, .m=0x8053b2cf78cccc0a, .l=0xec9cca5c876dbba2, .ex=-14, .sgn=0},
   {.h=0x9ead2030fc299795, .m=0xf9c5b7618c5f96db, .l=0xb7ce31287e6a9d07, .ex=-17, .sgn=0},
   {.h=0xfc991b64738b5260, .m=0xa73f14419f5d2b2c, .l=0xa4bbc783cd04bdc4, .ex=-21, .sgn=0},
   {.h=0xcc34200f9f416932, .m=0xaec370d48de4c4d, .l=0xc99f8131a7284e26, .ex=-25, .sgn=0},
};

static const tint_t Q[30] = {
   {.h=0xca74ffd1b7a0e038, .m=0x80e14310f1692537, .l=0x8fd8a1cab6bf3d9a, .ex=-11, .sgn=0},
   {.h=0x95436fa854afbd45, .m=0x6b8bbf4c71be0863, .l=0x1ef581a18ac2a7a4, .ex=-8, .sgn=0},
   {.h=0x983bc91e39fa0ddf, .m=0x209c9f6eaca9ba46, .l=0x69f9b1f058cbb10c, .ex=-6, .sgn=0},
   {.h=0xe416acf4ff72039e, .m=0xf464a3fbb5997562, .l=0x8ea35a87ca5f5911, .ex=-5, .sgn=0},
   {.h=0x8ca259a268edac0d, .m=0x9e7dc0ae37fb116a, .l=0xf3b625887ed85bdc, .ex=-3, .sgn=0},
   {.h=0x92a263adb3bc1fb3, .m=0x5ddfc3925f97f914, .l=0xddc60035512d3ec8, .ex=-2, .sgn=0},
   {.h=0x8558a6df682834b8, .m=0x6e7306daa5cb9211, .l=0x869f2d75d4b4e8d7, .ex=-1, .sgn=0},
   {.h=0xd680c5c60cd60e2a, .m=0x2c6b8644a65526f9, .l=0x6e7257bdd6f1c279, .ex=-1, .sgn=0},
   {.h=0x9a8f3af8a77de0d1, .m=0x4007f1f31ef89fd5, .l=0xa01248168e49f6e7, .ex=0, .sgn=0},
   {.h=0xc9100122c50ab765, .m=0x41f598d816d2c08c, .l=0x62e7a68a8a299eea, .ex=0, .sgn=0},
   {.h=0xed932ad9ecaf8225, .m=0x674a27d78957658, .l=0x43db41d03bee0276, .ex=0, .sgn=0},
   {.h=0x8000000000000000, .m=0x0, .l=0x0, .ex=1, .sgn=0},
   {.h=0xfc49ed94a7267647, .m=0x763d93f0eaaa6a5a, .l=0xdafaa4b63a84cdbf, .ex=0, .sgn=0},
   {.h=0xe3c2aa93d24dbbbf, .m=0x29c439b479305c86, .l=0xa6191280d431d255, .ex=0, .sgn=0},
   {.h=0xbc7cb990cce68d06, .m=0x92b52316098fa98e, .l=0x441344660792ac10, .ex=0, .sgn=0},
   {.h=0x8ef40e19a39b3eb5, .m=0xc76b49475360dad2, .l=0xa1709b63693dd9ad, .ex=0, .sgn=0},
   {.h=0xc67e17c3194a3de4, .m=0xaada2a5192c571c6, .l=0x296bb27f63625e14, .ex=-1, .sgn=0},
   {.h=0xfbc252ef61ca366e, .m=0xac736dcda7eefd06, .l=0xac9c0ca21ba9dd70, .ex=-2, .sgn=0},
   {.h=0x9162b58ac819843d, .m=0x4ca9fac19a5914a5, .l=0xefc0c6864a8eea7c, .ex=-2, .sgn=0},
   {.h=0x983d160034bebe13, .m=0x1f3d6609a6de952, .l=0xae48eae220fd60cd, .ex=-3, .sgn=0},
   {.h=0x8fb28045098521eb, .m=0xc56c497d8ee5866b, .l=0x732abc4da0841659, .ex=-4, .sgn=0},
   {.h=0xf2aaa7539713a64b, .m=0x4060d2b86ce150c0, .l=0xd31971cc518f17a7, .ex=-6, .sgn=0},
   {.h=0xb576a605b8977b9c, .m=0xd841ceff6218a28, .l=0x6e5149ab29271d9c, .ex=-7, .sgn=0},
   {.h=0xed28bcf730306a8f, .m=0xfb3e6adfa61c86f1, .l=0x36bdab5f925fff2e, .ex=-9, .sgn=0},
   {.h=0x84fb9c3fb1239620, .m=0x3537460a995229ad, .l=0xb9ae71a1bc5217b6, .ex=-10, .sgn=0},
   {.h=0xf974c936c0529e0a, .m=0xcfdeca0747462bd7, .l=0xa7ce4b50d3d1f1a4, .ex=-13, .sgn=0},
   {.h=0xbc58581e4302394e, .m=0x8ccf45f7588b566c, .l=0xc1c19e23a0be96ac, .ex=-15, .sgn=0},
   {.h=0xd73dc53307ab75e9, .m=0xd3efb85cc4679485, .l=0xa5fe763b56ca35c4, .ex=-18, .sgn=0},
   {.h=0xa5fb41638b988db8, .m=0x272b06a51b4b8217, .l=0x38ae57504269a7be, .ex=-21, .sgn=0},
   {.h=0x820003d3a19b27d3, .m=0xbd0fa88cdaa91406, .l=0xc799aad8bc01a3cc, .ex=-25, .sgn=0},
};

// use a type [29,29] rational approximation of atan(z) for 0.000425 <= z <= 1
static double
atan2_accurate_rminimax (double y, double x)
{
  //if (y == TRACEY && x == TRACEX) printf ("atan2_accurate_rminimax: y=%la x=%la\n", y, x);
  int inv = __builtin_fabs (y) > __builtin_fabs (x);
  tint_t z[1], p[1], q[1];
  if (inv)
    div_tint_d (z, x, y);
  else
    div_tint_d (z, y, x);
  // the rational approximation is only for z > 0, it is not antisymmetric
  int sz = z->sgn;
  z->sgn = 0;
  // if (y == TRACEY && x == TRACEX) { printf ("z="); print_tint (z); }
  cp_tint (p, P + 29);
  cp_tint (q, Q + 29);
  //if (y == TRACEY && x == TRACEX) { printf ("p29="); print_tint (p); }
  //if (y == TRACEY && x == TRACEX) { printf ("q="); print_tint (q); }
  for (int i = 28; i >= 0; i--)
  {
    mul_tint (p, p, z);
    mul_tint (q, q, z);
    //if (i==0 && y == TRACEY && x == TRACEX) { printf ("p%d=", i); print_tint (p); }
    add_tint (p, p, P + i);
    add_tint (q, q, Q + i);
    //if (i==0 && y == TRACEY && x == TRACEX) { printf ("p%d=", i); print_tint (p); }
  }
  // if (y == TRACEY && x == TRACEX) { printf ("p="); print_tint (p); }
  // if (y == TRACEY && x == TRACEX) { printf ("q="); print_tint (q); }
  // divide p by q
  div_tint (z, p, q);
  z->sgn = sz; // restore sign
  // if (y == TRACEY && x == TRACEX) { printf ("z1="); print_tint (z); }
  /* Now z approximates atan(y/x) for inv=0, and atan(x/y) for inv=1,
     with -pi/4 < z < pi/4.
  */
  if (inv)
  {
    // if (y == TRACEY && x == TRACEX) printf ("+/-pi/2-atan(x/y)\n");
    // if x/y > 0 thus atan(x/y) > 0 we apply pi/2 - atan(x/y)
    // if x/y < 0 thus atan(x/y) < 0 we apply -pi/2 - atan(x/y)
    if (z->sgn == 0) { // 0 < atan(x/y) < pi/4
      z->sgn = 1;
      add_tint (z, &PI2, z);
      // now pi/4 < z < pi/2
    }
    else // -pi/4 < atan(x/y) < 0
    {
      add_tint (z, &PI2, z);
      z->sgn = 1;
      // now -pi/2 < z < -pi/4
    }
    // if (y == TRACEY && x == TRACEX) { printf ("z2="); print_tint (z); }
  }
  // if x is negative we go to the opposite quadrant
  if (x < 0) {
    if (z->sgn == 0) { // 1st quadrant -> 3rd quadrant (subtract pi)
      z->sgn = 1;
      add_tint (z, &PI, z);
      z->sgn = 1;
    }
    else // 4th quadrant -> 2nd quadrant (add pi)
      add_tint (z, &PI, z);
  }
  return tint_tod (z);
}

// accurate path, assumes both y and x are neither NaN, nor +/-Inf, nor +/-0
static double
atan2_accurate (double y, double x)
{
  double z = y / x;

#define Z0 0x1.bdb8cdadbe12p-12
#define Z1 0x1.2611186bae675p+11
  if (__builtin_fabs (z) <= Z0 || Z1 <= __builtin_fabs (z))
    // |z| < 2^-11.2 or 2^11.2 < |z|
    return atan2_accurate_small_or_large (y, x);

  return atan2_accurate_rminimax (y, x);
}

// atan(y/x)
double cr_atan2 (double y, double x)
{
  d64u64 uy = {.f = y}, ux = {.f = x};
  uint64_t ay = uy.u & MASK, ax = ux.u & MASK;

  if (__builtin_expect (inf_or_nan (ay) || inf_or_nan (ax), 0))
  {
    if (is_nan (ay) || is_nan (ax))
      return y + x; // if y or x is sNaN, returns qNaN are raises invalid
    // Now neither y nor x is NaN, but at least one is +Inf or -Inf
    if (is_inf (ay) && is_inf (ax)) // both y and x are +/-Inf
    {
      // atan2 (+/-Inf,-Inf) = +/-3pi/4
      if (x < 0)
        return (y > 0) ? 3 * PI_OVER4_H + 3 * PI_OVER4_L
          : -3 * PI_OVER4_H - 3 * PI_OVER4_L;
      // atan2 (+/-Inf,+Inf) = +/-pi/4
      return (y > 0) ? PI_OVER4_H + PI_OVER4_L : -PI_OVER4_H - PI_OVER4_L;
    }
    // now only one of y and x is +/-Inf
    if (is_inf (ax))
    {
      if (x < 0)
        return (uy.u >> 63) ? -PI_H - PI_L : PI_H + PI_L;
      // atan2(+/-0,x) = +/-0 for x > 0
      // atan2(+/-y,+Inf) = +/-0 for finite y>0
      return __builtin_copysign (0, y);
    }
    // now y = +/-Inf
    // atan2(+/-Inf,x) = +/-pi/2 for finite x
    return (y > 0) ? PI_OVER2_H + PI_OVER2_L : -PI_OVER2_H - PI_OVER2_L;
  }

  if (__builtin_expect (y == 0 || x == 0, 0))
  {
    if (y == 0 && x == 0)
    {
      if (ux.u == 0) // atan2(+/-0, +0) = +/-0
        return y;
      // atan2(+/-0, +0) = +/-pi
      return (uy.u == 0) ? PI_H + PI_L : -PI_H - PI_L;
    }
    // only one of y and x is zero
    if (y == 0)
    {
      // atan2(+/-0,x) = +/-0 for x>0
      if (x > 0) return y;
      // atan2(+/-0,x) = +/-pi for x<0
      return (uy.u == 0) ? PI_H + PI_L : -PI_H - PI_L;
    }
    // now only x is zero
    // atan2(y,+/-0) = -pi/2 for y<0
    // atan2(y,+/-0) = +pi/2 for y>0
    return (y > 0) ? PI_OVER2_H + PI_OVER2_L : -PI_OVER2_H - PI_OVER2_L;
  }

  // now both y and x are neither NaN, nor +/-Inf, nor +/-0

  return atan2_accurate (y, x);

  return 0.0;
}
