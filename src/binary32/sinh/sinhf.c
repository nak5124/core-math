/* Correctly-rounded hyperbolic sine function for binary32 value.

Copyright (c) 2022 Alexei Sibidanov.

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

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;

float cr_sinhf(float x){
  static const double c[] =
    {0x1.62e42fefa398bp-5, 0x1.ebfbdff84555ap-11, 0x1.c6b08d4ad86d3p-17,
     0x1.3b2ad1b1716a2p-23, 0x1.5d7472718ce9dp-30, 0x1.4a1d7f457ac56p-37};
  static const double tb[] =
    {0x1p+0, 0x1.0b5586cf9890fp+0, 0x1.172b83c7d517bp+0, 0x1.2387a6e756238p+0,
     0x1.306fe0a31b715p+0, 0x1.3dea64c123422p+0, 0x1.4bfdad5362a27p+0, 0x1.5ab07dd485429p+0,
     0x1.6a09e667f3bcdp+0, 0x1.7a11473eb0187p+0, 0x1.8ace5422aa0dbp+0, 0x1.9c49182a3f09p+0,
     0x1.ae89f995ad3adp+0, 0x1.c199bdd85529cp+0, 0x1.d5818dcfba487p+0, 0x1.ea4afa2a490dap+0};
  static const struct {union{float arg; uint32_t uarg;}; float rh, rl;} st[] = {
    {{0x1.250bfep-11}, 0x1.250bfep-11f, 0x1.fffffep-36f}
  };
  const double iln2h = 0x1.7154765p+0*16, iln2l = 0x1.5c17f0bbbe88p-31*16;
  b32u32_u t = {.f = x};
  double z = x;
  uint32_t ux = t.u, ex = (ux>>23)&0xff;
  if (__builtin_expect(ex>127+7, 0)){
    static const float ir[] = {__builtin_inff(),-__builtin_inff()};
    if(ex==0xff) {
      if(ux<<9) return x; // nan
      return ir[ux>>31]; // +-inf
    }
    static const float q[] = {0x1.fffffep127f,-0x1.fffffep127f};
    return q[ux>>31]*2.0f;
  }
  if (__builtin_expect(ex<127-3, 0)){
    if (__builtin_expect(ex<127-12, 0)){
      if (__builtin_expect(ex<127-25, 0))
	return __builtin_fmaf(x, __builtin_fabsf(x), x);
      float x2 = x*x;
      return __builtin_fmaf(x, 0x1.555556p-3f*x2, x);
    }
    if(__builtin_expect(st[0].uarg == (ux&(~0u>>1)), 0))
      return __builtin_copysignf(st[0].rh,x) + __builtin_copysignf(st[0].rl,x);
    static const double c[] = {0x1.5555555555555p-3, 0x1.11111111146e1p-7, 0x1.a01a00930dda6p-13, 0x1.71f92198aa6e9p-19};
    double z2 = z*z, z4 = z2*z2;
    return z + (z2*z)*((c[0] + z2*c[1]) + z4*(c[2] + z2*(c[3])));
  } else {
    double a = iln2h*z, ia = __builtin_floor(a), hp = (a - ia) + iln2l*z, hm = 1 - hp;
    long i = ia, jp = i&0xf, ep = i - jp, jm = ~i&0xf, em = ~i - jm;

    ep >>= 4;
    double sp = tb[jp];
    b64u64_u spu = {.u = (ep + 0x3feul)<<52};
    sp *= spu.f;
    double hp2 = hp*hp;
    double cp0 = c[0] + hp*c[1];
    double cp2 = c[2] + hp*c[3];
    double cp4 = c[4] + hp*c[5];
    cp0 += hp2*(cp2 + hp2*cp4);
    double wp = sp*hp;
    double rp = sp + wp*cp0;

    em >>= 4;
    double sm = tb[jm];
    b64u64_u smu = {.u = (em + 0x3feul)<<52};
    sm *= smu.f;
    double hm2 = hm*hm;
    double cm0 = c[0] + hm*c[1];
    double cm2 = c[2] + hm*c[3];
    double cm4 = c[4] + hm*c[5];
    cm0 += hm2*(cm2 + hm2*cm4);
    double wm = sm*hm;
    double rm = sm + wm*cm0;

    return rp - rm;
  }
}
