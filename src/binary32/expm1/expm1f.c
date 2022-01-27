/* Correctly-rounded natural exponent function biased by 1 for binary32 value.

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

typedef union {float f; unsigned u;} b32u32_u;
typedef union {double f; unsigned long u;} b64u64_u;

float cr_expm1f(float x){
  static const double c[] =
    {0x1.62e42fefa398bp-5, 0x1.ebfbdff84555ap-11, 0x1.c6b08d4ad86d3p-17,
     0x1.3b2ad1b1716a2p-23, 0x1.5d7472718ce9dp-30, 0x1.4a1d7f457ac56p-37};
  static const double tb[] =
    {0x1p+0, 0x1.0b5586cf9890fp+0, 0x1.172b83c7d517bp+0, 0x1.2387a6e756238p+0,
     0x1.306fe0a31b715p+0, 0x1.3dea64c123422p+0, 0x1.4bfdad5362a27p+0, 0x1.5ab07dd485429p+0,
     0x1.6a09e667f3bcdp+0, 0x1.7a11473eb0187p+0, 0x1.8ace5422aa0dbp+0, 0x1.9c49182a3f09p+0,
     0x1.ae89f995ad3adp+0, 0x1.c199bdd85529cp+0, 0x1.d5818dcfba487p+0, 0x1.ea4afa2a490dap+0};
  static const float q[][2] = {{0x1.fffffep127f, 0x1.fffffep127f}, {-1.0f, 0x1p-26f}};
  const double iln2h = 0x1.7154765p+0*16, iln2l = 0x1.5c17f0bbbe88p-31*16;
  b32u32_u t = {.f = x};
  double z = x;
  unsigned ux = t.u, ax = ux&(~0u>>1);
  if(__builtin_expect(ux>0xc18aa123u, 0)){ // x < -17.32
    if(ax>(0xffu<<23)) return x; // nan
    return q[1][0] + q[1][1];
  } else if(__builtin_expect(ax>0x42b17218u, 0)){  // x > 88.72
    if(ax>(0xffu<<23)) return x; // nan
    return q[0][0] + q[0][1];
  } else if (__builtin_expect(ax<0x3e000000u, 1)){ // x < 0.125
    if (__builtin_expect(ax<0x32000000u, 0)) // x < 2^-25
      return __builtin_fmaf(x,x,x);
    static const double p[] =
      {0x1.ffffffffffff6p-2, 0x1.5555555555572p-3, 0x1.5555555566a8fp-5, 0x1.11111110f18aep-7,
       0x1.6c16bf78e5645p-10, 0x1.a01a03fd7c6cdp-13, 0x1.a0439d78f6d66p-16, 0x1.71de38ef84d8cp-19};
    double z2 = z*z, z4 = z2*z2;
    double c0 = p[0] + z*p[1];
    double c2 = p[2] + z*p[3];
    double c4 = p[4] + z*p[5];
    double c6 = p[6] + z*p[7];
    c0 += z2*c2;
    c4 += z2*c6;
    c0 += z4*c4;
    return z + z2*c0;
  } else {
    double a = iln2h*z, ia = __builtin_floor(a), h = (a - ia) + iln2l*z;
    long i = ia, j = i&0xf, e = i - j;
    e >>= 4;
    double s = tb[j];
    b64u64_u su = {.u = (e + 0x3fful)<<52};
    s *= su.f;
    double h2 = h*h;
    double c0 = c[0] + h*c[1];
    double c2 = c[2] + h*c[3];
    double c4 = c[4] + h*c[5];
    c0 += h2*(c2 + h2*c4);
    double w = s*h;
    return (s-1.0) + w*c0;
  }
}
