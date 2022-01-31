/* Correctly-rounded arctangent function of two binary32 values.

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

/*

  To compile:
  gcc -O3 -march=native -ffinite-math-only -frounding-math -fno-math-errno -W -Wall -c atan2f.c

*/

#include <stdint.h>

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;
typedef uint64_t u64;

static inline double muldd(double xh, double xl, double ch, double cl, double *l){
  double ahlh = ch*xl, alhh = cl*xh, ahhh = ch*xh, ahhl = __builtin_fma(ch, xh, -ahhh);
  ahhl += alhh + ahlh;
  ch = ahhh + ahhl;
  *l = (ahhh - ch) + ahhl;
  return ch;
}

static double polydd(double xh, double xl, int n, const double c[][2], double *l){
  int i = n-1;
  double ch = c[i][1], cl = c[i][0];
  while(--i>=0){
    ch = muldd(xh,xl,ch,cl,&cl);
    double th = ch + c[i][1], tl = (c[i][1] - th) + ch;
    ch = th;
    cl += tl + c[i][0];
  }
  *l = cl;
  return ch;
}

float cr_atan2f(float y, float x){
  static const double cn[] =
    {0x1p+0, 0x1.40e0698f94c35p+1, 0x1.248c5da347f0dp+1, 0x1.d873386572976p-1, 0x1.46fa40b20f1dp-3,
     0x1.33f5e041eed0fp-7, 0x1.546bbf28667c5p-14};
  static const double cd[] =
    {0x1p+0, 0x1.6b8b143a3f6dap+1, 0x1.8421201d18ed5p+1, 0x1.8221d086914ebp+0, 0x1.670657e3a07bap-2,
     0x1.0f4951fd1e72dp-5, 0x1.b3874b8798286p-11};
  static const double m[] = {0, 1};
  const double pi = 0x1.921fb54442d18p+1, pi2 = 0x1.921fb54442d18p+0, pi2l = 0x1.1a62633145c07p-54;
  static const double off[] = {0.0f, pi2, pi, pi2, -0.0f, -pi2, -pi, -pi2};
  static const double offl[] = {0.0f, pi2l, 2*pi2l, pi2l, -0.0f, -pi2l, -2*pi2l, -pi2l};
  static const double sgn[] = {1,-1};
  b32u32_u tx = {.f = x}, ty = {.f = y};
  uint32_t ux = tx.u, uy = ty.u, ax = ux&(~0u>>1), ay = uy&(~0u>>1);
  if(__builtin_expect(ay >= (0xff<<23)||ax >= (0xff<<23), 0)){
    if(ay > (0xff<<23)) return y; // nan
    if(ax > (0xff<<23)) return x; // nan
    uint32_t yinf = ay==(0xff<<23), xinf = ax==(0xff<<23);
    if(yinf&xinf){
      if(ux>>31)
	return 0x1.2d97c7f3321d2p+1*sgn[uy>>31];
      else
	return 0x1.921fb54442d18p-1*sgn[uy>>31];
    }
    if(xinf){
      if(ux>>31)
	return pi*sgn[uy>>31];
      else
	return 0.0f*sgn[uy>>31];
    }
    if(yinf){
      return pi2*sgn[uy>>31];
    }
  }
  if(__builtin_expect(ay==0, 0)){
    if(__builtin_expect(!(ay|ax),0)){
      uint32_t i = (uy>>31)*4 + (ux>>31)*2;
      if(ux>>31)
	return off[i] + offl[i];
      else
	return off[i];
    }
    if(!(ux>>31)) return 0.0f*sgn[uy>>31];
  }
  uint32_t gt = ay>ax, i = (uy>>31)*4 + (ux>>31)*2 + gt;
  
  double zx = x, zy = y;
  double z = (m[gt]*zx + m[1-gt]*zy)/(m[gt]*zy + m[1-gt]*zx);
  double z2 = z*z, z4 = z2*z2, z8 = z4*z4;
  double cn0 = cn[0] + z2*cn[1];
  double cn2 = cn[2] + z2*cn[3];
  double cn4 = cn[4] + z2*cn[5];
  double cn6 = cn[6];
  cn0 += z4*cn2;
  cn4 += z4*cn6;
  cn0 += z8*cn4;
  z *= sgn[gt];
  double cd0 = cd[0] + z2*cd[1];
  double cd2 = cd[2] + z2*cd[3];
  double cd4 = cd[4] + z2*cd[5];
  double cd6 = cd[6];
  cd0 += z4*cd2;
  cd4 += z4*cd6;
  cd0 += z8*cd4;
  double r = cn0/cd0;
  r = z*r + off[i];
  b64u64_u res = {.f = r};
  if(__builtin_expect(((res.u + 5)&0xfffffff) <= 10, 0)){
    double zh,zl;
    if(!gt){
      zh = zy/zx;
      zl = __builtin_fma(zh,-zx,zy)/zx;
    } else {
      zh = zx/zy;
      zl = __builtin_fma(zh,-zy,zx)/zy;
    }
    double z2l, z2h = muldd(zh,zl,zh,zl,&z2l);
    static const double c[30][2] =
      {{-0x1.bfdf6472p-82, 0x1p+0}, {-0x1.55522cf051bb7p-56, -0x1.5555555555555p-2},
       {-0x1.a13119a775722p-57, 0x1.999999999999ap-3}, {-0x1.80dd3b0eb53dap-57, -0x1.2492492492491p-3},
       {0x1.961c71122022fp-58, 0x1.c71c71c71c6a5p-4}, {0x1.d8873ae6474bfp-58, -0x1.745d1745d047ap-4},
       {0x1.47bd8f2f1877p-58, 0x1.3b13b13af39a1p-4}, {0x1.e7bda3f460852p-61, -0x1.1111110e9c5bbp-4},
       {0x1.0c07246705383p-59, 0x1.e1e1e199dd2adp-5}, {0x1.ae1ccf560cc5cp-60, -0x1.af28689a8395cp-5},
       {0x1.f3c877ef088b2p-60, 0x1.861844f9bb71fp-5}, {0x1.d686cb108e152p-59, -0x1.642bb7467eb59p-5},
       {-0x1.c8628a6b73a35p-61, 0x1.47a9501596294p-5}, {-0x1.c0c8a2f7773c8p-60, -0x1.2f50ec063dcc2p-5},
       {0x1.feb7021a2783cp-59, 0x1.1a1ba245d6116p-5}, {-0x1.95760e5ea6ff6p-60, -0x1.06f580c2b3b3cp-5},
       {-0x1.bac06658805ddp-62, 0x1.e8d3b0aa7e342p-6}, {-0x1.81be614231dep-61, -0x1.c0cba92af0035p-6},
       {-0x1.54e9ec905c7dcp-64, 0x1.90d85bf533d6p-6}, {-0x1.1dd5347f9d701p-63, -0x1.551ce2de13b14p-6},
       {0x1.5490a51372d33p-60, 0x1.0ddbdd787f62fp-6}, {-0x1.c766eb3ed3487p-62, -0x1.82b28ae9a24bbp-7},
       {0x1.e0282c6640316p-62, 0x1.e8c90da74be8dp-8}, {0x1.ec23b693ec582p-63, -0x1.094c35a3c5f4ap-8},
       {0x1.1bd48d253a2d1p-64, 0x1.e0ab2f3b33e79p-10}, {-0x1.c23a6acda6b24p-66, -0x1.5f500b1b46c96p-11},
       {0x1.2292eb52e1fd9p-67, 0x1.8c46d90303f2p-13}, {-0x1.1b39adf3ea87ap-69, -0x1.42a057ec505f1p-15},
       {0x1.a1f8c235de9f8p-72, 0x1.50986e7b11a12p-18}, {0x1.05f04cb8b6abfp-82, -0x1.514e4943fe90dp-22}};
    double pl, ph = polydd(z2h, z2l, 30, c, &pl);
    zh *= sgn[gt];
    zl *= sgn[gt];
    ph = muldd(zh,zl,ph,pl,&pl);
    double sh = ph + off[i], sl = ((off[i] - sh) + ph) + pl + offl[i];
    float rf = sh;
    double th = rf, dh = sh - th, tm = dh + sl;
    b64u64_u tth = {.f = th};
    if(th + th*0x1p-60 == th - th*0x1p-60){
      tth.u &= 0x7fful<<52;
      tth.u -= 24l<<52;
      if(__builtin_fabs(tm)>tth.f)
	tm *= 1.25;
      else
	tm *= 0.75;
    }
    r = th + tm;
  }
  return r;
}
