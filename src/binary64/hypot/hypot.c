/* Correctly-rounded Euclidean distance function (hypot) for binary64 values.

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
#include <errno.h>
#include <x86intrin.h>

typedef uint64_t u64;
typedef int64_t i64;
typedef unsigned __int128 u128;
typedef union {double f; u64 u;} b64u64_u;

static inline double fasttwosum(double x, double y, double *e){
  double s = x + y;
  double z = s - x;
  *e = y - z;
  return s;
}

static double __attribute__((noinline)) as_hypot_denorm(u64 a, u64 b){
  double op = 1.0 + 0x1p-54, om = 1.0 - 0x1p-54;
  double af = (i64)a, bf = (i64)b;
  a <<= 1;
  b <<= 1;
  u64 rm = __builtin_sqrt(af*af + bf*bf);
  u64 tm = rm << 1;
  i64 D = a*a - tm*tm + b*b, sD = D>>63, um = (rm^sD) - sD, drm = sD + 1;
  i64 dd = (um<<3) + 4, pD;
  rm -= drm;
  drm += sD;
  do {
    pD = D;
    rm += drm;
    D -= dd;
    dd += 8;
  } while(__builtin_expect((D^pD)>0, 0));
  pD = (sD&D)+(~sD&pD);
  if(__builtin_expect(pD != 0, 1)){
    if(__builtin_expect(op == om, 1)){
      i64 sum = pD - 4*rm - 1;
      if(__builtin_expect(sum, 1))
	rm += (sum>>63) + 1;
      else
	rm += rm&1;
    } else {
      rm += op>1.0;
    }
  }
  b64u64_u xi = {.u = rm};
  return xi.f;
}

static double  __attribute__((noinline)) as_hypot_hard(double x, double y){
  double op = 1.0 + 0x1p-54, om = 1.0 - 0x1p-54;
  b64u64_u xi = {.f = x}, yi = {.f = y};
  u64 bm = (xi.u&(~0ul>>12))|1l<<52;
  u64 lm = (yi.u&(~0ul>>12))|1l<<52;
  int be = xi.u>>52;
  int le = yi.u>>52;
  b64u64_u ri = {.f = __builtin_sqrt(x*x + y*y)};
  const int bs = 3;
  u64 rm = (ri.u&(~0ul>>12))|1l<<52; int re = (ri.u>>52)-0x3ff;
  bm <<= bs;
  u64 m2 = bm*bm;
  int de = be-le;
  int ls = bs-de;
  if(__builtin_expect(ls>=0, 1)){
    lm <<= ls;
    m2 += lm*lm;
  } else {
    u128 lm2 = (u128)lm*lm;
    ls *= 2;
    m2 += lm2 >> -ls;
    m2 |= !!(lm2 << (128 + ls));
  }
  int k = bs+re, sr = 2*k+1;
  u64 tm = rm << k, rm2 = tm*tm;
  i64 D = m2 - rm2, pD, sD = D>>63, um = (rm^sD) - sD, rms = um<<sr, drm = sD + 1;
  u64 off = 1<<2*k, off2 = off<<1, dd = rms + off;
  rm -= drm;
  drm += sD;
  do {
    pD = D;
    rm += drm;
    D -= dd;
    dd += off2;
  } while(__builtin_expect((D^pD)>0, 0));
  pD = (sD&D)+(~sD&pD);
  if(__builtin_expect(pD != 0, 1)){
    if(__builtin_expect(op == om, 1)){
      u64 t = (rm<<2) + 1;
      unsigned s = 2*k - 2;
      i64 rms1 = t<<s;
      i64 sum = pD - rms1;
      if(__builtin_expect(sum, 1))
	rm += (sum>>63) + 1;
      else
	rm += rm&1;
    } else {
      rm += op>1.0;
    }
  }

  i64 e = be - 1 + re;
  xi.u = (e<<52) + rm;
  return xi.f;
}

double cr_hypot(double x, double y){
  b64u64_u xi = {.f = x}, yi = {.f = y};
  u64 emsk = 0x7ffl<<52, ex = xi.u&emsk, ey = yi.u&emsk;
  /* emsk corresponds to the upper bits of NaN and Inf (apart the sign bit) */
  x = __builtin_fabs(x), y = __builtin_fabs(y);
  if(__builtin_expect(ex==emsk||ey==emsk, 0)){
    /* Either x or y is NaN or Inf */
    u64 wx = xi.u<<1, wy = yi.u<<1, wm = emsk<<1;
    int ninf = (wx==wm) ^ (wy==wm);
    int nqnn = ((wx>>52)==0xfff) ^ ((wy>>52)==0xfff);
    /* ninf is 1 when only one of x and y is +/-Inf
       nqnn is 1 when only one of x and y is qNaN
       IEEE 754 says that hypot(+/-Inf,qNaN)=hypot(qNaN,+/-Inf)=+Inf. */
    if(ninf && nqnn) return __builtin_inf ();
    return x + y; /* inf, nan */
  }
  double u = __builtin_fmax(x,y), v = __builtin_fmin(x,y);
  __m128d xd, yd;
  asm("":"=x"(xd):"0"(u):);
  asm("":"=x"(yd):"0"(v):);
  __m128i de = (__m128i)xd - (__m128i)yd;
  if(__builtin_expect(de[0]>(27l<<52),0)) return __builtin_fma(0x1p-27, u, v);
  ey = ((__m128i)yd)[0];
  if(__builtin_expect(!(ey>>52),0)){
    ex = ((__m128i)xd)[0];
    if(__builtin_expect(!(ex>>52),0)){
      if(!ex) return 0;
      return as_hypot_denorm(ex,ey);
    }
    int nz = __builtin_clzll(ey);
    ey <<= nz-11;
    ey &= ~0ul>>12;
    ey -= (nz-12l)<<52;
    b64u64_u t = {.u = ey};
    yd[0] = t.f;
  }
  __m128i msk = {emsk, emsk}, off = {0x3ffl<<52, 0x3ffl<<52};
  off -= ((__m128i)xd & msk);
  xd = (__m128d)((__m128i)xd + off);
  yd = (__m128d)((__m128i)yd + off);
  x = xd[0];
  y = yd[0];
  double x2 = x*x, dx2 = __builtin_fma(x,x,-x2);
  double y2 = y*y, dy2 = __builtin_fma(y,y,-y2);
  double r2 = x2 + y2, ir2 = 0.5/r2, dr2 = ((x2 - r2) + y2) + (dx2 + dy2);
  double th = __builtin_sqrt(r2), rsqrt = th*ir2;
  double dz = dr2 - __builtin_fma(th,th,-r2), tl = rsqrt*dz;
  th = fasttwosum(th, tl, &tl);
  __m128d thd; asm("":"=x"(thd):"0"(th):);
  __m128d tld; asm("":"=x"(tld):"0"(__builtin_fabs(tl)):);
  ex = ((__m128i)thd)[0];
  ey = ((__m128i)tld)[0];
  ex &= 0x7ffl<<52;
  u64 aidr = ey + (0x3fel<<52) - ex;
  u64 mid = (aidr - 0x3c90000000000000 + 16)>>5;
  if(__builtin_expect( mid==0 || aidr<0x39b0000000000000ul || aidr>0x3c9fffffffffff80ul, 0)) 
    thd[0] = as_hypot_hard(x,y);
  thd = (__m128d)((__m128i)thd - off);
  if(__builtin_expect((u64)(((__m128i)thd)[0])>=(0x7fful<<52), 0)){
    errno = ERANGE;
    return 0x1.fffffffffffffp1023 + 0x1.fffffffffffffp1023;
  }
  return thd[0];
}
