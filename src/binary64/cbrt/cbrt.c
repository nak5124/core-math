/* Correctly-rounded cubic root of binary64 value.

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
#include <x86intrin.h>

typedef union {double f; uint64_t u;} b64u64_u;
typedef uint64_t u64;
typedef unsigned __int128 u128;

static double __attribute__((noinline)) as_cbrt_refine(double x, double r){
  unsigned flag = _mm_getcsr();
  b64u64_u cvt0 = {.f = x};
  u64 isgn = cvt0.u>>63;
  u64 hx = cvt0.u;
  unsigned e = (hx>>52)&0x7ff;
  if(__builtin_expect(e==0, 0)){
    u64 ix = hx&((~0ul)>>1);
    int nz = __builtin_clzll(ix);  /* denormals */
    hx <<= nz - 11;
    e -= nz - 12;
  }
  e += 2046;
  unsigned it = e%3, rm = (flag>>13)&3;
  b64u64_u ir = {.f = r};
  u64 z = (ir.u&(~0ul>>11)) | (1l<<52);
  u128 z2 = (u128)z*z;
  __int128 z3 = z2*z;
  u64 t = hx<<(40+it);
  long z3h = z3>>64; u64 z3l = z3;
  z3h -= t;
  z3 = (u128)z3h<<64|z3l;
  if(__builtin_expect(z3==0, 0)) return ir.f;
  long zs = ~(z3h>>63), dr = (zs<<1)+1;
  long d0 = (zs^6l) - zs, d1 = 6*z + d0, z3u = z3h;
  u128 d2 = ((__int128)(3*z2 + 1)^zs) - zs;
  d2 += 3*z;
  z3 += d2;
  ir.u += dr;
  if(__builtin_expect(z3==0, 0)) return ir.f;
  z3u = z3>>64;
  if( __builtin_expect(!((z3u^z3h) < 0), 1)){
    z3 += d2 + d1;
    ir.u += dr;
    if(__builtin_expect(z3==0, 0)) return ir.f;
  }
  ir.u -= zs;
  z3 -= d2&zs;
  if(rm==0){
    z = (ir.u&(~0ul>>11)) | 1l<<52;
    z2 = (u128)z*z;
    z3 -= (3*((z2+z2)-z))>>2;
    long d = z3>>127;
    ir.u -= 1 + d;
  } else {
    if (rm==1) {
      ir.u -= 1 - isgn;
    } else if(rm==2){
      ir.u -= isgn;
    } else {
      ir.u -= 1;
    }
  }
  return ir.f;
}

double cr_cbrt(double x){
  static const double escale[3] = {1.0, 0x1.428a2f98d728bp+0, 0x1.965fea53d6e3dp+0};
  static const double su[3] = {1.0, 2, 4};
  static const double sd[3] = {1.0, 0.5, 0.25};
  static const double eps[3] = {0.899e-18, 1.133e-18, 1.427e-18};
  static const double c[] = {
    0x1.22fe0d2edda62p-1, 0x1.67f254bb67748p-1, -0x1.9403dfa7453c5p-2, 0x1.b787fa3ff961ep-3,
    -0x1.6174462425c15p-4, 0x1.7d0352230cd22p-6, -0x1.e86777682f2dcp-9, 0x1.18ae3c4e5c285p-12,
    -0x1.98961922f4f6dp-6};
  b64u64_u cvt0 = {.f = x};
  u64 isgn = cvt0.u>>63;
  u64 hx = cvt0.u;
  cvt0.u = (hx&(~0ul>>12))|(0x3ffl<<52); double z = cvt0.f;
  unsigned e = (hx>>52)&0x7ff;
  if(__builtin_expect((unsigned)(e-1)>=0x7feu, 0)){
    u64 ix = hx&((~0ul)>>1);
    if((e==0x7ff)||ix==0) return x + x; /* 0, inf, nan */
    int nz = __builtin_clzll(ix);  /* denormals */
    hx <<= nz - 11;
    e -= nz - 12;
    b64u64_u cvt1 = {.u = hx|(0x3ffl<<52)};
    z = cvt1.f;
  }
  double rz = 1/z, z2 = z*z, z4 = z2*z2;
  e += 2046;
  unsigned et = e/3, it = e%3;
  cvt0.u = (et|isgn<<11)<<52;
  double c0 = c[0] + z*c[1], c2 = c[2] + z*c[3], c4 = c[4] + z*c[5], c6 = c[6] + z*c[7];
  z *= su[it];
  c0 += z2*c2;
  c4 += z2*c6;
  double y = escale[it]*((c0 + z4*c4) + c[8]*rz);
  double y2 = y*y, y2l = __builtin_fma(y, y, -y2);
  double h = __builtin_fma(y2, y, -z) + y2l*y, dy = (-0x1.5555555555555p-2)*sd[it]*rz*y*h;
  double dyp = dy + eps[it], rp = y + dyp, rm = y + dy;
  y *= cvt0.f;
  dy *= cvt0.f;
  y += dy;
  if(__builtin_expect(rp != rm, 0)) return as_cbrt_refine(x,y);
  return y;
}
