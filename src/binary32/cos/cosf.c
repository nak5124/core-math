/* Correctly-rounded cosine of binary32 value.

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

Tested on x86_64-linux with and without FMA (-march=native).
*/

#include <stdint.h>
#include <errno.h>
#include <fenv.h>
#include <math.h>

typedef union {float f; unsigned u;} b32u32_u;
typedef union {double f; unsigned long u;} b64u64_u;
typedef unsigned __int128 u128;
typedef unsigned long u64;

static inline double rltl(float z, int *q){
  double x = z;
  double idl = -0x1.b1bbead603d8bp-32*x, idh = 0x1.45f306ep-1*x, id = __builtin_roundeven(idh);
  *q = (long)id;
  return (idh - id) + idl;
}

static double __attribute__((noinline)) rbig(unsigned u, int *q){
  static const u64 ipi[] = {0xfe5163abdebbc562, 0xdb6295993c439041, 0xfc2757d1f534ddc0, 0xa2f9836e4e441529};
  int e = (u>>23)&0xff, i;
  u64 m = (u&(~0u>>9))|1<<23;
  u128 p0 = (u128)m*ipi[0];
  u128 p1 = (u128)m*ipi[1]; p1 += p0>>64;
  u128 p2 = (u128)m*ipi[2]; p2 += p1>>64;
  u128 p3 = (u128)m*ipi[3]; p3 += p2>>64;
  u64 p3h = p3>>64, p3l = p3, p2l = p2, p1l = p1;
  long a;
  int k = e-127, s = k-23;
  if(s<0){
    i =        p3h>>-s;
    a = p3h<<s|p3l>>-s;
  } else if(s==0) {
    i = p3h;
    a = p3l;
  } else if(s<64) {
    i = p3h<<s|p3l>>-s;
    a = p3l<<s|p2l>>-s;
  } else if(s==64) {
    i = p3l;
    a = p2l;
  } else {
    i = p3l<<s|p2l>>-s;
    a = p2l<<s|p1l>>-s;
  }
  int sgn = u; sgn >>= 31;
  long sm = a>>63;
  i -= sm;
  double z = (a^sgn)*0x1p-64;
  i = (i^sgn) - sgn;
  *q = i;
  return z;
}

float cr_cosf(float x){
  static const double
    cs[] = {
    -0x1.a51a6625307bdp-2, 0x1.9f9cb402b97ffp-5, -0x1.86a8e46de2fc4p-9,
    0x1.ac67ffda7836cp-14, -0x1.337d0b88f8cb7p-19, 0x1.3417d9749e139p-25},
    cc[] = {
      -0x1.3bd3cc9be458bp+0, 0x1.03c1f081b078ep-2, -0x1.55d3c7dbfe042p-6,
      0x1.e1f4fb610f151p-11, -0x1.a6c9c224d18abp-16, 0x1.f3dbf0909677fp-22},
    q[] = {1.0, 0, -1.0, 0};

  b32u32_u t = {.f = x};
  int e = (t.u>>23)&0xff, i;
  double z;
  if (__builtin_expect(e<127+28, 1)){
    if (__builtin_expect(e<114, 0))
      return __builtin_fmaf(-x, x, 1.0f);
    z = rltl(x, &i);
  } else if (e<0xff) {
    z = rbig(t.u, &i);
  } else {
    if(t.u<<9) return x; // nan
    errno = EDOM;
    feraiseexcept(FE_INVALID);
    return __builtin_nanf("cinf"); // inf
  }
  double z2 = z*z, z4 = z2*z2;
  double ms = q[(i+1)&3], mc = q[i&3];
  z *= 0x1.921fb54442d18p+0*ms;
  if(__builtin_expect(z2<0x1p-25, 0)){
    if(!(i&1)){
      float a = z2;
      int j = (i>>1)&1;
      static const float tb[] = {1, -1};
      return __builtin_fmaf(((float)cc[0])*tb[j], a, tb[j]);
    }
  }
  double s0 = cs[0] + z2*cs[1];
  double s2 = cs[2] + z2*cs[3];
  double s4 = cs[4] + z2*cs[5];
  double rs = (z*z2)*(s0 + z4*(s2 + z4*s4)) + z;
  double c0 = cc[0] + z2*cc[1];
  double c2 = cc[2] + z2*cc[3];
  double c4 = cc[4] + z2*cc[5];
  double rc = (z2*mc)*(c0 + z4*(c2 + z4*c4)) + mc;
  double r = rs + rc;
  float or = r;
  b64u64_u tr = {.f = r}; u64 tail = (tr.u + 6)&(~0ul>>36);
  if(__builtin_expect(tail<=12, 0)){
    static const struct {union{float arg; unsigned uarg;}; float rh, rl;} st[] = {
      {{0x1.4555p+51f}, 0x1.115d7ep-1f, -0x1.fffffep-26f},
      {{0x1.3170fp+63f}, 0x1.fe2976p-1f, 0x1.fffffep-26f},
      {{0x1.119ae6p+115f}, 0x1.f3176ap-1f, 0x1.fffffep-26f},
      {{0x1.96344ep+117f}, 0x1.f8a9b4p-1f, 0x1.fffffep-26f},
      {{0x1.24f2eep+85f}, 0x1.af5c6ap-2f, -0x1.91995ap-54f},
      {{0x1.ddebdep+120f}, 0x1.114438p-1f, 0x1.0b776ep-53f},
      {{0x1.ba5a98p+95f}, -0x1.ac093cp-1f, -0x1.7b13cep-53f},
      {{0x1.a13cb2p+1f}, -0x1.fc6f6cp-1f, 0x1.28b55ap-53f},
      {{0x1.417a5cp+8f}, 0x1.052988p-1f, 0x1.2f43fp-53f},
      {{0x1.34486p+19f}, -0x1.edfe3p-1f, 0x1.4ff1aep-55f},
      {{0x1.48a858p+54f}, 0x1.f48148p-2f, 0x1.e09b5p-56f},
      {{0x1.f4a7f4p+106f}, -0x1.b503dap-1f, 0x1.0c60e4p-53f},
      {{0x1.0a1f74p+58f}, -0x1.404ecep-2f, 0x1.58808ap-54f},
      {{0x1.ea56e2p+73f}, -0x1.d74dfap-3f, 0x1.7fa9c2p-55f},
    };
    unsigned ax = t.u&(~0u>>1);
    for(int i=0;i<14;i++) {
      if(__builtin_expect(st[i].uarg == ax, 0))
	return st[i].rh + st[i].rl;
    }
  }
  return or;
}
