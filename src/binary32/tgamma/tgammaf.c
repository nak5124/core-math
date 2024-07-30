/* Correctly-rounded true gamma function for binary32 value.

Copyright (c) 2023-2024 Alexei Sibidanov.

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

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;

float cr_tgammaf(float x){
  static const struct {b32u32_u x; float f, df;} tb[] = {
    {{.u = 0x27de86a9u}, 0x1.268266p+47f, 0x1p22f},
    {{.u = 0x27e05475u}, 0x1.242422p+47f, 0x1p22f},
    {{.u = 0xb63befb3u}, -0x1.5cb6e4p+18f, 0x1p-7f},
    {{.u = 0x3c7bb570u}, 0x1.021d9p+6f, 0x1p-19f},
    {{.u = 0x41e886d1u}, 0x1.33136ap+98f, 0x1p73f},
    {{.u = 0xc067d177u}, 0x1.f6850cp-3f, 0x1p-28f},
    {{.f = -0x1.33b462p-4}, -0x1.befe66p+3, -0x1p-22f},
    {{.f = -0x1.a988b4p-1}, -0x1.a6b4ecp+2, +0x1p-23f},
    {{.f = 0x1.dceffcp+4}, 0x1.d3631cp+101, -0x1p-76f},
    {{.f = 0x1.0874c8p+0}, 0x1.f6c638p-1, 0x1p-26f},
  };

  b32u32_u t = {.f = x};
  uint32_t ax = t.u<<1;
  if(__builtin_expect(ax>=(0xffu<<24), 0)){
    if(ax==(0xffu<<24)){
      if(t.u>>31){
	errno = EDOM;
	return __builtin_nanf("12");
      }
      return x;
    }
    return x + x; /* case x=NaN, where x+x ensures the invalid exception is
                     set if x=sNaN */
  }
  double z = x;
  if(__builtin_expect(ax<0x6d000000u, 0)){
    volatile double d = (0x1.fa658c23b1578p-1 - 0x1.d0a118f324b63p-1*z)*z - 0x1.2788cfc6fb619p-1;
    double f = 1.0/z + d;
    float r = f;
    if(__builtin_fabsf(r)>0x1.fffffep+127f) errno = ERANGE;
    b64u64_u rt = {.f = f};
    if(((rt.u+2)&0xfffffff) < 4){
      for(unsigned i=0;i<sizeof(tb)/sizeof(tb[0]);i++)
	if(t.u==tb[i].x.u) return tb[i].f + tb[i].df;
    }
    return r;
  }
  float fx = __builtin_floorf(x);
  if(__builtin_expect(x >= 0x1.18522p+5f, 0)){
    /* The C standard says that if the function overflows,
       errno is set to ERANGE. */
    errno = ERANGE;
    return 0x1p127f * 0x1p127f;
  }
  /* compute k only after the overflow check, otherwise the case to integer
     might overflow */
  int k = fx;
  if(__builtin_expect(fx==x, 0)){
    if(x == 0.0f){
      errno = ERANGE;
      return 1.0f/x;
    }
    if(x < 0.0f) {
      errno = EDOM;
      return __builtin_nanf("12");
    }
    double t0 = 1, x0 = 1;
    for(int i=1; i<k; i++, x0 += 1.0) t0 *= x0;
    return t0;
  }
  if(__builtin_expect(x<-42.0f, 0)){
    // for x < -42, x non-integer, |gamma(x)| < 2^-151
    static const float sgn[2] = {0x1p-127, -0x1p-127};
    /* The C standard says that if the function underflows,
       errno is set to ERANGE. */
    errno = ERANGE;
    return 0x1p-127f * sgn[k&1];
  }
  static const double c[] =
    {0x1.c9a76be577123p+0, 0x1.8f2754ddcf90dp+0, 0x1.0d1191949419bp+0, 0x1.e1f42cf0ae4a1p-2,
     0x1.82b358a3ab638p-3, 0x1.e1f2b30cd907bp-5, 0x1.240f6d4071bd8p-6, 0x1.1522c9f3cd012p-8,
     0x1.1fd0051a0525bp-10, 0x1.9808a8b96c37ep-13, 0x1.b3f78e01152b5p-15, 0x1.49c85a7e1fd04p-18,
     0x1.471ca49184475p-19, -0x1.368f0b7ed9e36p-23, 0x1.882222f9049efp-23, -0x1.a69ed2042842cp-25};

  double m = z - 0x1.7p+1, i = __builtin_roundeven(m), step = __builtin_copysign(1.0,i);
  double d = m - i, d2 = d*d, d4 = d2*d2, d8 = d4*d4;
  double f = (c[0] + d*c[1]) + d2*(c[2] + d*c[3]) + d4*((c[4] + d*c[5]) + d2*(c[6] + d*c[7]))
    + d8*((c[8] + d*c[9]) + d2*(c[10] + d*c[11]) + d4*((c[12] + d*c[13]) + d2*(c[14] + d*c[15])));
  int jm = __builtin_fabs(i);
  double w = 1;
  if(jm){
    z -= 0.5 + step*0.5;
    w = z;
    for(int j=jm-1; j; j--) {z -= step; w *= z;}
  }
  if(i<=-0.5) w = 1/w;
  f *= w;
  b64u64_u rt = {.f = f};
  float r = f;
  if(__builtin_expect(r==0.0f, 0)) errno = ERANGE;
  // deal with exceptional cases
  if(__builtin_expect(((rt.u+2)&0xfffffff) < 8, 0)){
    for(unsigned j=0;j<sizeof(tb)/sizeof(tb[0]);j++) {
      if(t.u==tb[j].x.u) return tb[j].f + tb[j].df;
    }
  }
  return r;
}
