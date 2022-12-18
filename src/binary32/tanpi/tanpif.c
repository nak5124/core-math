/* Correctly-rounded tangent of binary32 value for angles in half-revolutions

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
#include <fenv.h>

typedef union {float f; uint32_t u;} b32u32_u;

float cr_tanpif(float x){
  b32u32_u ix = {.f = x};
  int32_t e = (ix.u>>23)&0xff;
  if(__builtin_expect(e == 0xff, 0)){
    if(!(ix.u << 9)){
      errno = EDOM;
      feraiseexcept (FE_INVALID);
      return __builtin_nanf("inf");
    }
    return x;
  }
  static const double tn[] = {0x1.921fb54442d18p-37, 0x1.4abbce6260fdbp-111, 0x1.466bc4dfb775dp-185, 0x1.4628bbec47de1p-259};
  static const double T[] =
    {0x0p+0, 0x1.927278a3b1162p-5, 0x1.936bb8c5b2da2p-4, 0x1.2fcac73a6064p-3, 0x1.975f5e0553158p-3, 0x1.007fa758626aep-2,
    0x1.36a08355c63dcp-2, 0x1.6e649f7d78649p-2, 0x1.a827999fcef32p-2, 0x1.e450e0d273e7ap-2, 0x1.11ab7190834ecp-1,
    0x1.32e1889047ffdp-1, 0x1.561b82ab7f99p-1, 0x1.7bb99ed2990cfp-1, 0x1.a43002ae4285p-1, 0x1.d00cbc7384d2ep-1,0x1p+0,
    0x1.1a73d55278c4bp+0, 0x1.37efd8d87607ep+0, 0x1.592d11142fa55p+0, 0x1.7f218e25a7461p+0, 0x1.ab1c35d8a74eap+0,
    0x1.def13b73c1406p+0, 0x1.0ea21d716fbf7p+1, 0x1.3504f333f9de6p+1, 0x1.65bc6cc825147p+1, 0x1.a5f59e90600ddp+1,
    0x1.ff01305ecd8dcp+1, 0x1.41bfee2424771p+2, 0x1.af73f4ca3310fp+2, 0x1.44e6c595afdccp+3, 0x1.45affed201b55p+4,
    0x1p989, -0x1.45affed201b55p+4, -0x1.44e6c595afdccp+3, -0x1.af73f4ca3310fp+2, -0x1.41bfee2424771p+2,
    -0x1.ff01305ecd8dcp+1, -0x1.a5f59e90600ddp+1, -0x1.65bc6cc825147p+1, -0x1.3504f333f9de6p+1, -0x1.0ea21d716fbf7p+1,
    -0x1.def13b73c1406p+0, -0x1.ab1c35d8a74eap+0, -0x1.7f218e25a7461p+0, -0x1.592d11142fa55p+0, -0x1.37efd8d87607ep+0,
    -0x1.1a73d55278c4bp+0, -0x1p+0, -0x1.d00cbc7384d2ep-1, -0x1.a43002ae4285p-1, -0x1.7bb99ed2990cfp-1,
    -0x1.561b82ab7f99p-1, -0x1.32e1889047ffdp-1, -0x1.11ab7190834ecp-1, -0x1.e450e0d273e7ap-2, -0x1.a827999fcef32p-2,
    -0x1.6e649f7d78649p-2, -0x1.36a08355c63dcp-2, -0x1.007fa758626aep-2, -0x1.975f5e0553158p-3, -0x1.2fcac73a6064p-3,
    -0x1.936bb8c5b2da2p-4, -0x1.927278a3b1162p-5};
  int32_t m = (ix.u&~0u>>9)|1<<23, sgn = ix.u; sgn >>= 31;
  m = (m^sgn) - sgn;
  int32_t s = 143 - e;
  if(__builtin_expect(s<0,0)){
    if(__builtin_expect(s<-7,0)) return __builtin_copysign(0.0f, x);
    int32_t iq = m<<(-s-1);
    double ts = T[iq&63];
    if(ts>1e3) return ((iq&64)?-1.0f:1.0f)/0.0f;
    if(ts==0.0f) return ((iq&64)?-ts:ts)*__builtin_copysignf(1.0f,x);
    return ts;
  } else if(__builtin_expect(s>29, 0)){
    double z = x, z2 = z*z, r = z * (0x1.921fb54442d17p+1 + z2 *  0x1.4abbcfa8c7acfp+3);
    return r;
  }
  int32_t si = 26 - s;
  if(__builtin_expect(si>=0, 1)){
    if(__builtin_expect((m<<si)==0, 0)){
      b32u32_u o = {.u = ((m << (si - 2)) & (1<<31)) | 0x3f800000};
      if(!((m >> (32 - si))&1)) {
	return o.f*__builtin_copysignf(0.0f, x);
      } else {
	return o.f/0.0f;
      }
    }
  }
  int32_t k = m<<(31-s);
  double z = k, z2 = z*z, z4 = z2*z2;
  unsigned iq = m>>s; iq = (iq + 1)>>1;
  double ts = T[iq&63], fs = (tn[0] + z2*tn[1]) + z4*(tn[2] + z2*tn[3]);
  double r = (ts + z*fs)/(1 - (ts*z)*fs);
  return r;
}

float tanpif(float x){
  return cr_tanpif(x);
}
