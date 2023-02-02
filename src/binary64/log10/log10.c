/* Correctly rounded base-10 logarithm of binary64 values.

Copyright (c) 2022-2023 INRIA and CERN.
Authors: Paul Zimmermann and Tom Hubrecht.

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
#include "dint.h"

#define TRACE 0x1.89d948a94fe17p+2

typedef union { double f; uint64_t u; } d64u64;

/* Add a + b, such that *hi + *lo approximates a + b.
   Assumes |a| >= |b|.  */
static void
fast_two_sum (double *hi, double *lo, double a, double b)
{
  double e;

  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
  /* Now hi + lo = a + b exactly for rounding to nearest.
     For directed rounding modes, this is not always true.
     Take for example a = 1, b = 2^-200, and rounding up,
     then hi = 1 + 2^-52, e = 2^-52 (it can be proven that
     e is always exact), and lo = -2^52 + 2^-105, thus
     hi + lo = 1 + 2^-105 <> a + b = 1 + 2^-200.
     A bound on the error is given
     in "Note on FastTwoSum with Directed Roundings"
     by Paul Zimmermann, https://hal.inria.fr/hal-03798376, 2022.
     Theorem 1 says that
     the difference between a+b and hi+lo is bounded by 2u^2|a+b|
     and also by 2u^2|hi|. Here u=2^-53, thus we get:
     |(a+b)-(hi+lo)| <= 2^-105 min(|a+b|,|hi|) */
}

/* For 362 <= i <= 724, r[i] = _INVERSE[i-362] is a 10-bit approximation of
   1/x[i], where i*2^-9 <= x[i] < (i+1)*2^-9.
   More precisely r[i] is a 10-bit value such that r[i]*y-1 is representable
   exactly on 53 bits for for any y, i*2^-9 <= y < (i+1)*2^-9.
   Moreover |r[i]*y-1| <= 0.00212097167968735. */
static const double _INVERSE[363]= {
    0x1.698p+0, 0x1.688p+0, 0x1.678p+0, 0x1.668p+0, 0x1.658p+0, 0x1.648p+0, 0x1.638p+0,
    0x1.63p+0, 0x1.62p+0, 0x1.61p+0, 0x1.6p+0, 0x1.5fp+0, 0x1.5ep+0, 0x1.5dp+0,
    0x1.5cp+0, 0x1.5bp+0, 0x1.5a8p+0, 0x1.598p+0, 0x1.588p+0, 0x1.578p+0, 0x1.568p+0,
    0x1.56p+0, 0x1.55p+0, 0x1.54p+0, 0x1.53p+0, 0x1.52p+0, 0x1.518p+0, 0x1.508p+0,
    0x1.4f8p+0, 0x1.4fp+0, 0x1.4ep+0, 0x1.4dp+0, 0x1.4cp+0, 0x1.4b8p+0, 0x1.4a8p+0,
    0x1.4ap+0, 0x1.49p+0, 0x1.48p+0, 0x1.478p+0, 0x1.468p+0, 0x1.458p+0, 0x1.45p+0,
    0x1.44p+0, 0x1.43p+0, 0x1.428p+0, 0x1.418p+0, 0x1.41p+0, 0x1.4p+0, 0x1.3f8p+0,
    0x1.3e8p+0, 0x1.3ep+0, 0x1.3dp+0, 0x1.3cp+0, 0x1.3b8p+0, 0x1.3a8p+0, 0x1.3ap+0,
    0x1.39p+0, 0x1.388p+0, 0x1.378p+0, 0x1.37p+0, 0x1.36p+0, 0x1.358p+0, 0x1.35p+0,
    0x1.34p+0, 0x1.338p+0, 0x1.328p+0, 0x1.32p+0, 0x1.31p+0, 0x1.308p+0, 0x1.3p+0,
    0x1.2fp+0, 0x1.2e8p+0, 0x1.2d8p+0, 0x1.2dp+0, 0x1.2c8p+0, 0x1.2b8p+0, 0x1.2bp+0,
    0x1.2ap+0, 0x1.298p+0, 0x1.29p+0, 0x1.28p+0, 0x1.278p+0, 0x1.27p+0, 0x1.26p+0,
    0x1.258p+0, 0x1.25p+0, 0x1.24p+0, 0x1.238p+0, 0x1.23p+0, 0x1.228p+0, 0x1.218p+0,
    0x1.21p+0, 0x1.208p+0, 0x1.2p+0, 0x1.1fp+0, 0x1.1e8p+0, 0x1.1ep+0, 0x1.1dp+0,
    0x1.1c8p+0, 0x1.1cp+0, 0x1.1b8p+0, 0x1.1bp+0, 0x1.1ap+0, 0x1.198p+0, 0x1.19p+0,
    0x1.188p+0, 0x1.18p+0, 0x1.17p+0, 0x1.168p+0, 0x1.16p+0, 0x1.158p+0, 0x1.15p+0,
    0x1.14p+0, 0x1.138p+0, 0x1.13p+0, 0x1.128p+0, 0x1.12p+0, 0x1.118p+0, 0x1.11p+0,
    0x1.1p+0, 0x1.0f8p+0, 0x1.0fp+0, 0x1.0e8p+0, 0x1.0ep+0, 0x1.0d8p+0, 0x1.0dp+0,
    0x1.0c8p+0, 0x1.0cp+0, 0x1.0bp+0, 0x1.0a8p+0, 0x1.0ap+0, 0x1.098p+0, 0x1.09p+0,
    0x1.088p+0, 0x1.08p+0, 0x1.078p+0, 0x1.07p+0, 0x1.068p+0, 0x1.06p+0, 0x1.058p+0,
    0x1.05p+0, 0x1.048p+0, 0x1.04p+0, 0x1.038p+0, 0x1.03p+0, 0x1.028p+0, 0x1.02p+0,
    0x1.018p+0, 0x1.01p+0, 0x1.008p+0, 0x1.ff8p-1, 0x1.fe8p-1, 0x1.fd8p-1, 0x1.fc8p-1,
    0x1.fb8p-1, 0x1.fa8p-1, 0x1.f98p-1, 0x1.f88p-1, 0x1.f78p-1, 0x1.f68p-1, 0x1.f58p-1,
    0x1.f5p-1, 0x1.f4p-1, 0x1.f3p-1, 0x1.f2p-1, 0x1.f1p-1, 0x1.fp-1, 0x1.efp-1,
    0x1.eep-1, 0x1.edp-1, 0x1.ec8p-1, 0x1.eb8p-1, 0x1.ea8p-1, 0x1.e98p-1, 0x1.e88p-1,
    0x1.e78p-1, 0x1.e7p-1, 0x1.e6p-1, 0x1.e5p-1, 0x1.e4p-1, 0x1.e3p-1, 0x1.e28p-1,
    0x1.e18p-1, 0x1.e08p-1, 0x1.df8p-1, 0x1.dfp-1, 0x1.dep-1, 0x1.ddp-1, 0x1.dcp-1,
    0x1.db8p-1, 0x1.da8p-1, 0x1.d98p-1, 0x1.d9p-1, 0x1.d8p-1, 0x1.d7p-1, 0x1.d6p-1,
    0x1.d58p-1, 0x1.d48p-1, 0x1.d38p-1, 0x1.d3p-1, 0x1.d2p-1, 0x1.d1p-1, 0x1.d08p-1,
    0x1.cf8p-1, 0x1.ce8p-1, 0x1.cep-1, 0x1.cdp-1, 0x1.cc8p-1, 0x1.cb8p-1, 0x1.ca8p-1,
    0x1.cap-1, 0x1.c9p-1, 0x1.c88p-1, 0x1.c78p-1, 0x1.c68p-1, 0x1.c6p-1, 0x1.c5p-1,
    0x1.c48p-1, 0x1.c38p-1, 0x1.c3p-1, 0x1.c2p-1, 0x1.c18p-1, 0x1.c08p-1, 0x1.bf8p-1,
    0x1.bfp-1, 0x1.bep-1, 0x1.bd8p-1, 0x1.bc8p-1, 0x1.bcp-1, 0x1.bbp-1, 0x1.ba8p-1,
    0x1.b98p-1, 0x1.b9p-1, 0x1.b8p-1, 0x1.b78p-1, 0x1.b68p-1, 0x1.b6p-1, 0x1.b58p-1,
    0x1.b48p-1, 0x1.b4p-1, 0x1.b3p-1, 0x1.b28p-1, 0x1.b18p-1, 0x1.b1p-1, 0x1.bp-1,
    0x1.af8p-1, 0x1.afp-1, 0x1.aep-1, 0x1.ad8p-1, 0x1.ac8p-1, 0x1.acp-1, 0x1.ab8p-1,
    0x1.aa8p-1, 0x1.aap-1, 0x1.a9p-1, 0x1.a88p-1, 0x1.a8p-1, 0x1.a7p-1, 0x1.a68p-1,
    0x1.a6p-1, 0x1.a5p-1, 0x1.a48p-1, 0x1.a4p-1, 0x1.a3p-1, 0x1.a28p-1, 0x1.a2p-1,
    0x1.a1p-1, 0x1.a08p-1, 0x1.ap-1, 0x1.9fp-1, 0x1.9e8p-1, 0x1.9ep-1, 0x1.9dp-1,
    0x1.9c8p-1, 0x1.9cp-1, 0x1.9bp-1, 0x1.9a8p-1, 0x1.9ap-1, 0x1.998p-1, 0x1.988p-1,
    0x1.98p-1, 0x1.978p-1, 0x1.968p-1, 0x1.96p-1, 0x1.958p-1, 0x1.95p-1, 0x1.94p-1,
    0x1.938p-1, 0x1.93p-1, 0x1.928p-1, 0x1.92p-1, 0x1.91p-1, 0x1.908p-1, 0x1.9p-1,
    0x1.8f8p-1, 0x1.8e8p-1, 0x1.8ep-1, 0x1.8d8p-1, 0x1.8dp-1, 0x1.8c8p-1, 0x1.8b8p-1,
    0x1.8bp-1, 0x1.8a8p-1, 0x1.8ap-1, 0x1.898p-1, 0x1.888p-1, 0x1.88p-1, 0x1.878p-1,
    0x1.87p-1, 0x1.868p-1, 0x1.86p-1, 0x1.85p-1, 0x1.848p-1, 0x1.84p-1, 0x1.838p-1,
    0x1.83p-1, 0x1.828p-1, 0x1.82p-1, 0x1.81p-1, 0x1.808p-1, 0x1.8p-1, 0x1.7f8p-1,
    0x1.7fp-1, 0x1.7e8p-1, 0x1.7ep-1, 0x1.7d8p-1, 0x1.7c8p-1, 0x1.7cp-1, 0x1.7b8p-1,
    0x1.7bp-1, 0x1.7a8p-1, 0x1.7ap-1, 0x1.798p-1, 0x1.79p-1, 0x1.788p-1, 0x1.78p-1,
    0x1.778p-1, 0x1.77p-1, 0x1.76p-1, 0x1.758p-1, 0x1.75p-1, 0x1.748p-1, 0x1.74p-1,
    0x1.738p-1, 0x1.73p-1, 0x1.728p-1, 0x1.72p-1, 0x1.718p-1, 0x1.71p-1, 0x1.708p-1,
    0x1.7p-1, 0x1.6f8p-1, 0x1.6fp-1, 0x1.6e8p-1, 0x1.6ep-1, 0x1.6d8p-1, 0x1.6dp-1,
    0x1.6c8p-1, 0x1.6cp-1, 0x1.6b8p-1, 0x1.6bp-1, 0x1.6a8p-1, 0x1.6ap-1,
};

/* For 362 <= i <= 724, (h,l) = _LOG10_INV[i-362] is a double-double
   approximation of -log10(r) with r=INVERSE[i-362]), with h an integer
   multiple of 2^-42, and |l| < 2^-43. The maximal difference between
   -log10(r) and h+l is bounded by 1/2 ulp(l) < 2^-97. */
static const double _LOG10_INV[363][2] = {
    {-0x1.32ee2b998ep-3, -0x1.adc8525d9f1b1p-44}, /* i=362*/                   
    {-0x1.30776f4d08p-3, 0x1.012fd64d71996p-44}, /* i=363*/                    
    {-0x1.2dfef27a74p-3, 0x1.7a050c54c60c1p-45}, /* i=364*/                    
    {-0x1.2b84b2a226p-3, 0x1.b74a059912c45p-45}, /* i=365*/                    
    {-0x1.2908ad3f14p-3, 0x1.3f99affe412a7p-44}, /* i=366*/                    
    {-0x1.268adfc6c8p-3, 0x1.1848822a01d4bp-44}, /* i=367*/                    
    {-0x1.240b47a95p-3, -0x1.ef1762dbe5406p-44}, /* i=368*/                    
    {-0x1.22cacece26p-3, -0x1.d590237ba79aep-44}, /* i=369*/                   
    {-0x1.204881dee8p-3, -0x1.ddd54b04da9d8p-45}, /* i=370*/                   
    {-0x1.1dc463ca42p-3, 0x1.03dc97d1e9c72p-46}, /* i=371*/                    
    {-0x1.1b3e71ec94p-3, -0x1.ef57776664942p-44}, /* i=372*/                   
    {-0x1.18b6a99c8p-3, 0x1.30e7599e7dcddp-44}, /* i=373*/                     
    {-0x1.162d082acap-3, 0x1.7838c72e86b79p-46}, /* i=374*/                    
    {-0x1.13a18ae256p-3, -0x1.7313a1aeda20cp-44}, /* i=375*/                   
    {-0x1.11142f0812p-3, 0x1.95237189e3611p-44}, /* i=376*/                    
    {-0x1.0e84f1dadcp-3, 0x1.5b3545950d96ep-44}, /* i=377*/                    
    {-0x1.0d3c9de722p-3, -0x1.e07e9228a44dep-49}, /* i=378*/                   
    {-0x1.0aaa89860cp-3, 0x1.0c332653781d9p-45}, /* i=379*/                    
    {-0x1.08168cd45ep-3, 0x1.eaeb5a9fb365ep-44}, /* i=380*/                    
    {-0x1.0580a4fb4ap-3, -0x1.ef8beb3d7a4dap-46}, /* i=381*/                   
    {-0x1.02e8cf1dacp-3, -0x1.2e0bd0e5ec1cdp-45}, /* i=382*/                   
    {-0x1.019c2a064cp-3, 0x1.6f31d0b132f8fp-44}, /* i=383*/                    
    {-0x1.fe02d36a38p-4, -0x1.56d63de758ea6p-44}, /* i=384*/                   
    {-0x1.f8c9683468p-4, -0x1.9084e03494e7dp-48}, /* i=385*/                   
    {-0x1.f38c0c8324p-4, -0x1.d85ad659b2175p-44}, /* i=386*/                   
    {-0x1.ee4aba611p-4, 0x1.bf71dec67fbefp-45}, /* i=387*/                     
    {-0x1.eba893055p-4, -0x1.00aa20eb3beafp-44}, /* i=388*/                    
    {-0x1.e66143f02cp-4, -0x1.8b9b32f90d3dap-45}, /* i=389*/                   
    {-0x1.e115ef491p-4, 0x1.1f90b2dcf62b3p-44}, /* i=390*/                     
    {-0x1.de6ec0f394p-4, 0x1.4fb46fcd56eddp-44}, /* i=391*/                    
    {-0x1.d91d5866acp-4, 0x1.66473a1327abcp-44}, /* i=392*/                    
    {-0x1.d3c7dacf58p-4, 0x1.fd4d14dc4b2d3p-46}, /* i=393*/                    
    {-0x1.ce6e41e464p-4, 0x1.2d85bc180e427p-47}, /* i=394*/                    
    {-0x1.cbbfe934b4p-4, 0x1.17f6d23315832p-44}, /* i=395*/                    
    {-0x1.c6601b6324p-4, 0x1.934b7876f38d2p-44}, /* i=396*/                    
    {-0x1.c3aea4a5c8p-4, 0x1.10162e464084cp-44}, /* i=397*/                    
    {-0x1.be4893762cp-4, -0x1.7ed341fd1b85ep-45}, /* i=398*/                   
    {-0x1.b8de4d3ab4p-4, 0x1.340511b402e0ap-47}, /* i=399*/                    
    {-0x1.b627942aecp-4, -0x1.dba29d8601f93p-45}, /* i=400*/                   
    {-0x1.b0b6f203acp-4, 0x1.f50fc58c04622p-46}, /* i=401*/                    
    {-0x1.ab420a41fp-4, -0x1.07629be399955p-44}, /* i=402*/                    
    {-0x1.a885fa2d6p-4, -0x1.51d88b2679d95p-44}, /* i=403*/                    
    {-0x1.a30a9d60ap-4, 0x1.0163d7e67d282p-44}, /* i=404*/                     
    {-0x1.9d8aea084cp-4, 0x1.56414945caae4p-44}, /* i=405*/                    
    {-0x1.9ac96dc174p-4, -0x1.775ce56bfce22p-44}, /* i=406*/                   
    {-0x1.95432ba8d4p-4, 0x1.1d142bd45c754p-44}, /* i=407*/                    
    {-0x1.927e64181p-4, 0x1.0e034db19292ep-45}, /* i=408*/                     
    {-0x1.8cf1838864p-4, -0x1.019365163f2fcp-45}, /* i=409*/                   
    {-0x1.8a2968c438p-4, -0x1.a82909fa92ddep-45}, /* i=410*/                   
    {-0x1.8495d9ce0cp-4, -0x1.0bebcc8c3bbb1p-48}, /* i=411*/                   
    {-0x1.81ca63d05cp-4, 0x1.bb67d8e7b2c03p-44}, /* i=412*/                    
    {-0x1.7c30164a6p-4, -0x1.06c11064a1f6ep-45}, /* i=413*/                    
    {-0x1.769140a254p-4, 0x1.9036b130dc2b9p-44}, /* i=414*/                    
    {-0x1.73c02075ep-4, 0x1.8354653757401p-45}, /* i=415*/                     
    {-0x1.6e1a70cb0cp-4, 0x1.998af72eaa028p-46}, /* i=416*/                    
    {-0x1.6b45df6f4p-4, 0x1.d36a6f1f2ab39p-44}, /* i=417*/                     
    {-0x1.659944f8bcp-4, 0x1.fd3676723abe2p-44}, /* i=418*/                    
    {-0x1.62c139f9b4p-4, 0x1.f257e93f31b79p-46}, /* i=419*/                    
    {-0x1.62c139f9b4p-4, 0x1.f257e93f31b79p-46}, /* i=419*/                    
    {-0x1.5d0da3b09cp-4, 0x1.b81de7b6c0dfdp-44}, /* i=420*/                    
    {-0x1.5a32167b34p-4, 0x1.0fdce69e470e8p-44}, /* i=421*/                    
    {-0x1.5477731974p-4, 0x1.7b786f3ec11dcp-48}, /* i=422*/                    
    {-0x1.51985afaap-4, 0x1.017470a82ba2p-47}, /* i=423*/                      
    {-0x1.48f3ed1df4p-4, -0x1.1f6bc109076dp-45}, /* i=425*/                    
    {-0x1.46100e075p-4, -0x1.b4237039b5162p-45}, /* i=426*/                  
    {-0x1.4044b2285cp-4, -0x1.9249c88a644fbp-44}, /* i=427*/                   
    {-0x1.3d5d335c54p-4, 0x1.d0e6a0f629baep-45}, /* i=428*/                    
    {-0x1.378a8ef848p-4, -0x1.71e18847c4d89p-44}, /* i=429*/                   
    {-0x1.349f6754ecp-4, -0x1.0b3af7943bc15p-44}, /* i=430*/                   
    {-0x1.31b3055c48p-4, 0x1.dcffc97be8c9cp-45}, /* i=431*/                    
    {-0x1.2bd68e462p-4, -0x1.37167d7ea5dffp-44}, /* i=432*/                    
    {-0x1.28e67712d8p-4, -0x1.dfbbdcdb9df2cp-44}, /* i=433*/                   
    {-0x1.23028c1b7cp-4, -0x1.aecff1b44ae64p-44}, /* i=434*/                   
    {-0x1.200eb639a4p-4, 0x1.d1ae17927dc54p-45}, /* i=435*/                    
    {-0x1.1d199ea83cp-4, 0x1.07165a9c258d4p-50}, /* i=436*/                    
    {-0x1.172ba62c4cp-4, -0x1.6ddee0592d983p-44}, /* i=437*/                   
    {-0x1.1432c31918p-4, 0x1.7c1114d5dc0b3p-47}, /* i=438*/                    
    {-0x1.0e3d29d81p-4, -0x1.65e62559618f2p-44}, /* i=439*/                    
    {-0x1.0b40717934p-4, 0x1.469948de978f2p-44}, /* i=440*/                    
    {-0x1.08426fcdbp-4, -0x1.ee6e333b614f5p-44}, /* i=441*/                    
    {-0x1.02428c1f08p-4, -0x1.5ea6bc2bc8c2cp-52}, /* i=442*/                   
    {-0x1.fe814fbec8p-5, 0x1.2c44f9cb9781cp-44}, /* i=443*/                    
    {-0x1.f87aebb44p-5, 0x1.8fd11436f9361p-44}, /* i=444*/                     
    {-0x1.ec6647eb58p-5, -0x1.0108fa031185ap-46}, /* i=445*/                   
    {-0x1.e658039c88p-5, -0x1.a3f5067fd6fabp-47}, /* i=446*/                   
    {-0x1.e0471aa188p-5, 0x1.70b753590c5d3p-45}, /* i=447*/                    
    {-0x1.d41d5164f8p-5, -0x1.659d00c475908p-44}, /* i=448*/                   
    {-0x1.ce046c7ad8p-5, -0x1.346a323a2bb75p-44}, /* i=449*/                   
    {-0x1.c7e8d9935p-5, -0x1.3f17c624bd312p-46}, /* i=450*/                    
    {-0x1.c1ca96526p-5, 0x1.85524e64e2a5fp-44}, /* i=451*/                     
    {-0x1.b585f54498p-5, 0x1.72e382fd9e528p-44}, /* i=452*/                    
    {-0x1.af5f92b01p-5, 0x1.9f05921f59258p-45}, /* i=453*/                     
    {-0x1.a9367632bp-5, 0x1.3b49bc8394523p-44}, /* i=454*/                     
    {-0x1.a30a9d60ap-5, 0x1.0163d7e67d282p-45}, /* i=455*/                     
    {-0x1.96aaacffp-5, 0x1.862239186139bp-44}, /* i=456*/                      
    {-0x1.907690878p-5, -0x1.b71d9fe1dd5c2p-44}, /* i=457*/                    
    {-0x1.8a3fadeb88p-5, 0x1.c063628960c25p-44}, /* i=458*/                    
    {-0x1.7dc98c51c8p-5, -0x1.2127595668247p-48}, /* i=459*/                   
    {-0x1.778a48519p-5, -0x1.bcd51444ab4fp-47}, /* i=460*/                     
    {-0x1.71483427dp-5, -0x1.54c670f08803p-44}, /* i=461*/                     
    {-0x1.6b034d4adp-5, -0x1.9efa9137a1fa4p-44}, /* i=462*/                    
    {-0x1.64bb912d68p-5, 0x1.1fc708031099dp-44}, /* i=463*/                    
    {-0x1.58238eeb38p-5, 0x1.612c205756103p-44}, /* i=464*/                    
    {-0x1.51d3439abp-5, 0x1.95ae2836c3efap-44}, /* i=465*/                     
    {-0x1.4b8018b22p-5, 0x1.2b0af21a91024p-45}, /* i=466*/                     
    {-0x1.452a0b92dp-5, 0x1.f8a145593666p-44}, /* i=467*/                      
    {-0x1.3ed1199a6p-5, 0x1.bdafc8ad828b8p-45}, /* i=468*/                     
    {-0x1.32167c82cp-5, 0x1.193325724a3ddp-44}, /* i=469*/                     
    {-0x1.2bb4cc0cbp-5, 0x1.72ffd7ffe731bp-49}, /* i=470*/                     
    {-0x1.25502c0fcp-5, -0x1.8a5c00ed6bef7p-44}, /* i=471*/                    
    {-0x1.1ee899d748p-5, -0x1.01f1010f86affp-44}, /* i=472*/                   
    {-0x1.187e12aad8p-5, -0x1.dd9adc1c7f97fp-51}, /* i=473*/                   
    {-0x1.0ba01a817p-5, -0x1.5f1d45244f437p-60}, /* i=474*/                    
    {-0x1.052ca400a8p-5, 0x1.8328c28c5de2ap-44}, /* i=475*/                    
    {-0x1.fd6c5b085p-6, -0x1.c4b8600163d9fp-46}, /* i=476*/                    
    {-0x1.f0796880dp-6, -0x1.c197a7259ab27p-46}, /* i=477*/                    
    {-0x1.e3806acbdp-6, -0x1.63c35e7d67688p-48}, /* i=478*/                    
    {-0x1.d6815c427p-6, -0x1.774e831b960aap-46}, /* i=479*/                    
    {-0x1.c97c3735ep-6, -0x1.f028fb72dfe8ap-44}, /* i=480*/                    
    {-0x1.af5f92b01p-6, 0x1.9f05921f59258p-46}, /* i=481*/                     
    {-0x1.a24807b0ep-6, -0x1.ad70142ffbdb7p-44}, /* i=482*/                    
    {-0x1.952a4f22cp-6, -0x1.6ba3837b618d8p-44}, /* i=483*/                    
    {-0x1.8806632e4p-6, -0x1.473fae5a8e918p-44}, /* i=484*/                    
    {-0x1.7adc3df3bp-6, -0x1.ff81b980714c6p-46}, /* i=485*/                    
    {-0x1.6dabd98bp-6, 0x1.9fd76cbaba4e4p-45}, /* i=486*/                      
    {-0x1.60753003bp-6, 0x1.ac42915db7ec8p-44}, /* i=487*/                     
    {-0x1.53383b64cp-6, 0x1.d82564e0ad4cep-47}, /* i=488*/                     
    {-0x1.45f4f5accp-6, 0x1.d07e22587685bp-44}, /* i=489*/                     
    {-0x1.2b5b5ec02p-6, -0x1.3a62b79ada68bp-47}, /* i=490*/                    
    {-0x1.1e05015d3p-6, -0x1.4e241ecdd26bfp-47}, /* i=491*/                  
    {-0x1.10a83a844p-6, -0x1.b1de84602abd8p-44}, /* i=492*/                    
    {-0x1.034504082p-6, -0x1.6eb34d4aa89ep-44}, /* i=493*/                     
    {-0x1.ebb6af654p-7, 0x1.d11c9508ca27ap-47}, /* i=494*/                     
    {-0x1.d0d65e89p-7, -0x1.0169545b91b96p-45}, /* i=495*/                     
    {-0x1.b5e908eb2p-7, 0x1.90dfe0d1601cap-44}, /* i=496*/                     
    {-0x1.9aeea1e8ap-7, 0x1.00429e06e1f32p-44}, /* i=497*/                     
    {-0x1.7fe71ccc4p-7, -0x1.cd60c6a5111d3p-44}, /* i=498*/                    
    {-0x1.64d26cce6p-7, -0x1.0dd3afe4cd0b2p-47}, /* i=499*/                    
    {-0x1.49b085144p-7, -0x1.b41e70df8592fp-46}, /* i=500*/                    
    {-0x1.2e8158b08p-7, -0x1.b2ae8a612cfe6p-44}, /* i=501*/                    
    {-0x1.1344daa2ep-7, 0x1.155ac9c1a811dp-44}, /* i=502*/                     
    {-0x1.eff5fbaf4p-8, 0x1.a87eeb5600788p-44}, /* i=503*/                     
    {-0x1.b9476a4fcp-8, -0x1.a21db136b482ep-45}, /* i=504*/                    
    {-0x1.827de6b3p-8, -0x1.03502db0555dp-44}, /* i=505*/                      
    {-0x1.4b99563d4p-8, 0x1.5e42f025b514ep-44}, /* i=506*/                     
    {-0x1.14999e2acp-8, -0x1.8ea5eaca88675p-44}, /* i=507*/                    
    {-0x1.bafd4722p-9, -0x1.ed2665c1ba949p-45}, /* i=508*/                     
    {-0x1.4c9096b98p-9, 0x1.a081515cbbf14p-44}, /* i=509*/                     
    {-0x1.bbd9e948p-10, -0x1.5784564411e7p-45}, /* i=510*/                     
    {-0x1.bc48a868p-11, 0x1.ded251d1ef535p-45}, /* i=511*/                     
    {0x1.bcef519p-12, -0x1.d69ee5af9439bp-44}, /* i=512*/                      
    {0x1.4e071755p-10, -0x1.b3a32c67b1bd2p-45}, /* i=513*/                     
    {0x1.16a117e1p-9, -0x1.5a680cacbe55ep-44}, /* i=514*/                      
    {0x1.8676c714p-9, 0x1.a7dc81de12997p-44}, /* i=515*/                       
    {0x1.f684d1d88p-9, 0x1.7380c2252ac38p-44}, /* i=516*/                      
    {0x1.3365b88cp-8, 0x1.e68e94df71746p-45}, /* i=517*/                       
    {0x1.6ba56f098p-8, -0x1.de43b4dad193bp-44}, /* i=518*/                     
    {0x1.a401a93p-8, -0x1.f6073154765d7p-46}, /* i=519*/                       
    {0x1.dc7a83f74p-8, 0x1.a96d29ea2d41bp-44}, /* i=520*/                      
    {0x1.0a880e41ap-7, 0x1.9fd7900a373dp-45}, /* i=521*/                       
    {0x1.26e148122p-7, -0x1.87a3d4b610e8bp-44}, /* i=522*/                     
    {0x1.351352a8ep-7, 0x1.ccfd2495d8b8ep-45}, /* i=523*/                      
    {0x1.51824c758p-7, 0x1.fabf59b5d80b8p-45}, /* i=524*/                      
    {0x1.6dffd8d3cp-7, -0x1.023f21feb5c47p-45}, /* i=525*/                     
    {0x1.8a8c06bb2p-7, -0x1.685fc114e61bfp-46}, /* i=526*/                     
    {0x1.a726e53a6p-7, 0x1.5b64be2b1b54p-49}, /* i=527*/                       
    {0x1.c3d083778p-7, 0x1.310272fe17537p-45}, /* i=528*/                      
    {0x1.e088f0bp-7, 0x1.209b0cfc0a6aep-45}, /* i=529*/                        
    {0x1.fd503c39p-7, 0x1.3c757d5b7376ap-45}, /* i=530*/                       
    {0x1.0d133abfcp-6, 0x1.f8d484ac7f8e4p-45}, /* i=531*/                      
    {0x1.144b98114p-6, -0x1.53fd31ec07e5fp-46}, /* i=532*/                     
    {0x1.22c1f5933p-6, -0x1.8392da13bf183p-44}, /* i=533*/                     
    {0x1.313fdd70fp-6, -0x1.a187fc7242c16p-46}, /* i=534*/                     
    {0x1.3fc5578b9p-6, 0x1.d1c33bd58c76ep-45}, /* i=535*/                      
    {0x1.4e526bd08p-6, 0x1.4d0cadec0287cp-46}, /* i=536*/                      
    {0x1.5ce72239ap-6, -0x1.d1a132eeb5289p-44}, /* i=537*/                     
    {0x1.64345cbd4p-6, -0x1.6dbd82e7594f2p-44}, /* i=538*/                     
    {0x1.72d4956cap-6, 0x1.0336603c83b1bp-45}, /* i=539*/                      
    {0x1.817c84683p-6, -0x1.d10a37cfc8c78p-44}, /* i=540*/                     
    {0x1.902c31d63p-6, -0x1.5ef22a507e5aep-44}, /* i=541*/                     
    {0x1.9ee3a5e9fp-6, 0x1.5fa037c49bb95p-44}, /* i=542*/                      
    {0x1.a6424d05ap-6, -0x1.cc18691f65161p-49}, /* i=543*/                     
    {0x1.b5057a8eep-6, 0x1.ce3fc4237e8a5p-50}, /* i=544*/                      
    {0x1.c3d083778p-6, 0x1.310272fe17537p-44}, /* i=545*/                      
    {0x1.d2a37021p-6, 0x1.496560fa0a672p-44}, /* i=546*/                       
    {0x1.da0fde804p-6, -0x1.c85d4a65cdeb6p-44}, /* i=547*/                     
    {0x1.e8eeb09f3p-6, -0x1.26ac877784097p-47}, /* i=548*/                     
    {0x1.c3d083778p-6, 0x1.310272fe17537p-44}, /* i=545*/                      
    {0x1.d2a37021p-6, 0x1.496560fa0a672p-44}, /* i=546*/                       
    {0x1.071f58e2dp-5, -0x1.03292ffa8c4abp-44}, /* i=551*/                     
    {0x1.0e9cc861d8p-5, -0x1.b5e226d4b03a1p-44}, /* i=552*/                    
    {0x1.161e4374cp-5, 0x1.bf9c87e354a83p-44}, /* i=553*/                      
    {0x1.19e086b3b8p-5, 0x1.99ac1fd443e4p-48}, /* i=554*/                      
    {0x1.21681b5c9p-5, -0x1.ef65393de7321p-44}, /* i=555*/                  
    {0x1.28f3c6991p-5, 0x1.53a1c756ef644p-44}, /* i=556*/                      
    {0x1.30838cdc3p-5, -0x1.00c12f7a1b586p-47}, /* i=557*/                     
    {0x1.344cfb8618p-5, 0x1.7d6dcea9ff6eep-44}, /* i=558*/                     
    {0x1.3be2f2ba78p-5, -0x1.7f0abe264207cp-44}, /* i=559*/                    
    {0x1.437d103498p-5, 0x1.dd73fc4a73ab1p-46}, /* i=560*/                     
    {0x1.474baeb78p-5, -0x1.6fc1197beb0bep-45}, /* i=561*/                     
    {0x1.4eec0e2458p-5, 0x1.e5ff3439d368dp-46}, /* i=562*/                     
    {0x1.56909f44a8p-5, -0x1.a04483513cf5bp-46}, /* i=563*/                    
    {0x1.5a647be9a8p-5, -0x1.584d21687f44dp-44}, /* i=564*/                    
    {0x1.620f604498p-5, -0x1.5e031024e5e2dp-44}, /* i=565*/                    
    {0x1.69be81f018p-5, 0x1.7ccab8a9dfd4ep-44}, /* i=566*/                     
    {0x1.6d97ab3ba8p-5, -0x1.0f7e884cbf05cp-44}, /* i=567*/                    
    {0x1.754d31b1bp-5, 0x1.79c3d52199ef2p-45}, /* i=568*/                      
    {0x1.7929900bd8p-5, -0x1.bfb442450bd02p-44}, /* i=569*/                    
    {0x1.80e585f92p-5, 0x1.8fc42622cabb9p-45}, /* i=570*/                      
    {0x1.88a5cc3158p-5, 0x1.e533f553fef1cp-45}, /* i=571*/                     
    {0x1.8c878eeb08p-5, -0x1.7c5d8e8ad876cp-44}, /* i=572*/                    
    {0x1.944e56a0dp-5, 0x1.a27d124156d9ap-44}, /* i=573*/                      
    {0x1.98335cd4ap-5, 0x1.6c2c0931d0032p-45}, /* i=574*/                      
    {0x1.a000b0fd1p-5, -0x1.4a4db81d2bcdfp-48}, /* i=575*/                     
    {0x1.a7d268eb5p-5, -0x1.b6df2546e12f3p-44}, /* i=576*/                     
    {0x1.abbcebd85p-5, -0x1.b0197d2cb982ep-48}, /* i=577*/                     
    {0x1.b39542ba2p-5, 0x1.eb996591c96ap-44}, /* i=578*/                       
    {0x1.b78317eefp-5, 0x1.a2974b12d552bp-45}, /* i=579*/                      
    {0x1.bf62190448p-5, 0x1.a43990040e4d6p-46}, /* i=580*/                     
    {0x1.c3534628p-5, 0x1.6dcbde98cd2abp-45}, /* i=581*/                       
    {0x1.cb38fccd88p-5, 0x1.fedb4b594a31bp-44}, /* i=582*/                     
    {0x1.cf2d8795c8p-5, 0x1.2d2a174d75277p-44}, /* i=583*/                     
    {0x1.d719ff456p-5, -0x1.411a14b5ff378p-46}, /* i=584*/                     
    {0x1.df0afe1508p-5, -0x1.463f99aa83d1ap-44}, /* i=585*/                    
    {0x1.e30531c77p-5, -0x1.e5aecd219cde7p-44}, /* i=586*/                     
    {0x1.eafd050358p-5, 0x1.e9d92d38dc40cp-44}, /* i=587*/                     
    {0x1.eefaa5dc28p-5, 0x1.91c924fc63616p-44}, /* i=588*/                     
    {0x1.f6f9594de8p-5, -0x1.f0f9c4cf885d9p-45}, /* i=589*/                    
    {0x1.fafa6d398p-5, -0x1.024e9d08ce301p-45}, /* i=590*/                     
    {0x1.0180066494p-4, -0x1.7e8c80c5d00f2p-44}, /* i=591*/                    
    {0x1.03824ce1a8p-4, 0x1.10086c6ebde6cp-44}, /* i=592*/                     
    {0x1.078898bc04p-4, 0x1.bf44e72402874p-44}, /* i=593*/                     
    {0x1.098c9ec61cp-4, -0x1.8015cc91ff616p-45}, /* i=594*/                    
    {0x1.0d966cc65p-4, 0x1.f3735158d42c3p-49}, /* i=595*/                      
    {0x1.0f9c356b04p-4, 0x1.c4980f2256fa6p-47}, /* i=596*/                     
    {0x1.13a98bb45p-4, 0x1.f0168af5d72d9p-45}, /* i=597*/                      
    {0x1.15b11a094cp-4, -0x1.e565a88cb0dfdp-44}, /* i=598*/                    
    {0x1.17b94049e8p-4, -0x1.a2febc43331cap-44}, /* i=599*/                    
    {0x1.1bcb55f22p-4, -0x1.9f1d32a2d5372p-45}, /* i=600*/                     
    {0x1.1dd5460c8cp-4, -0x1.d227d61f9e88dp-45}, /* i=601*/                    
    {0x1.21eaf28f58p-4, -0x1.afa5f213c5a5bp-46}, /* i=602*/                    
    {0x1.23f6afac64p-4, -0x1.df6401b93d9b3p-44}, /* i=603*/                    
    {0x1.280ff963cp-4, 0x1.3ee8851dc50dfp-46}, /* i=604*/                      
    {0x1.2a1d86b4ap-4, -0x1.c3b3d2f55224dp-45}, /* i=605*/                     
    {0x1.2e3a740b78p-4, 0x1.d288560689912p-53}, /* i=606*/                     
    {0x1.3049d4c9e4p-4, 0x1.29ff8b0900a32p-44}, /* i=607*/                     
    {0x1.3259d2107cp-4, 0x1.b5474ae667c99p-44}, /* i=608*/                     
    {0x1.367ba3aaap-4, 0x1.882886d8893a8p-44}, /* i=609*/                      
    {0x1.388d78b934p-4, 0x1.0feadd5856604p-44}, /* i=610*/                     
    {0x1.3cb2fd2f68p-4, -0x1.0f6d74f54f373p-48}, /* i=611*/                    
    {0x1.3ec6ad5408p-4, -0x1.e5e3b38ac267ap-46}, /* i=612*/                    
    {0x1.40dafc92e4p-4, -0x1.23b63abb43a32p-45}, /* i=613*/                    
    {0x1.450579dcf8p-4, 0x1.186399c574613p-44}, /* i=614*/                     
    {0x1.3cb2fd2f68p-4, -0x1.0f6d74f54f373p-48}, /* i=611*/                    
    {0x1.3ec6ad5408p-4, -0x1.e5e3b38ac267ap-46}, /* i=612*/                    
    {0x1.4d61fa2514p-4, 0x1.3ffcc075aa95p-45}, /* i=617*/                      
    {0x1.4f7aad9bbcp-4, 0x1.75da8a5871b9ap-45}, /* i=618*/                     
    {0x1.53adfb462cp-4, 0x1.c2c6f11e3581cp-45}, /* i=619*/                     
    {0x1.55c8963e6cp-4, -0x1.2153feab94ebp-44}, /* i=620*/                     
    {0x1.57e3d47c3cp-4, -0x1.085061f7b3786p-44}, /* i=621*/                 
    {0x1.5c1c3c5558p-4, -0x1.0cd9f826e0577p-45}, /* i=622*/                    
    {0x1.5e3966b7e8p-4, 0x1.2951bb9cd2fb7p-44}, /* i=623*/                     
    {0x1.605735ee98p-4, 0x1.7c3cf23a17d9fp-46}, /* i=624*/                     
    {0x1.6494c46ac8p-4, -0x1.1b294ba31c9a5p-44}, /* i=625*/                    
    {0x1.66b4847a68p-4, 0x1.32d4fb541dddp-45}, /* i=626*/                      
    {0x1.68d4eaf26cp-4, 0x1.7ee531d3da9e2p-44}, /* i=627*/                     
    {0x1.6d17acb3e4p-4, 0x1.f5db574a58c15p-44}, /* i=628*/                     
    {0x1.6f3a08ca68p-4, -0x1.b1f6db68cff6cp-45}, /* i=629*/                    
    {0x1.715d0ce368p-4, -0x1.41149840eaa65p-46}, /* i=630*/                    
    {0x1.75a50ebb18p-4, -0x1.db5c368bd3023p-44}, /* i=631*/                    
    {0x1.77ca0d49ccp-4, -0x1.3ef7e60371c2ap-45}, /* i=632*/                    
    {0x1.79efb57b1p-4, -0x1.ff281b9601ce6p-46}, /* i=633*/                     
    {0x1.7e3d04697cp-4, -0x1.1f8744b80ca8fp-45}, /* i=634*/                    
    {0x1.8064abf9b4p-4, -0x1.e1ed7f91288b9p-45}, /* i=635*/                    
    {0x1.828cfed29cp-4, -0x1.deb4fc182476ep-44}, /* i=636*/                    
    {0x1.86dfa808d4p-4, -0x1.2c09bb60238bap-45}, /* i=637*/                    
    {0x1.8909ff3c4cp-4, 0x1.19097bd5ee8f4p-44}, /* i=638*/                     
    {0x1.8b350364c8p-4, -0x1.da8c4bd8546e7p-44}, /* i=639*/                    
    {0x1.8d60b4ee4cp-4, 0x1.900e5cc4e4c82p-44}, /* i=640*/                     
    {0x1.91ba21d6cp-4, -0x1.088de1de62c8cp-44}, /* i=641*/                     
    {0x1.93e7de0fc4p-4, -0x1.80743406505e6p-48}, /* i=642*/                    
    {0x1.9616495e1p-4, -0x1.e17cc57d7b696p-44}, /* i=643*/                     
    {0x1.9a752ef318p-4, -0x1.eb6c2b439cdccp-44}, /* i=644*/                    
    {0x1.9ca5aa1728p-4, 0x1.f44a74b04a5b6p-44}, /* i=645*/                     
    {0x1.9ed6d60b3p-4, 0x1.846ac6badba1cp-46}, /* i=646*/                      
    {0x1.a108b33edcp-4, -0x1.ff522c50af44cp-45}, /* i=647*/                    
    {0x1.a56e8325f4p-4, 0x1.c86eeec5e03ccp-44}, /* i=648*/                     
    {0x1.a7a276badcp-4, 0x1.2c79e9957c9d4p-44}, /* i=649*/                     
    {0x1.a9d71d5258p-4, 0x1.20f04dbb4400ap-46}, /* i=650*/                     
    {0x1.ac0c775e3p-4, -0x1.0d20740fee0f3p-47}, /* i=651*/                     
    {0x1.ae42855098p-4, 0x1.50b57a903b366p-44}, /* i=652*/                     
    {0x1.b2b0beb418p-4, 0x1.ad02ad13bc4d7p-44}, /* i=653*/                     
    {0x1.b4e8eb0bcp-4, -0x1.4e29fd772d38cp-45}, /* i=654*/                     
    {0x1.b721cd1714p-4, 0x1.7e295f660b9dap-44}, /* i=655*/                     
    {0x1.b95b654a78p-4, 0x1.90e5e24764ec7p-45}, /* i=656*/                     
    {0x1.bdd0b9fd08p-4, 0x1.a0fc2eafba507p-44}, /* i=657*/                     
    {0x1.c00c776724p-4, -0x1.a356c78b99edcp-44}, /* i=658*/
    {0x1.c248eccf2p-4, -0x1.1694549c88295p-46}, /* i=659*/
    {0x1.c4861aab94p-4, -0x1.775b6b51fca7bp-46}, /* i=660*/
    {0x1.c6c4017384p-4, -0x1.415cfbccdfea4p-44}, /* i=661*/
    {0x1.cb41fba428p-4, -0x1.7933d3334b1dp-44}, /* i=662*/
    {0x1.cd820ffd28p-4, -0x1.c358f377e27bcp-46}, /* i=663*/
    {0x1.cfc2df223cp-4, 0x1.b2d69192b3939p-44}, /* i=664*/
    {0x1.d204698cb4p-4, 0x1.5e533080ecf32p-47}, /* i=665*/
    {0x1.d446afb64cp-4, 0x1.f95fec4153145p-46}, /* i=666*/
    {0x1.d8cd71303cp-4, -0x1.6d13f77d50b7ap-47}, /* i=667*/
    {0x1.db11ed766cp-4, -0x1.40bcd23c3e44cp-44}, /* i=668*/
    {0x1.dd5727676cp-4, 0x1.2b2443f2d985p-45}, /* i=669*/
    {0x1.df9d1f7f5cp-4, -0x1.31751ca1d17c9p-45}, /* i=670*/
    {0x1.e1e3d63accp-4, -0x1.88c247b543938p-45}, /* i=671*/
    {0x1.e42b4c16ccp-4, -0x1.50d780639590cp-44}, /* i=672*/
    {0x1.e8bc77271cp-4, -0x1.a197240569ddfp-46}, /* i=673*/
    {0x1.eb062d57f4p-4, 0x1.bcf490baf38b3p-45}, /* i=674*/
    {0x1.ed50a4a27p-4, -0x1.50408544a92fap-44}, /* i=675*/
    {0x1.ef9bdd8608p-4, -0x1.e955671ae8b8ap-44}, /* i=676*/
    {0x1.f1e7d882b8p-4, -0x1.765bdaa918999p-44}, /* i=677*/
    {0x1.f4349618fcp-4, -0x1.6e6791f17a5c9p-44}, /* i=678*/
    {0x1.f68216c9ccp-4, 0x1.c9a3bd0891bccp-46}, /* i=679*/
    {0x1.fb1f638184p-4, 0x1.6f3d316ca77fp-44}, /* i=680*/
    {0x1.fd6f308ce4p-4, 0x1.b51fe006d8435p-44}, /* i=681*/
    {0x1.ffbfc2bbc8p-4, -0x1.ff229f20ed3d2p-46}, /* i=682*/
    {0x1.01088d48d6p-3, 0x1.c055a0b1de245p-44}, /* i=683*/
    {0x1.02319c495p-3, -0x1.abb841c89d23p-45}, /* i=684*/
    {0x1.035b0ea194p-3, -0x1.18999d93bfed1p-44}, /* i=685*/
    {0x1.0484e4942ap-3, 0x1.4867cc62a8c08p-44}, /* i=686*/
    {0x1.05af1e63ep-3, 0x1.da36af484664ep-46}, /* i=687*/
    {0x1.0804bea724p-3, -0x1.55d7cb736f965p-45}, /* i=688*/
    {0x1.093025a19ap-3, -0x1.128d0950e065ap-44}, /* i=689*/
    {0x1.0a5bf186fep-3, 0x1.20c7dbf14229ep-45}, /* i=690*/
    {0x1.0b88229b72p-3, -0x1.bb284c008ba7cp-44}, /* i=691*/
    {0x1.0cb4b92356p-3, 0x1.9a02fb2cd8eb5p-47}, /* i=692*/
    {0x1.0de1b56356p-3, 0x1.608adb0ce4227p-44}, /* i=693*/
    {0x1.0f0f17a062p-3, 0x1.a9547d2cbfbb2p-46}, /* i=694*/
    {0x1.103ce01faep-3, 0x1.118edef8bb50ap-46}, /* i=695*/
    {0x1.116b0f26b6p-3, 0x1.eeb24143ef26bp-45}, /* i=696*/
    {0x1.1299a4fb3ep-3, 0x1.82c6326f70b35p-46}, /* i=697*/
    {0x1.13c8a1e35p-3, -0x1.0bfa55f697578p-44}, /* i=698*/
    {0x1.14f806253cp-3, 0x1.f65c144ac2d9fp-46}, /* i=699*/
    {0x1.175805d158p-3, 0x1.f04d633b79054p-45}, /* i=700*/
    {0x1.1888a1c996p-3, -0x1.b3c1301ae9ac8p-45}, /* i=701*/
    {0x1.19b9a637cap-3, 0x1.4a430f4988ed7p-46}, /* i=702*/
    {0x1.1aeb1363b4p-3, 0x1.3219d92f934ccp-45}, /* i=703*/
    {0x1.1c1ce9955cp-3, 0x1.8b891b6d05a73p-48}, /* i=704*/
    {0x1.1d4f291514p-3, -0x1.fed10b114338bp-47}, /* i=705*/
    {0x1.1e81d22b7ap-3, -0x1.e57123f1e6459p-44}, /* i=706*/
    {0x1.1fb4e52174p-3, 0x1.e7216e57e6c06p-48}, /* i=707*/
    {0x1.20e8624038p-3, 0x1.fd946b34ff6ccp-44}, /* i=708*/
    {0x1.221c49d148p-3, -0x1.32a4d18be7546p-49}, /* i=709*/
    {0x1.23509c1e6ep-3, -0x1.b2215ab3ec84fp-45}, /* i=710*/
    {0x1.24855971c4p-3, -0x1.9f1c9091a3411p-44}, /* i=711*/
    {0x1.25ba8215bp-3, -0x1.007bf2d8bf07ap-44}, /* i=712*/
    {0x1.26f01654e6p-3, 0x1.beb0d37066d42p-44}, /* i=713*/
    {0x1.2826167a6cp-3, -0x1.b1d55056642p-46}, /* i=714*/
    {0x1.295c82d19p-3, -0x1.7973c2b40dd9dp-44}, /* i=715*/
    {0x1.2a935ba5f2p-3, -0x1.70e07d84e08ffp-44}, /* i=716*/
    {0x1.2bcaa14382p-3, -0x1.8f4c621480e44p-44}, /* i=717*/
    {0x1.2d0253f67ep-3, 0x1.32ac22596b4ap-45}, /* i=718*/
    {0x1.2e3a740b78p-3, 0x1.d288560689912p-52}, /* i=719*/
    {0x1.2f7301cf4ep-3, 0x1.0f5c70d1a6341p-44}, /* i=720*/
    {0x1.30abfd8f34p-3, -0x1.893a6508d5aa5p-44}, /* i=721*/
    {0x1.31e56798aap-3, -0x1.debad218b9aa2p-44}, /* i=722*/
    {0x1.331f403986p-3, -0x1.ed2dc2cb8ae4ap-44}, /* i=723*/
    {0x1.345987bfeep-3, 0x1.521558148a413p-44}, /* i=724*/
};

/* The following is a degree-6 polynomial generated by Sollya over
   [-0.00202941894531250,0.00212097167968735],
   with absolute error < 2^-70.278.
   The polynomial is P[0]*x + P[1]*x^2 + ... + P[5]*x^6.
   The algorithm assumes that P[0]=1. */
static const double P[6] = {0x1p0,                 /* degree 1 */
                            -0x1.ffffffffffffap-2, /* degree 2 */
                            0x1.555555554f4d8p-2,  /* degree 3 */
                            -0x1.0000000537df6p-2, /* degree 4 */
                            0x1.999a14758b084p-3,  /* degree 5 */
                            -0x1.55362255e0f63p-3, /* degree 6 */
};

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma(a, b, -*hi);
}

// Returns (ah + al) * (bh + bl) - (al * bl)
// We can ignore al * bl when assuming al <= ulp(ah) and bl <= ulp(bh)
static inline void d_mul(double *hi, double *lo, double ah, double al,
                         double bh, double bl) {
  double s, t;

  a_mul(hi, &s, ah, bh);
  t = __builtin_fma(al, bh, s);
  *lo = __builtin_fma(ah, bl, t);
}

/* Given 1 <= x := v.f < 2, where x = v.f, put in h+l a double-double
   approximation of log10(2^e*x), with absolute error bounded by 2^-68.18
   (details below).
*/
static void
cr_log10_fast (double *h, double *l, int e, d64u64 v)
{
  uint64_t m = 0x10000000000000 + (v.u & 0xfffffffffffff);
  /* x = m/2^52 */
  /* if x > sqrt(2), we divide it by 2 to avoid cancellation */
  int c = m >= 0x16a09e667f3bcd;
  e += c; /* now -1074 <= e <= 1024 */
  static const double cy[] = {1.0, 0.5};
  static const uint64_t cm[] = {43, 44};

  int i = m >> cm[c];
  double y = v.f * cy[c];
#define OFFSET 362
  double r = (_INVERSE - OFFSET)[i];
  double l1 = (_LOG10_INV - OFFSET)[i][0];
  double l2 = (_LOG10_INV - OFFSET)[i][1];
  double z = __builtin_fma (r, y, -1.0); /* exact */
  /* evaluate P(z), for |z| < 0.00212097167968735 */
  double ph, pl; /* will hold the approximation of log(1+z) */
  double z2 = z * z; /* |z2| < 4.5e-6 thus the rounding error on z2 is
                        bounded by ulp(4.5e-6) = 2^-70. */
  double p45 = __builtin_fma (P[5], z, P[4]);
  /* |P[5]| < 0.167, |z| < 0.0022, |P[4]| < 0.21 thus |p45| < 0.22:
     the rounding (and total) error on p45 is bounded by ulp(0.22) = 2^-55 */
  double p23 = __builtin_fma (P[3], z, P[2]);
  /* |P[3]| < 0.26, |z| < 0.0022, |P[2]| < 0.34 thus |p23| < 0.35:
     the rounding (and total) error on p23 is bounded by ulp(0.35) = 2^-54 */
  ph = __builtin_fma (p45, z2, p23);
  /* |p45| < 0.22, |z2| < 4.5e-6, |p23| < 0.35 thus |ph| < 0.36:
     the rounding error on ph is bounded by ulp(0.36) = 2^-54.
     Adding the error on p45 multiplied by z2, that on z2 multiplied by p45,
     and that on p23 (ignoring low order errors), we get for the total error
     on ph the following bound:
     2^-54 + err(p45)*4.5e-6 + 0.22*err(z2) + err(p23) <
     2^-54 + 2^-55*4.5e-6 + 0.22*2^-70 + 2^-54 < 2^-52.99 */
  ph = __builtin_fma (ph, z, P[1]);
  /* let ph0 be the value at input, and ph1 the value at output:
     |ph0| < 0.36, |z| < 0.0022, |P[1]| < 0.5 thus |ph1| < 0.501:
     the rounding error on ph1 is bounded by ulp(0.501) = 2^-53.
     Adding the error on ph0 multiplied by z, we get for the total error
     on ph1 the following bound:
     2^-53 + err(ph0)*0.0022 < 2^-53 + 2^-52.99*0.0022 < 2^-52.99 */
  ph *= z2;
  /* let ph2 be the value at output of the above instruction:
     |ph2| < |z2| * |ph1| < 4.5e-6 * 0.501 < 2.26e-6 thus the
     rounding error on ph2 is bounded by ulp(2.26e-6) = 2^-71.
     Adding the error on ph1 multiplied by z2, and the error on z2
     multiplied by ph1, we get for the total error on ph2 the following bound:
     2^-71 + err(ph1)*z2 + ph1*err(z2) <
     2^-71 + 2^-52.99*4.5e-6 + 0.501*2^-70 < 2^-69.32. */

  fast_two_sum (&ph, &pl, z, ph);
  /* let ph2 be the input value, and ph3,pl3 the output values.
     We have |z| < 0.0022 and |ph2| < 2.26e-6, thus |ph3| < 0.0023
     and |pl3| < ulp(ph3) <= 2^-61.
     The rounding error is bounded by 2^-105 |ph3| < 2^-113.76. */

  /* now ph+pl approximates log(1+z), with error bounded by:
     2^-70.278 from the error in the Sollya polynomial
     2^-69.32 for the rounding error on ph2
     2^-113.76 for the rounding error in fast_two_sum()
     This gives a total of 2^-70.278+2^-69.32+2^-113.76 < 2^-68.72. */
  
  /* we divide ph+pl by log(10) */
#define ONE_OVER_LOG10_H 0x1.bcb7b1526e50ep-2
#define ONE_OVER_LOG10_L 0x1.95355baaafad3p-57
  d_mul (&ph, &pl, ph, pl, ONE_OVER_LOG10_H, ONE_OVER_LOG10_L);
  /* Let ph3,pl3 be the input values, and ph4,pl4 the output values.
     The d_mul instruction decomposes into:
     a_mul (ph4, s, ph3, ONE_OVER_log10_H)
     t = __builtin_fma (pl3, ONE_OVER_log10_H, s)
     pl4 = __builtin_fma (ph3, ONE_OVER_log10_L, t)
     where the a_mul() instruction is exact.
     The d_mul() error is bounded by ulp(t)+ulp(pl4)+|pl3|*ONE_OVER_log10_L
     where |pl3|*ONE_OVER_log10_L is the ignored term.
     Since |ph3| < 0.0023 and |ONE_OVER_log10_H| < 1.45, we have |ph4| < 0.0034
     thus |s| < ulp(0.0034) <= 2^-61.
     Then since |pl3| < 2^-61, |ONE_OVER_log10_H| < 1.45, and |s| < 2^-61,
     we have |t| < 2^-61*1.45+2^-61 < 2^-59.7 thus ulp(t) < 2^-112.
     Then since |ph3| < 0.0023, |ONE_OVER_log10_L| < 2^-55.44, |t| < 2^-59.7,
     we have |pl4| < 0.0023*2^-55.44+2^-59.7 < 2^-59.6 thus ulp(pl4) < 2^-112.
     This gives a d_mul() error bounded by:
     ulp(t)+ulp(pl4)+|pl3|*ONE_OVER_log10_L < 2^-112 + 2^-112 + 2^-61*2^-55.44
     < 2^-110.96.

     Adding the total error on ph3+pl3 multiplied by
     ONE_OVER_log10_H+ONE_OVER_log10_L, we get:
     2^-110.96 + 2^-68.72*(ONE_OVER_log10_H+ONE_OVER_log10_L) < 2^-68.19.
  */

  /* We have -1074 <= e <= 1023, l1 is an integer multiple of 2^-42,
     with |l1| <= 2199701051328*2^-42.
     Thus e+l1 is an integer multiple of 2^-42 too, with
     2^42*|e+l1| <= 1074*2^42+2199701051328 < 2^53, thus e+l1 is exactly
     representable below. */

  /* now add e+l1+l2 and ph+pl */
  fast_two_sum (h, l, e + l1, ph);
  /* here |e+l1| < 1074.51 and |ph| < 0.0034 thus |h| < 1074.52 and the
     error of fast_two_sum() is bounded by 2^-105*|h| < 2^-94.93. */
  *l += l2 + pl;
  /* here |l| < ulp(h) <= 2^-42, |l2| < 2^-43 and |pl| < 2^-59.6 thus
     |l_out| < 2^-42+2^-43+2^-59.6 < 2^-41.41, and the rounding error
     is bounded by ulp(l2+pl) + ulp(l) < 2^-95 + 2^-94 < 2^-93.41. */

  /* The absolute error on h + l is bounded by:
     2^-68.19 for the error on ph+pl wrt log10(1+z)
     2^-94.93 for the error in the 2nd fast_two_sum()
     2^-93.41 for the error in *l += l2 + pl
     This gives a total error < 2^-68.19+2^-94.93+2^-93.41 < 2^-68.18.
  */
}

static inline void dint_fromd (dint64_t *a, double b);
static void log_2 (dint64_t *r, dint64_t *x);
static inline double dint_tod (dint64_t *a);

/* accurate path, using Tom Hubrecht's code below */
static double
cr_log10_accurate (double x)
{
  dint64_t X, Y;

#define EXCEPTIONS 27
  static double T[EXCEPTIONS][3] = {
    { 0x1p0, 0x0p0, 0x0p0 },
    { 0x1.fffffffffff7p-1, -0x1.2000000000029p-46, 0x1.fffffffffe1ap-100 },
    { 0x1.fffffffffff5p-1, -0x1.600000000003dp-46, 0x1.fffffffffc88bp-100 },
    { 0x1.fffffffffff3p-1, -0x1.a000000000055p-46, 0x1.fffffffffa475p-100 },
    { 0x1.fffffffffff1p-1, -0x1.e000000000071p-46, 0x1.fffffffff736p-100 },
    { 0x1.ffffffffffffep-1, -0x1.0000000000001p-52, 0x1.fffffffffffffp-106 },
    { 0x1.fffffffffff6p-1, -0x1.4000000000032p-46, -0x1.4d555555555a3p-139 },
    { 0x1.fffffffffffp-1, -0x1.000000000004p-45, -0x1.55555555555d5p-137 },
    { 0x1.ffffffffffeep-1, -0x1.2000000000051p-45, -0x1.e6000000000cdp-137 },
    { 0x1.fffffffffff4p-1, -0x1.8000000000048p-46, -0x1.2000000000051p-138 },
    { 0x1.fffffffffff2p-1, -0x1.c000000000062p-46, -0x1.c9555555555ebp-138 },
    { 0x1.ffffffffffeap-1, -0x1.6000000000079p-45, -0x1.bbaaaaaaaab8fp-136 },
    { 0x1.ffffffffffe8p-1, -0x1.800000000009p-45, -0x1.20000000000a2p-135 },
    { 0x1.ffffffffffff8p-1, -0x1.0000000000002p-50, -0x1.5555555555559p-152 },
    { 0x1.ffffffffffffcp-1, -0x1.0000000000001p-51, -0x1.5555555555557p-155 },
    { 0x1.fffffffffffc8p-1, -0x1.c000000000019p-48, 0x1.ffffffffff8dbp-102 },
    { 0x1.fffffffffffd8p-1, -0x1.400000000000dp-48, 0x1.ffffffffffd65p-102 },
    { 0x1.ffffffffffff4p-1, -0x1.8000000000005p-50, 0x1.fffffffffffb8p-104 },
    { 0x1.fffffffffff9p-1, -0x1.c000000000031p-47, -0x1.c9555555555ap-141 },
    { 0x1.fffffffffffap-1, -0x1.8000000000024p-47, -0x1.2000000000029p-141 },
    { 0x1.fffffffffffbp-1, -0x1.4000000000019p-47, -0x1.4d5555555557cp-142 },
    { 0x1.fffffffffffcp-1, -0x1.000000000001p-47, -0x1.5555555555575p-143 },
    { 0x1.fffffffffffdp-1, -0x1.8000000000012p-48, -0x1.2000000000014p-144 },
    { 0x1.fffffffffffep-1, -0x1.0000000000008p-48, -0x1.5555555555565p-146 },
    { 0x1.fffffffffffe8p-1, -0x1.8000000000009p-49, -0x1.200000000000ap-147 },
    { 0x1.ffffffffffffp-1, -0x1.0000000000004p-49, -0x1.555555555555dp-149 },
    { 0x1.fffffffffff8p-1, -0x1.000000000002p-46, -0x1.5555555555595p-140 },
  };
  for (int i = 0; i < EXCEPTIONS; i++)
    if (x == T[i][0])
      return T[i][1] + T[i][2];

  dint_fromd (&X, x);
  /* x = (-1)^sgn*2^ex*(hi/2^63+lo/2^127) */
  log_2 (&Y, &X);
  return dint_tod (&Y);
}

double
cr_log10 (double x)
{
  d64u64 v = {.f = x};
  int e = (v.u >> 52) - 0x3ff;
  if (e >= 0x400 || e == -0x3ff) /* x <= 0 or NaN/Inf or subnormal */
  {
    if (x <= 0.0)
    {
      /* log10(x<0) is NaN, f(+/-0) is -Inf and raises DivByZero */
      if (x < 0)
        return 0.0 / 0.0;
      else
        return 1.0 / -0.0;
    }
    if (e == 0x400) /* +Inf or NaN */
      return x;
    if (e == -0x3ff) /* subnormal */
    {
      v.f *= 0x1p52;
      e = (v.u >> 52) - 0x3ff - 52;
    }
  }
  /* now x > 0 */
  /* normalize v in [1,2) */
  v.u = (0x3fful << 52) | (v.u & 0xfffffffffffff);
  /* now x = m*2^e with 1 <= m < 2 (m = v.f) and -1074 <= e <= 1023 */
  double h, l;
  cr_log10_fast (&h, &l, e, v);

  // if (x == TRACE) printf ("h=%la l=%la\n", h, l);

  /* maximal absolute error from cr_log10_fast: 2^-68.18 < 1.c4p-69 */
  static double err = 0x1.c4p-69;

  double left = h + (l - err), right = h + (l + err);
  if (left == right)
  {
    // if (x == TRACE) printf ("rounding test succedded\n");
    return left;
  }
  // if (x == TRACE) printf ("rounding test failed\n");
  /* the probability of failure of the fast path is about 2^-11.5 */
  return 0;
  return cr_log10_accurate (x);
}

/* the following code was copied from Tom Hubrecht's implementation of
   correctly rounded pow for CORE-MATH */

// Approximation for the second iteration
static inline void p_2(dint64_t *r, dint64_t *z) {
  cp_dint(r, &P_2[0]);

  mul_dint(r, z, r);
  add_dint(r, &P_2[1], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[2], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[3], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[4], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[5], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[6], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[7], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[8], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[9], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[10], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[11], r);

  mul_dint(r, z, r);
  add_dint(r, &P_2[12], r);

  mul_dint(r, z, r);
}

static void log_2(dint64_t *r, dint64_t *x) {
#if DEBUG > 0
  printf("Calcul du logarithme :\n");
  printf("  x := ");
  print_dint(x);
  printf("\n");
#endif

  int64_t E = x->ex;

  // Find the lookup index
  uint16_t i = x->hi >> 55;

  if (x->hi > 0xb504f333f9de6484) {
    E++;
    i = i >> 1;
  }

  x->ex = x->ex - E;

#if DEBUG > 0
  printf("  E := %ld\n\n", E);
#endif

  dint64_t z;
  mul_dint(&z, x, &_INVERSE_2[i - 128]);

#if DEBUG > 0
  printf("  y := ");
  print_dint(x);
  printf("  i := %d\n", i);
  printf("  r_i := ");
  print_dint(&_INVERSE_2[i - 128]);
  printf("  y·r_i := ");
  print_dint(&z);
  printf("\n");
#endif

  add_dint(&z, &M_ONE, &z);

#if DEBUG > 0
  printf("  z := ");
  print_dint(&z);
  printf("\n");
#endif

  // E·log10(2)
  mul_dint_2(r, E, &LOG2);

#if DEBUG > 0
  printf("  E·log(2) := ");
  print_dint(r);
  printf("\n");
#endif

#if DEBUG > 0
  printf("  -log(r_i) := ");
  print_dint(&_LOG_INV_2[i - 128]);
  printf("  E·log(2) - log(r_i) := ");
  print_dint(r);
  printf("\n");
#endif

  dint64_t p;

  p_2(&p, &z);

  add_dint(&p, &_LOG_INV_2[i - 128], &p);

#if DEBUG > 0
  printf("  log(1 + z) := ");
  print_dint(&p);
  printf("\n");
#endif

  add_dint(r, &p, r);

#if DEBUG > 0
  printf("  log(x) := ");
  print_dint(r);
  printf("\n");
#endif
}

typedef union {
  double f;
  uint64_t u;
} f64_u;

// Extract both the mantissa and exponent of a double
static inline void fast_extract(int64_t *e, uint64_t *m, double x) {
  f64_u _x = {.f = x};

  *e = (_x.u >> 52) & 0x7ff;
  *m = (_x.u & (~0ul >> 12)) + (*e ? (1ul << 52) : 0);
  *e = *e - 0x3ff;
}

// Convert a double to the corresponding dint64_t value
static inline void dint_fromd(dint64_t *a, double b) {
  fast_extract(&a->ex, &a->hi, b);

  uint32_t t = __builtin_clzl(a->hi);

  a->sgn = b < 0.0;
  a->hi = a->hi << t;
  a->ex = a->ex - (t > 11 ? t - 12 : 0);
  a->lo = 0;
}

// Convert a dint64_t value to a double
// assuming the input is not in the subnormal range
static inline double dint_tod(dint64_t *a) {

  f64_u r = {.u = (a->hi >> 11) | (0x3ffl << 52)};
  /* r contains the upper 53 bits of a->hi, 1 <= r < 2 */

  double rd = 0.0;
  /* if round bit is 1, add 2^-53 */
  if ((a->hi >> 10) & 0x1)
    rd += 0x1p-53;

  /* if trailing bits after the rounding bit are non zero, add 2^-54 */
  if (a->hi & 0x3ff || a->lo)
    rd += 0x1p-54;

  r.u = r.u | a->sgn << 63;
  r.f += (a->sgn == 0) ? rd : -rd;

  f64_u e;

  /* For log, the result is always in the normal range,
     thus a->ex > -1023. Similarly, we cannot have a->ex > 1023. */

  e.u = ((a->ex + 1023) & 0x7ff) << 52;

  return r.f * e.f;
}
