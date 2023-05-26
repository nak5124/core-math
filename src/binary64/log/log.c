/* Correctly rounded logarithm of binary64 values.

Copyright (c) 2022 INRIA and CERN.
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

// Add a + (bh + bl), assuming |a| >= |bh|
// Code copied from pow.h
static inline void fast_sum(double *hi, double *lo, double a, double bh,
                            double bl) {
  fast_two_sum(hi, lo, a, bh);
  /* |(a+bh)-(hi+lo)| <= 2^-105 |hi| and |lo| < ulp(hi) */
  *lo += bl;
  /* |(a+bh+bl)-(hi+lo)| <= 2^-105 |hi| + ulp(lo),
     where |lo| <= ulp(hi) + |bl|. */
}

// Multiply exactly a and b, such that *hi + *lo = a * b.
// Code copied from pow.h
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma (a, b, -*hi);
}

/* for 181 <= i <= 362, r[i] = _INVERSE[i-181] is a 9-bit approximation of
   1/x[i], where i*2^-8 <= x[i] < (i+1)*2^-8.
   More precisely r[i] is a 9-bit value such that r[i]*y-1 is representable
   exactly on 53 bits for for any y, i*2^-8 <= y < (i+1)*2^-8.
   Moreover |r[i]*y-1| < 0.0040283203125.
   Table generated with the accompanying pow.sage file,
   with l=inverse_centered(k=8,prec=9,maxbits=53,verbose=false).
   Coied from pow.h. */
static const double _INVERSE[182]= {
    0x1.69p+0, 0x1.67p+0, 0x1.65p+0, 0x1.63p+0, 0x1.61p+0, 0x1.5fp+0, 0x1.5ep+0,
    0x1.5cp+0, 0x1.5ap+0, 0x1.58p+0, 0x1.56p+0, 0x1.54p+0, 0x1.53p+0, 0x1.51p+0,
    0x1.4fp+0, 0x1.4ep+0, 0x1.4cp+0, 0x1.4ap+0, 0x1.48p+0, 0x1.47p+0, 0x1.45p+0,
    0x1.44p+0, 0x1.42p+0, 0x1.4p+0, 0x1.3fp+0, 0x1.3dp+0, 0x1.3cp+0, 0x1.3ap+0,
    0x1.39p+0, 0x1.37p+0, 0x1.36p+0, 0x1.34p+0, 0x1.33p+0, 0x1.32p+0, 0x1.3p+0,
    0x1.2fp+0, 0x1.2dp+0, 0x1.2cp+0, 0x1.2bp+0, 0x1.29p+0, 0x1.28p+0, 0x1.27p+0,
    0x1.25p+0, 0x1.24p+0, 0x1.23p+0, 0x1.21p+0, 0x1.2p+0, 0x1.1fp+0, 0x1.1ep+0,
    0x1.1cp+0, 0x1.1bp+0, 0x1.1ap+0, 0x1.19p+0, 0x1.17p+0, 0x1.16p+0, 0x1.15p+0,
    0x1.14p+0, 0x1.13p+0, 0x1.12p+0, 0x1.1p+0, 0x1.0fp+0, 0x1.0ep+0, 0x1.0dp+0,
    0x1.0cp+0, 0x1.0bp+0, 0x1.0ap+0, 0x1.09p+0, 0x1.08p+0, 0x1.07p+0, 0x1.06p+0,
    0x1.05p+0, 0x1.04p+0, 0x1.03p+0, 0x1.02p+0, 0x1.00p+0, 0x1.00p+0, 0x1.fdp-1,
    0x1.fbp-1, 0x1.f9p-1, 0x1.f7p-1, 0x1.f5p-1, 0x1.f3p-1, 0x1.f1p-1, 0x1.fp-1,
    0x1.eep-1, 0x1.ecp-1, 0x1.eap-1, 0x1.e8p-1, 0x1.e6p-1, 0x1.e5p-1, 0x1.e3p-1,
    0x1.e1p-1, 0x1.dfp-1, 0x1.ddp-1, 0x1.dcp-1, 0x1.dap-1, 0x1.d8p-1, 0x1.d7p-1,
    0x1.d5p-1, 0x1.d3p-1, 0x1.d2p-1, 0x1.dp-1, 0x1.cep-1, 0x1.cdp-1, 0x1.cbp-1,
    0x1.c9p-1, 0x1.c8p-1, 0x1.c6p-1, 0x1.c5p-1, 0x1.c3p-1, 0x1.c2p-1, 0x1.cp-1,
    0x1.bfp-1, 0x1.bdp-1, 0x1.bcp-1, 0x1.bap-1, 0x1.b9p-1, 0x1.b7p-1, 0x1.b6p-1,
    0x1.b4p-1, 0x1.b3p-1, 0x1.b1p-1, 0x1.bp-1, 0x1.aep-1, 0x1.adp-1, 0x1.acp-1,
    0x1.aap-1, 0x1.a9p-1, 0x1.a7p-1, 0x1.a6p-1, 0x1.a5p-1, 0x1.a3p-1, 0x1.a2p-1,
    0x1.a1p-1, 0x1.9fp-1, 0x1.9ep-1, 0x1.9dp-1, 0x1.9cp-1, 0x1.9ap-1, 0x1.99p-1,
    0x1.98p-1, 0x1.96p-1, 0x1.95p-1, 0x1.94p-1, 0x1.93p-1, 0x1.91p-1, 0x1.9p-1,
    0x1.8fp-1, 0x1.8ep-1, 0x1.8dp-1, 0x1.8bp-1, 0x1.8ap-1, 0x1.89p-1, 0x1.88p-1,
    0x1.87p-1, 0x1.86p-1, 0x1.84p-1, 0x1.83p-1, 0x1.82p-1, 0x1.81p-1, 0x1.8p-1,
    0x1.7fp-1, 0x1.7ep-1, 0x1.7cp-1, 0x1.7bp-1, 0x1.7ap-1, 0x1.79p-1, 0x1.78p-1,
    0x1.77p-1, 0x1.76p-1, 0x1.75p-1, 0x1.74p-1, 0x1.73p-1, 0x1.72p-1, 0x1.71p-1,
    0x1.7p-1, 0x1.6fp-1, 0x1.6ep-1, 0x1.6dp-1, 0x1.6cp-1, 0x1.6bp-1, 0x1.6ap-1,
};

/* For 181 <= i <= 362, (h,l) = _LOG_INV[i-181] is a double-double nearest
   approximation of -log(r) for r=_INVERSE[i-181], h being an integer
   multiple of 2^-42.
   Since |l| < 2^-43, the maximal error is 1/2 ulp(l) <= 2^-97.
   Copied from pow.h. */
static const double _LOG_INV[182][2] = {
    {-0x1.5ff3070a79p-2, -0x1.e9e439f105039p-45},
    {-0x1.5a42ab0f4dp-2, 0x1.e63af2df7ba69p-50},
    {-0x1.548a2c3addp-2, -0x1.3167e63081cf7p-45},
    {-0x1.4ec97326p-2, -0x1.34d7aaf04d104p-45},
    {-0x1.4900680401p-2, 0x1.8bccffe1a0f8cp-44},
    {-0x1.432ef2a04fp-2, 0x1.fb129931715adp-44},
    {-0x1.404308686ap-2, -0x1.f8ef43049f7d3p-44},
    {-0x1.3a64c55694p-2, -0x1.7a71cbcd735dp-44},
    {-0x1.347dd9a988p-2, 0x1.5594dd4c58092p-45},
    {-0x1.2e8e2bae12p-2, 0x1.67b1e99b72bd8p-45},
    {-0x1.2895a13de8p-2, -0x1.a8d7ad24c13fp-44},
    {-0x1.22941fbcf8p-2, 0x1.a6976f5eb0963p-44},
    {-0x1.1f8ff9e48ap-2, -0x1.7946c040cbe77p-45},
    {-0x1.1980d2dd42p-2, -0x1.b7b3a7a361c9ap-45},
    {-0x1.136870293bp-2, 0x1.d3e8499d67123p-44},
    {-0x1.1058bf9ae5p-2, 0x1.4ab9d817d52cdp-44},
    {-0x1.0a324e2739p-2, -0x1.c6bee7ef4030ep-47},
    {-0x1.0402594b4dp-2, -0x1.036b89ef42d7fp-48},
    {-0x1.fb9186d5e4p-3, 0x1.d572aab993c87p-47},
    {-0x1.f550a564b8p-3, 0x1.323e3a09202fep-45},
    {-0x1.e8c0252aa6p-3, 0x1.6805b80e8e6ffp-45},
    {-0x1.e27076e2bp-3, 0x1.a342c2af0003cp-44},
    {-0x1.d5c216b4fcp-3, 0x1.1ba91bbca681bp-45},
    {-0x1.c8ff7c79aap-3, 0x1.7794f689f8434p-45},
    {-0x1.c2968558c2p-3, 0x1.cfd73dee38a4p-45},
    {-0x1.b5b519e8fcp-3, 0x1.4b722ec011f31p-44},
    {-0x1.af3c94e80cp-3, 0x1.a4e633fcd9066p-52},
    {-0x1.a23bc1fe2cp-3, 0x1.539cd91dc9f0bp-44},
    {-0x1.9bb362e7ep-3, 0x1.1f2a8a1ce0ffcp-45},
    {-0x1.8e928de886p-3, -0x1.a8154b13d72d5p-44},
    {-0x1.87fa06520cp-3, -0x1.22120401202fcp-44},
    {-0x1.7ab890210ep-3, 0x1.bdb9072534a58p-45},
    {-0x1.740f8f5404p-3, 0x1.0b66c99018aa1p-44},
    {-0x1.6d60fe719ep-3, 0x1.bc6e557134767p-44},
    {-0x1.5ff3070a7ap-3, 0x1.8586f183bebf2p-44},
    {-0x1.59338d9982p-3, -0x1.0ba68b7555d4ap-48},
    {-0x1.4ba36f39a6p-3, 0x1.4354bb3f219e5p-44},
    {-0x1.44d2b6ccb8p-3, 0x1.70cc16135783cp-46},
    {-0x1.3dfc2b0eccp-3, -0x1.8a72a62b8c13fp-45},
    {-0x1.303d718e48p-3, 0x1.680b5ce3ecb05p-50},
    {-0x1.29552f82p-3, 0x1.5b967f4471dfcp-44},
    {-0x1.2266f190a6p-3, 0x1.4d20ab840e7f6p-45},
    {-0x1.1478584674p-3, -0x1.563451027c75p-46},
    {-0x1.0d77e7cd08p-3, -0x1.cb2cd2ee2f482p-44},
    {-0x1.0671512ca6p-3, 0x1.a47579cdc0a3dp-45},
    {-0x1.f0a30c0118p-4, 0x1.d599e83368e91p-44},
    {-0x1.e27076e2bp-4, 0x1.a342c2af0003cp-45},
    {-0x1.d4313d66ccp-4, 0x1.9454379135713p-45},
    {-0x1.c5e548f5bcp-4, -0x1.d0c57585fbe06p-46},
    {-0x1.a926d3a4acp-4, -0x1.563650bd22a9cp-44},
    {-0x1.9ab4246204p-4, 0x1.8a64826787061p-45},
    {-0x1.8c345d6318p-4, -0x1.b20f5acb42a66p-44},
    {-0x1.7da766d7bp-4, -0x1.2cc844480c89bp-44},
    {-0x1.60658a9374p-4, -0x1.0c3b1dee9c4f8p-44},
    {-0x1.51b073f06p-4, -0x1.83f69278e686ap-44},
    {-0x1.42edcbea64p-4, -0x1.bc0eeea7c9acdp-46},
    {-0x1.341d7961bcp-4, -0x1.1d0929983761p-44},
    {-0x1.253f62f0ap-4, -0x1.416f8fb69a701p-44},
    {-0x1.16536eea38p-4, 0x1.47c5e768fa309p-46},
    {-0x1.f0a30c0118p-5, 0x1.d599e83368e91p-45},
    {-0x1.d276b8adbp-5, -0x1.6a423c78a64bp-46},
    {-0x1.b42dd71198p-5, 0x1.c827ae5d6704cp-46},
    {-0x1.95c830ec9p-5, 0x1.c148297c5feb8p-45},
    {-0x1.77458f633p-5, 0x1.181dce586af09p-44},
    {-0x1.58a5bafc9p-5, 0x1.b2b739570ad39p-45},
    {-0x1.39e87b9fe8p-5, -0x1.eafd480ad9015p-44},
    {-0x1.1b0d98924p-5, 0x1.3401e9ae889bbp-44},
    {-0x1.f829b0e78p-6, -0x1.980267c7e09e4p-45},
    {-0x1.b9fc027bp-6, 0x1.b9a010ae6922ap-44},
    {-0x1.7b91b07d6p-6, 0x1.3b955b602ace4p-44},
    {-0x1.3cea44347p-6, 0x1.6a2c432d6a40bp-44},
    {-0x1.fc0a8b0fcp-7, -0x1.f1e7cf6d3a69cp-50},
    {-0x1.7dc475f82p-7, 0x1.eb1245b5da1f5p-44},
    {-0x1.fe02a6b1p-8, -0x1.9e23f0dda40e4p-46},
    {0, 0},
    {0, 0},
    {0x1.812121458p-8, 0x1.ad50382973f27p-46},
    {0x1.41929f968p-7, 0x1.977c755d01368p-46},
    {0x1.c317384c8p-7, -0x1.41f33fcefb9fep-44},
    {0x1.228fb1feap-6, 0x1.713e3284991fep-45},
    {0x1.63d617869p-6, 0x1.7abf389596542p-47},
    {0x1.a55f548c6p-6, -0x1.de0709f2d03c9p-45},
    {0x1.e72bf2814p-6, -0x1.8d75149774d47p-45},
    {0x1.0415d89e78p-5, -0x1.dddc7f461c516p-44},
    {0x1.252f32f8dp-5, 0x1.83e9ae021b67bp-45},
    {0x1.466aed42ep-5, -0x1.c167375bdfd28p-45},
    {0x1.67c94f2d48p-5, 0x1.dac20827cca0cp-44},
    {0x1.894aa149f8p-5, 0x1.9a19a8be97661p-44},
    {0x1.aaef2d0fbp-5, 0x1.0fc1a353bb42ep-45},
    {0x1.bbcebfc69p-5, -0x1.7bf868c317c2ap-46},
    {0x1.dda8adc68p-5, -0x1.1b1ac64d9e42fp-45},
    {0x1.ffa6911ab8p-5, 0x1.3008c98381a8fp-45},
    {0x1.10e45b3cbp-4, -0x1.7cf69284a3465p-44},
    {0x1.2207b5c784p-4, 0x1.49d8cfc10c7bfp-44},
    {0x1.2aa04a447p-4, 0x1.7a48ba8b1cb41p-44},
    {0x1.3bdf5a7d2p-4, -0x1.19bd0ad125895p-44},
    {0x1.4d3115d208p-4, -0x1.53a2582f4e1efp-48},
    {0x1.55e10050ep-4, 0x1.c1d740c53c72ep-47},
    {0x1.674f089364p-4, 0x1.a79994c9d3302p-44},
    {0x1.78d02263d8p-4, 0x1.69b5794b69fb7p-47},
    {0x1.8197e2f41p-4, -0x1.c0fe460d20041p-44},
    {0x1.9335e5d594p-4, 0x1.3115c3abd47dap-45},
    {0x1.a4e7640b1cp-4, -0x1.e42b6b94407c8p-47},
    {0x1.adc77ee5bp-4, -0x1.573b209c31904p-44},
    {0x1.bf968769fcp-4, 0x1.4218c8d824283p-45},
    {0x1.d179788218p-4, 0x1.36433b5efbeedp-44},
    {0x1.da72763844p-4, 0x1.a89401fa71733p-46},
    {0x1.ec739830ap-4, 0x1.11fcba80cdd1p-44},
    {0x1.f57bc7d9p-4, 0x1.76a6c9ea8b04ep-46},
    {0x1.03cdc0a51ep-3, 0x1.81a9cf169fc5cp-44},
    {0x1.08598b59e4p-3, -0x1.7e5dd7009902cp-45},
    {0x1.1178e8227ep-3, 0x1.1ef78ce2d07f2p-45},
    {0x1.160c8024b2p-3, 0x1.ec2d2a9009e3dp-45},
    {0x1.1f3b925f26p-3, -0x1.5f74e9b083633p-46},
    {0x1.23d712a49cp-3, 0x1.00d238fd3df5cp-46},
    {0x1.2d1610c868p-3, 0x1.39d6ccb81b4a1p-47},
    {0x1.31b994d3a4p-3, 0x1.f098ee3a5081p-44},
    {0x1.3b08b6758p-3, -0x1.aade8f29320fbp-44},
    {0x1.3fb45a5992p-3, 0x1.19713c0cae559p-44},
    {0x1.4913d8333cp-3, -0x1.53e43558124c4p-44},
    {0x1.4dc7b897bcp-3, 0x1.c79b60ae1ff0fp-47},
    {0x1.5737cc9018p-3, 0x1.9baa7a6b887f6p-44},
    {0x1.5bf406b544p-3, -0x1.27023eb68981cp-46},
    {0x1.6574ebe8c2p-3, -0x1.98c1d34f0f462p-44},
    {0x1.6a399dabbep-3, -0x1.8f934e66a15a6p-44},
    {0x1.6f0128b756p-3, 0x1.577390d31ef0fp-44},
    {0x1.7898d85444p-3, 0x1.8e67be3dbaf3fp-44},
    {0x1.7d6903caf6p-3, -0x1.4c06b17c301d7p-45},
    {0x1.871213750ep-3, 0x1.328eb42f9af75p-44},
    {0x1.8beafeb39p-3, -0x1.73d54aae92cd1p-47},
    {0x1.90c6db9fccp-3, -0x1.935f57718d7cap-46},
    {0x1.9a8778debap-3, 0x1.470fa3efec39p-44},
    {0x1.9f6c40708ap-3, -0x1.337d94bcd3f43p-44},
    {0x1.a454082e6ap-3, 0x1.60a77c81f7171p-44},
    {0x1.ae2ca6f672p-3, 0x1.7a8d5ae54f55p-44},
    {0x1.b31d8575bcp-3, 0x1.c794e562a63cbp-44},
    {0x1.b811730b82p-3, 0x1.e90683b9cd768p-46},
    {0x1.bd087383bep-3, -0x1.d4bc4595412b6p-45},
    {0x1.c6ffbc6fp-3, 0x1.ee138d3a69d43p-44},
    {0x1.cc000c9db4p-3, -0x1.d6d585d57aff9p-46},
    {0x1.d1037f2656p-3, -0x1.84a7e75b6f6e4p-47},
    {0x1.db13db0d48p-3, 0x1.2806a847527e6p-44},
    {0x1.e020cc6236p-3, -0x1.52b00adb91424p-45},
    {0x1.e530effe72p-3, -0x1.fdbdbb13f7c18p-44},
    {0x1.ea4449f04ap-3, 0x1.5e91663732a36p-44},
    {0x1.f474b134ep-3, -0x1.bae49f1df7b5ep-44},
    {0x1.f991c6cb3cp-3, -0x1.90d04cd7cc834p-44},
    {0x1.feb2233eap-3, 0x1.f3418de00938bp-45},
    {0x1.01eae5626cp-2, 0x1.a43dcfade85aep-44},
    {0x1.047e60cde8p-2, 0x1.dbdf10d397f3cp-45},
    {0x1.09aa572e6cp-2, 0x1.b50a1e1734342p-44},
    {0x1.0c42d67616p-2, 0x1.7188b163ceae9p-45},
    {0x1.0edd060b78p-2, 0x1.019b52d8435f5p-47},
    {0x1.1178e8227ep-2, 0x1.1ef78ce2d07f2p-44},
    {0x1.14167ef367p-2, 0x1.e0c07824daaf5p-44},
    {0x1.16b5ccbadp-2, -0x1.23299042d74bfp-44},
    {0x1.1bf99635a7p-2, -0x1.1ac89575c2125p-44},
    {0x1.1e9e16788ap-2, -0x1.82eaed3c8b65ep-44},
    {0x1.214456d0ecp-2, -0x1.caf0428b728a3p-44},
    {0x1.23ec5991ecp-2, -0x1.6dbe448a2e522p-44},
    {0x1.269621134ep-2, -0x1.1b61f10522625p-44},
    {0x1.2941afb187p-2, -0x1.210c2b730e28bp-44},
    {0x1.2bef07cdc9p-2, 0x1.a9cfa4a5004f4p-45},
    {0x1.314f1e1d36p-2, -0x1.8e27ad3213cb8p-45},
    {0x1.3401e12aedp-2, -0x1.17c73556e291dp-44},
    {0x1.36b6776be1p-2, 0x1.16ecdb0f177c8p-46},
    {0x1.396ce359bcp-2, -0x1.5839c5663663dp-47},
    {0x1.3c25277333p-2, 0x1.83b54b606bd5cp-46},
    {0x1.3edf463c17p-2, -0x1.f067c297f2c3fp-44},
    {0x1.419b423d5fp-2, -0x1.ce379226de3ecp-44},
    {0x1.44591e053ap-2, -0x1.6e95892923d88p-47},
    {0x1.4718dc271cp-2, 0x1.06c18fb4c14c5p-44},
    {0x1.49da7f3bccp-2, 0x1.07b334daf4b9ap-44},
    {0x1.4c9e09e173p-2, -0x1.e20891b0ad8a4p-45},
    {0x1.4f637ebbaap-2, -0x1.fc158cb3124b9p-44},
    {0x1.522ae0738ap-2, 0x1.ebe708164c759p-45},
    {0x1.54f431b7bep-2, 0x1.a8954c0910952p-46},
    {0x1.57bf753c8dp-2, 0x1.fadedee5d40efp-46},
    {0x1.5a8cadbbeep-2, -0x1.7c79b0af7ecf8p-48},
    {0x1.5d5bddf596p-2, -0x1.a0b2a08a465dcp-47},
    {0x1.602d08af09p-2, 0x1.ebe9176df3f65p-46},
    {0x1.630030b3abp-2, -0x1.db623e731aep-45},
};

/* The following is a degree-8 polynomial generated by Sollya for
   log(1+x)-x+x^2/2 over [-0.0040283203125,0.0040283203125]
   with absolute error < 2^-81.63
   and relative error < 2^-72.423 (see P_1.sollya).
   The relative error is for x - x^2/2 + P(x) with respect to log(1+x).
   Copied from pow.h. */
static const double P_1[] = {0x1.5555555555558p-2,  /* degree 3 */
                             -0x1.0000000000003p-2, /* degree 4 */
                             0x1.999999981f535p-3,  /* degree 5 */
                             -0x1.55555553d1eb4p-3, /* degree 6 */
                             0x1.2494526fd4a06p-3,  /* degree 7 */
                             -0x1.0001f0c80e8cep-3, /* degree 8 */
};

/*
  Put in hi+lo an approximation of log(1 + z) - z,
  for |z| < 0.0040283203125, and z integer multiple of 2^-61,
  with maximal error 2^-81.63 from the Sollya polynomial.
  Let p2,...,p8 be the polynomial coefficients. Since |p3*z^3| < 2^-25.45,
  it suffices to compute p3*z^3+...+p8*z^8 in double precision, and we only
  need a double-double representation for p2*z^2.
  Note: we impose the degree-2 coefficient to be -1/2, which saves a
  double-double multiplication.

  Maximal absolute error: |hi + lo - (log(1+z) - z)| < 2^-75.492
  with |hi| < 8.11369e-6 and |lo| < 2.2e-8 (see pow.tex).

  Maximal relative error 2^-67.55 from analyze_p1_all(k=8,rel=true),
  taking into account the z term (which is not computed here), and
  assuming there is no rounding error in z + p_1(z).

  Code copied from pow.c, see comments in pow.c.
*/
static inline void p_1 (double *hi, double *lo, double z) {
  double wh, wl;
  a_mul (&wh, &wl, z, z);
  double t = __builtin_fma (P_1[5], z, P_1[4]);
  double u = __builtin_fma (P_1[3], z, P_1[2]);
  double v = __builtin_fma (P_1[1], z, P_1[0]);
  u = __builtin_fma (t, wh, u);
  v = __builtin_fma (u, wh, v);
  u = v * wh;
  *hi = -0.5 * wh;
  *lo = __builtin_fma (u, z, -0.5 * wl);
}

/* Approximation of log|x|, assuming x is not zero:
   log|x| is approximated by hi + lo.

  For E <> 0, i.e., for x outside [0x1.6a09e667f3bcdp-1, 0x1.6a09e667f3bcdp+0).
  the relative error is bounded by: |hi + lo - log|x||/|log|x|| < 2^-73.528
  (this bound is maybe not tight, the largest relative error we observe is for
  x=0x1.77f73e3aefee7p+0 with 2^-75.74).

  Return 0 if E<>0, and 1 if E=0.

  If E<>0, i.e., x < sqrt(2)/2 or sqrt(2) < x, then the relative error
  is bounded by 2^-73.528, and |lo/hi| < 2^-23.9.

  if E=0, i.e., sqrt(2)/2 < x < sqrt(2), then the relative error is bounded
  by 2^-67.052, and |lo/hi| < 2^-52.

  [Code copied from pow.c, see comments in pow.c]

  Assumes 1 <= x < 2, where x = v.f, and the considered input is 2^_e*x.
*/
static int cr_log_fast (double *hi, double *lo, int _e, d64u64 v) {
  double t = v.f;

  uint64_t i;

  uint64_t _m = 0x10000000000000 + (v.u & 0xfffffffffffff);
  uint64_t c = _m >= 0x16a09e667f3bcd;
  static const double cy[] = {1.0, 0.5};
  static const uint64_t cm[] = {44, 45};

  _e += c;
  double E = _e;
  i = _m >> cm[c];
  t *= cy[c];

  double r = (_INVERSE - 181)[i];
  double l1 = (_LOG_INV - 181)[i][0];
  double l2 = (_LOG_INV - 181)[i][1];

  double z = __builtin_fma (r, t, -1.0); /* exact */

#define LOG2_H 0x1.62e42fefa38p-1
#define LOG2_L 0x1.ef35793c7673p-45

  double th, tl;
  th = __builtin_fma (E, LOG2_H, l1);
  tl = __builtin_fma (E, LOG2_L, l2);

  fast_sum (hi, lo, th, z, tl);
  double ph, pl;
  p_1 (&ph, &pl, z);
  fast_sum (hi, lo, *hi, ph, *lo + pl);
  if (_e == 0)
  {
    fast_two_sum (hi, lo, *hi, *lo);
    return 1;
  }

  return 0;
}

static inline void dint_fromd (dint64_t *a, double b);
static void log_2 (dint64_t *r, dint64_t *x);
static inline double dint_tod (dint64_t *a);

/* accurate path, using Tom Hubrecht's code below */
static double
cr_log_accurate (double x)
{
  dint64_t X, Y;

  if (x == 1.0)
    return 0.0;

  dint_fromd (&X, x);
  /* x = (-1)^sgn*2^ex*(hi/2^63+lo/2^127) */
  log_2 (&Y, &X);
  return dint_tod (&Y);
}

double
cr_log (double x)
{
  d64u64 v = {.f = x};
  int e = (v.u >> 52) - 0x3ff;
  if (e >= 0x400 || e == -0x3ff) /* x <= 0 or NaN/Inf or subnormal */
  {
    if (x <= 0.0)
    {
      /* f(x<0) is NaN, f(+/-0) is -Inf and raises DivByZero */
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
  int cancel = cr_log_fast (&h, &l, e, v);
  /* When cancel=0, i.e., x outside [sqrt(2)/2, sqrt(2)], the relative error
     is bounded by 2^-73.528 and |l/h| < 2^-23.9.
     When cancel=1, i.e., x in [sqrt(2)/2, sqrt(2)], the relative error is
     bounded by 2^-67.052, and |l/h| < 2^-52.
     In both cases, the relative error is bounded by
     err*(1 + 2^-23.9)*|h|, where err = 2^-73.528 for cancel=0,
     and err = 2^-67.052 for cancel = 1:
     2^-73.528*(1 + 2^-23.9) < 0x1.64p-74
     2^-67.052*(1 + 2^-23.9) < 0x1.eep-68 */

  static double Err[2] = {0x1.64p-74, 0x1.eep-68};

  double err = Err[cancel]; /* maximal relative error from cr_log_fast */

  double left =  h + __builtin_fma (-h, err, l);
  double right = h + __builtin_fma (h,  err, l);
  if (left == right)
    return left;
  /* the probability of failure of the fast path is about 2^-11.5 */
  return cr_log_accurate (x);
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

  // E·log(2)
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

// Extract both the significand and exponent of a double
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
