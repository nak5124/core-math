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

// #define TRACE 0x1.c19bdd1656c31p+0
// #define TRACEM 0x1.c19bdd1656c31p+0

typedef union { double f; uint64_t u; } d64u64;

/* Add a + b exactly, such that *hi + *lo = a + b.
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
     in "Tight interval inclusions with compensated algorithms"
     by Stef Graillat and Fabienne Jézéquel,
     IEEE Transactions on Computers, 2020. Proposition 3.2 says that
     the difference between a+b and hi+lo is bounded by 4u^2|a+b|
     and also by 4u^2|hi|. Here u=2^-53, thus we get:
     |(a+b)-(hi+lo)| <= 2^-104 min(|a+b|,|hi|) */
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

/* For 362 <= i <= 724, (h,l) = _LOG_INV[i-362] is a double-double nearest
   approximation of -log(_INVERSE[i-362]) */
static const double _LOG_INV[363][2] = {
    {-0x1.615ddb4bec13cp-2, -0x1.f2a42f012c6e2p-56},
    {-0x1.5e87b20c2954ap-2, 0x1.738446382fc51p-59},
    {-0x1.5baf846aa1b19p-2, -0x1.17078c05bf27dp-56},
    {-0x1.58d54f86e02f2p-2, -0x1.f30a795214b66p-57},
    {-0x1.55f9107a43ee2p-2, 0x1.91de37d2989eep-56},
    {-0x1.531ac457ee77ep-2, -0x1.dbec98a8071bdp-59},
    {-0x1.503a682cb1cb3p-2, 0x1.91e2df36b99dep-58},
    {-0x1.4ec973260026ap-2, 0x1.42a87d977dc5ep-56},
    {-0x1.4be5f957778a1p-2, 0x1.259b35b04813dp-57},
    {-0x1.49006804009d1p-2, 0x1.9ffc341f177dcp-57},
    {-0x1.4618bc21c5ec2p-2, -0x1.f42decdeccf1dp-56},
    {-0x1.432ef2a04e814p-2, 0x1.29931715ac903p-56},
    {-0x1.404308686a7e4p-2, 0x1.0bcfb6082ce6dp-56},
    {-0x1.3d54fa5c1f71p-2, 0x1.e3265c6a1c98dp-56},
    {-0x1.3a64c556945eap-2, 0x1.c68651945f97cp-57},
    {-0x1.3772662bfd85bp-2, 0x1.b5629d8117de7p-59},
    {-0x1.35f865c93293ep-2, -0x1.0d8af2d5b0557p-59},
    {-0x1.3302c16586588p-2, -0x1.7dc2a3e08aa16p-56},
    {-0x1.300aead06350cp-2, 0x1.52e91406a8a04p-57},
    {-0x1.2d10dec508583p-2, -0x1.87dc220d4d94cp-58},
    {-0x1.2a1499f762bc9p-2, -0x1.02b831548f4a6p-58},
    {-0x1.2895a13de86a3p-2, -0x1.7ad24c13f040ep-56},
    {-0x1.2596010df763ap-2, 0x1.0f76c57075e9ep-58},
    {-0x1.22941fbcf7966p-2, 0x1.76f5eb09628afp-56},
    {-0x1.1f8ff9e48a2f3p-2, 0x1.c9fdf9a0c4b07p-56},
    {-0x1.1c898c16999fbp-2, 0x1.0e5c62aff1c44p-60},
    {-0x1.1b05791f07b49p-2, 0x1.466dc55e2d052p-56},
    {-0x1.17fb98e15095dp-2, -0x1.7458b5d97ba9dp-56},
    {-0x1.14ef67f88685ap-2, -0x1.66880da1b13e4p-58},
    {-0x1.136870293a8bp-2, -0x1.7b66298edd24ap-56},
    {-0x1.1058bf9ae4ad5p-2, -0x1.89fa0ab4cb31dp-58},
    {-0x1.0d46b579ab74bp-2, -0x1.03ec81c3cbd92p-57},
    {-0x1.0a324e27390e3p-2, -0x1.7dcfde8061c03p-56},
    {-0x1.08a73667c57afp-2, -0x1.d40c5a328e6c6p-60},
    {-0x1.058f3c703ebc6p-2, 0x1.9af348dab5c3cp-58},
    {-0x1.0402594b4d041p-2, 0x1.28ec217a5022dp-57},
    {-0x1.00e6c45ad501dp-2, 0x1.cb9568ff6feadp-57},
    {-0x1.fb9186d5e3e2bp-3, 0x1.caaae64f21acbp-57},
    {-0x1.f871b28955045p-3, -0x1.4ad6c8812d31ap-63},
    {-0x1.f22e5e72f105dp-3, -0x1.61d7d037c1899p-57},
    {-0x1.ebe61f4dd7b0bp-3, 0x1.4cc3f7293285cp-59},
    {-0x1.e8c0252aa5a6p-3, 0x1.6e03a39bfc89bp-59},
    {-0x1.e27076e2af2e6p-3, 0x1.61578001e0162p-59},
    {-0x1.dc1bca0abec7dp-3, -0x1.834c51998b6fcp-57},
    {-0x1.d8ef91af31d5ep-3, 0x1.7f0d931e0e2cap-60},
    {-0x1.d293581b6b3e7p-3, 0x1.c04a2aa97ac8ep-58},
    {-0x1.cf6354e09c5dcp-3, -0x1.239a07d55b695p-57},
    {-0x1.c8ff7c79a9a22p-3, 0x1.4f689f8434012p-57},
    {-0x1.c5cba543ae425p-3, 0x1.62134bab038d8p-57},
    {-0x1.bf601bb0e44e2p-3, 0x1.56b83c874aaf4p-57},
    {-0x1.bc286742d8cd6p-3, -0x1.4fce744870f55p-58},
    {-0x1.b5b519e8fb5a4p-3, -0x1.ba27fdc19e1ap-57},
    {-0x1.af3c94e80bff3p-3, 0x1.398cff3641985p-58},
    {-0x1.abfe5ae46124cp-3, 0x1.ea72be27390fp-57},
    {-0x1.a57df28244dcdp-3, 0x1.b9af132a24e39p-59},
    {-0x1.a23bc1fe2b563p-3, -0x1.93711b07a998cp-59},
    {-0x1.9bb362e7dfb83p-3, -0x1.575e31f003e0cp-57},
    {-0x1.986d3228180cap-3, 0x1.2a6c8af000189p-58},
    {-0x1.91dcc8c340bdep-3, -0x1.aaf77bfd17182p-58},
    {-0x1.8e928de886d41p-3, 0x1.569d851a5677p-57},
    {-0x1.87fa06520c911p-3, 0x1.bf7fdbfa08d9ap-57},
    {-0x1.84abb75865139p-3, -0x1.5482c750b9638p-58},
    {-0x1.815c0a14357ebp-3, 0x1.4be48073a0564p-58},
    {-0x1.7ab890210d909p-3, -0x1.be36b2d6a0608p-59},
    {-0x1.7764c128f2127p-3, -0x1.240d1e78f44cep-57},
    {-0x1.70b8f97a1aa75p-3, -0x1.93daa56fbafd6p-58},
    {-0x1.6d60fe719d21dp-3, 0x1.caae268ecd179p-57},
    {-0x1.66acd4272ad51p-3, 0x1.0900e4e1ea8b2p-58},
    {-0x1.6350a28aaa758p-3, 0x1.3f547e9c51633p-57},
    {-0x1.5ff3070a793d4p-3, 0x1.bc60efafc6f6ep-58},
    {-0x1.59338d9982086p-3, 0x1.65d22aa8ad7cfp-58},
    {-0x1.55d1ad4232d6fp-3, 0x1.ac8966e060839p-58},
    {-0x1.4f099f4a230b2p-3, -0x1.a0a02a1b24794p-61},
    {-0x1.4ba36f39a55e5p-3, -0x1.68981bcc36756p-57},
    {-0x1.483bccce6e3ddp-3, -0x1.29391fb1b4b22p-57},
    {-0x1.41682bf727bcp-3, 0x1.1c207e127261bp-59},
    {-0x1.3dfc2b0ecc62ap-3, 0x1.ab3a8e7d81017p-58},
    {-0x1.371fc201e8f74p-3, -0x1.de6cb62af18ap-58},
    {-0x1.33af575770e4fp-3, -0x1.9945fce5491eap-57},
    {-0x1.303d718e47fd3p-3, 0x1.6b9c7d96091fap-63},
    {-0x1.29552f81ff523p-3, -0x1.301771c407dbfp-57},
    {-0x1.25ded0abc6ad2p-3, 0x1.eac3a26edd19cp-58},
    {-0x1.2266f190a5acbp-3, -0x1.f547bf1809e88p-57},
    {-0x1.1b72ad52f67ap-3, -0x1.483023472cd74p-58},
    {-0x1.17f6458fca611p-3, 0x1.4bdb0dc8fdffap-63},
    {-0x1.14785846742acp-3, -0x1.a28813e3a7f07p-57},
    {-0x1.0d77e7cd08e59p-3, -0x1.9a5dc5e9030acp-57},
    {-0x1.09f561ee719c3p-3, -0x1.f51d505cb0b76p-58},
    {-0x1.0671512ca596ep-3, -0x1.50c647eb86499p-58},
    {-0x1.02ebb42bf3d4bp-3, 0x1.f4b9c01cb92c6p-59},
    {-0x1.f7b79fec37ddfp-4, 0x1.87e897ed01783p-59},
    {-0x1.f0a30c01162a6p-4, -0x1.85f325c5bbacdp-58},
    {-0x1.e98b549671467p-4, -0x1.d227143a5a998p-58},
    {-0x1.e27076e2af2e6p-4, 0x1.61578001e0162p-60},
    {-0x1.d4313d66cb35dp-4, -0x1.790dd951d90fap-58},
    {-0x1.cd0cdbf8c13e1p-4, -0x1.36d4375d0c271p-58},
    {-0x1.c5e548f5bc743p-4, -0x1.5d617ef8161b1p-60},
    {-0x1.b78c82bb0eda1p-4, -0x1.0878cf0327e21p-61},
    {-0x1.b05b49bee43fep-4, -0x1.160c7c252f298p-58},
    {-0x1.a926d3a4ad563p-4, -0x1.942f48aa70ea9p-58},
    {-0x1.a1ef1d8061cd4p-4, -0x1.76df97bcb177fp-60},
    {-0x1.9ab42462033adp-4, 0x1.2099e1c184e8ep-59},
    {-0x1.8c345d6319b21p-4, 0x1.4a697ab3424a9p-61},
    {-0x1.84ef898e8282ap-4, -0x1.96dcb441b9227p-59},
    {-0x1.7da766d7b12cdp-4, 0x1.eeedfcdd94131p-58},
    {-0x1.765bf23a6be13p-4, -0x1.0ff28ef6a592fp-58},
    {-0x1.6f0d28ae56b4cp-4, 0x1.906d99184b992p-58},
    {-0x1.60658a93750c4p-4, 0x1.388458ec21b6ap-58},
    {-0x1.590cafdf01c28p-4, -0x1.3d5c8aaea76d2p-58},
    {-0x1.51b073f06183fp-4, -0x1.a49e39a1a8be4p-58},
    {-0x1.4a50d3aa1b04p-4, -0x1.ecf768c1dd57bp-61},
    {-0x1.42edcbea646fp-4, -0x1.ddd4f935996c9p-59},
    {-0x1.341d7961bd1d1p-4, 0x1.b599f227becbbp-58},
    {-0x1.2cb0283f5de1fp-4, 0x1.d359a8fde8adep-60},
    {-0x1.253f62f0a1417p-4, 0x1.c125963fc4cfdp-62},
    {-0x1.1dcb263db1944p-4, -0x1.3d7a7a2605718p-58},
    {-0x1.16536eea37ae1p-4, 0x1.79da3e8c22cdap-60},
    {-0x1.0ed839b5526fep-4, -0x1.7256ea8988a68p-61},
    {-0x1.075983598e471p-4, -0x1.80da5333c45b8p-59},
    {-0x1.f0a30c01162a6p-5, -0x1.85f325c5bbacdp-59},
    {-0x1.e19070c276016p-5, 0x1.19918a7a17dc1p-59},
    {-0x1.d276b8adb0b52p-5, -0x1.1e3c53257fd47p-61},
    {-0x1.c355dd0921f2dp-5, 0x1.9b2a03e3be3a7p-60},
    {-0x1.b42dd711971bfp-5, 0x1.eb9759c130499p-60},
    {-0x1.a4fe9ffa3d235p-5, 0x1.28100a49366b4p-62},
    {-0x1.95c830ec8e3ebp-5, -0x1.f5a0e80520bf2p-59},
    {-0x1.868a83083f6cfp-5, 0x1.d09a5634943dbp-61},
    {-0x1.77458f632dcfcp-5, -0x1.18d3ca87b9296p-59},
    {-0x1.58a5bafc8e4d5p-5, 0x1.ce55c2b4e2b72p-59},
    {-0x1.494acc34d911cp-5, -0x1.e295bf491ccc5p-59},
    {-0x1.39e87b9febd6p-5, 0x1.5bfa937f551bbp-59},
    {-0x1.2a7ec2214e873p-5, -0x1.8856e9c01e6ddp-61},
    {-0x1.1b0d98923d98p-5, 0x1.e9ae889bac481p-60},
    {-0x1.0b94f7c196176p-5, -0x1.da43f761f4dc4p-59},
    {-0x1.f829b0e7833p-6, -0x1.33e3f04f1ef23p-60},
    {-0x1.d91a66c543cc4p-6, 0x1.d34e608cbdaabp-62},
    {-0x1.b9fc027af9198p-6, 0x1.0ae69229dc868p-64},
    {-0x1.9ace7551cc514p-6, -0x1.3409c1df8167fp-60},
    {-0x1.7b91b07d5b11bp-6, 0x1.5b602ace3a51p-60},
    {-0x1.5c45a51b8d389p-6, 0x1.b10b6c3ec21b4p-60},
    {-0x1.3cea44346a575p-6, 0x1.0cb5a902b3a1cp-62},
    {-0x1.1d7f7eb9eebe7p-6, 0x1.d41fe63d2dbf9p-61},
    {-0x1.fc0a8b0fc03e4p-7, 0x1.83092c59642a1p-62},
    {-0x1.bcf712c74384cp-7, 0x1.f6842688f499ap-62},
    {-0x1.7dc475f810a77p-7, 0x1.16d7687d3df21p-62},
    {-0x1.3e7295d25a7d9p-7, 0x1.ff29a11443a06p-65},
    {-0x1.fe02a6b106789p-8, 0x1.e44b7e3711ebfp-67},
    {-0x1.7ee11ebd82e94p-8, 0x1.61e96e2fc5d9p-62},
    {-0x1.ff00aa2b10bcp-9, -0x1.2821ad5a6d353p-63},
    {-0x1.ff802a9ab10e6p-10, -0x1.e29e3a153e3b2p-64},
    {0x1.0020055655889p-10, 0x1.9abe6bf0fa436p-65},
    {0x1.80904828985cp-9, 0x1.a5a9c30313fb6p-63},
    {0x1.40c8a747878e2p-8, -0x1.c41f1fffceda2p-63},
    {0x1.c189cbb0e27fbp-8, -0x1.44962158133afp-62},
    {0x1.2145e939ef1e9p-7, -0x1.b8e762c009922p-61},
    {0x1.61e77e8b53fc6p-7, 0x1.1507150cb4034p-64},
    {0x1.a2a9c6c170462p-7, 0x1.d94965c6f7205p-63},
    {0x1.e38ce3033310cp-7, -0x1.147b45033e1b2p-61},
    {0x1.12487a5507f7p-6, -0x1.a012b4c6d7e53p-60},
    {0x1.32db0ea132e22p-6, -0x1.a767bb50221f5p-60},
    {0x1.537e3f45f3565p-6, -0x1.3169406d66a6ap-60},
    {0x1.63d6178690bd6p-6, -0x1.8ed4d357c9c97p-64},
    {0x1.8492528c8cabfp-6, -0x1.d192d0619fa67p-60},
    {0x1.a55f548c5c43fp-6, 0x1.ec1a5f86d41f9p-62},
    {0x1.c63d2ec14aaf2p-6, -0x1.ce030a686bd86p-60},
    {0x1.e72bf2813ce51p-6, 0x1.75b44595cab18p-60},
    {0x1.0415d89e74444p-5, 0x1.c05cf1d753622p-59},
    {0x1.149e3e4005a8dp-5, -0x1.53482d1f9d7d7p-61},
    {0x1.252f32f8d183fp-5, -0x1.947f792615916p-59},
    {0x1.35c8bfaa1306bp-5, -0x1.50830a65543a4p-63},
    {0x1.3e18c1ca0ae92p-5, -0x1.2c091c871b7bap-60},
    {0x1.4ebf43349e26fp-5, -0x1.43dcef9dcdaeap-59},
    {0x1.5f6e73078efb8p-5, -0x1.6a0049252e1b9p-59},
    {0x1.70265a550e777p-5, 0x1.07701518c666p-60},
    {0x1.80e7023d8ccc4p-5, 0x1.f486ba05d8df6p-59},
    {0x1.91b073efd7314p-5, 0x1.d60449ab527bfp-61},
    {0x1.9a187b573de7cp-5, -0x1.727626c86b3abp-59},
    {0x1.aaef2d0fb10fcp-5, 0x1.a353bb42e0addp-61},
    {0x1.bbcebfc68f42p-5, 0x1.e5cf3a0f56f72p-60},
    {0x1.ccb73cdddb2ccp-5, -0x1.e48fb0500efd4p-59},
    {0x1.dda8adc67ee4ep-5, 0x1.4e6c986f44c55p-59},
    {0x1.e624c4a0b5e1bp-5, 0x1.1c5cc4cbb986dp-59},
    {0x1.f723b517fc523p-5, -0x1.86eab086959p-59},
    {0x1.0415d89e74444p-4, 0x1.c05cf1d753622p-58},
    {0x1.0c9e615ac4e17p-4, -0x1.5fda2c9a276b4p-59},
    {0x1.10e45b3cae831p-4, -0x1.a4a128d192686p-58},
    {0x1.1973bd1465567p-4, -0x1.7558367a6acf6p-59},
    {0x1.2207b5c78549ep-4, -0x1.cc0fbce104eaap-58},
    {0x1.2aa04a44717a5p-4, -0x1.d15d38d2fa3f7p-58},
    {0x1.2eee507b40301p-4, 0x1.edd77c85fad4p-63},
    {0x1.378dd7f749714p-4, 0x1.128f1faca0abep-60},
    {0x1.403207b414b7fp-4, -0x1.ed55faa100316p-59},
    {0x1.4485e03dbdfadp-4, 0x1.1ba349aadbc6ep-58},
    {0x1.4d3115d207eacp-4, 0x1.769f42c7842ccp-58},
    {0x1.55e10050e0384p-4, -0x1.45f9d61c68c1bp-58},
    {0x1.5e95a4d9791cbp-4, 0x1.f38745c5c450ap-58},
    {0x1.62f1be7d77743p-4, 0x1.aa04a70ec8925p-59},
    {0x1.6bad83c1883b6p-4, -0x1.867eedb2a4fafp-60},
    {0x1.746e100226ed9p-4, 0x1.748f0ef16ce8dp-59},
    {0x1.78d02263d82d3p-4, 0x1.abca5b4fdb88p-58},
    {0x1.8197e2f40e3fp-4, 0x1.b9f2dffbeed43p-60},
    {0x1.8a6477a91dc29p-4, -0x1.fa83214904842p-59},
    {0x1.8ecc933aeb6e8p-4, 0x1.8c833010ab15cp-58},
    {0x1.97a07024cbe74p-4, 0x1.9f466fb618c1dp-59},
    {0x1.a0792e9277cacp-4, 0x1.8c9b2957205c6p-58},
    {0x1.a4e7640b1bc38p-4, -0x1.5b5ca203e4259p-58},
    {0x1.adc77ee5aea8cp-4, 0x1.37d8f39bee659p-58},
    {0x1.b23965a52ffp-4, 0x1.3622bd91f0d86p-58},
    {0x1.bb20e936d6974p-4, 0x1.5f2ae991c8844p-58},
    {0x1.c40d6425a5cb1p-4, 0x1.21d1930dc8acdp-60},
    {0x1.c885801bc4b23p-4, 0x1.a38cb559a6706p-58},
    {0x1.d179788219364p-4, 0x1.9daf7df76ad2ap-59},
    {0x1.d5f55659210e2p-4, 0x1.ce60c2a34a8fbp-59},
    {0x1.def0d8d466db9p-4, -0x1.4104ba69d1fd8p-58},
    {0x1.e7f1691a32d3ep-4, 0x1.d321c420330f3p-59},
    {0x1.ec739830a112p-4, -0x1.a2bf991780d3fp-59},
    {0x1.f57bc7d9005dbp-4, -0x1.9361574fb24e2p-58},
    {0x1.fa01c9db57ce2p-4, -0x1.a8fd2453980aap-58},
    {0x1.0188d2ecf614p-3, -0x1.a6b8c01880701p-57},
    {0x1.03cdc0a51ec0dp-3, 0x1.39e2d3f8b7d1p-57},
    {0x1.08598b59e3a07p-3, -0x1.dd7009902bf32p-57},
    {0x1.0aa06912675d5p-3, 0x1.cc51f9bdae72dp-57},
    {0x1.0f301717cf0fbp-3, 0x1.2ef945e4ed0c2p-58},
    {0x1.13c2605c398c3p-3, -0x1.fdd94f6508b88p-57},
    {0x1.160c8024b27b1p-3, -0x1.2d56ff61c2bfbp-57},
    {0x1.1aa2b7e23f72ap-3, -0x1.c6ef1d9b2ef7ep-59},
    {0x1.1ceed09853752p-3, -0x1.c39192af8eb16p-60},
    {0x1.2188fd9807263p-3, -0x1.e7f50c701268fp-60},
    {0x1.23d712a49c202p-3, -0x1.6e38161051d69p-57},
    {0x1.28753bc11aba5p-3, -0x1.6394d9fa33311p-57},
    {0x1.2ac55095f5c59p-3, 0x1.cc979c9bac1cdp-57},
    {0x1.2f677cbbc0a96p-3, -0x1.9fbd3e17e5527p-57},
    {0x1.31b994d3a4f85p-3, -0x1.c4716bdfc0cc9p-58},
    {0x1.365fcb0159016p-3, 0x1.7d411a5b944adp-58},
    {0x1.38b3e9e027479p-3, -0x1.e0a8b82ff1eb2p-57},
    {0x1.3d5e3126bc27fp-3, 0x1.97c284b6258aap-57},
    {0x1.3fb45a59928ccp-3, -0x1.d87e6a354d056p-57},
    {0x1.420b32740fdd4p-3, -0x1.41642d0442295p-57},
    {0x1.46baf0f9f5db7p-3, -0x1.cd0790841a3eep-58},
    {0x1.4913d8333b561p-3, -0x1.0d5604930f135p-58},
    {0x1.4dc7b897bc1c8p-3, -0x1.927d47803c5f4p-57},
    {0x1.5022b292f6a45p-3, 0x1.7fcda896de0e9p-62},
    {0x1.54dabc26105d2p-3, -0x1.011a372f27d11p-57},
    {0x1.5737cc9018cddp-3, 0x1.4f4d710fec38ep-57},
    {0x1.5bf406b543db2p-3, -0x1.1f5b44c0df7e7p-61},
    {0x1.5e533144c1719p-3, -0x1.7e76521d4f8b4p-61},
    {0x1.60b3100b09476p-3, -0x1.5b2623e05016bp-58},
    {0x1.6574ebe8c133ap-3, -0x1.d34f0f4621bedp-60},
    {0x1.67d6e9d785771p-3, -0x1.10614e0da5fb8p-57},
    {0x1.6c9d07d203fc7p-3, 0x1.80a04c9a46c61p-59},
    {0x1.6f0128b756abcp-3, -0x1.8de59c21e166cp-57},
    {0x1.716600c914054p-3, 0x1.b157cec383873p-57},
    {0x1.7631d82935a86p-3, 0x1.047074183dfcbp-58},
    {0x1.7898d85444c73p-3, 0x1.ef8f6ebcfb201p-58},
    {0x1.7d6903caf5adp-3, -0x1.ac5f0c075b847p-59},
    {0x1.7fd22ff599d4fp-3, -0x1.f5fa2bdbe8c07p-57},
    {0x1.823c16551a3c2p-3, -0x1.1232ce70be781p-57},
    {0x1.871213750e994p-3, 0x1.d685f35eea2ap-57},
    {0x1.897e2b17b19a5p-3, 0x1.918fe682c4831p-58},
    {0x1.8beafeb38fe8cp-3, 0x1.55aa8b6997a4p-58},
    {0x1.90c6db9fcbcd9p-3, 0x1.054473941ad99p-57},
    {0x1.9335e5d594989p-3, -0x1.478a85704ccb7p-57},
    {0x1.95a5adcf7017fp-3, 0x1.142c507fb7a3dp-58},
    {0x1.9a8778debaa38p-3, 0x1.f47dfd871f87fp-57},
    {0x1.9cf97cdce0ec3p-3, 0x1.788620a71b828p-59},
    {0x1.9f6c407089664p-3, 0x1.35a19605e67efp-59},
    {0x1.a454082e6ab05p-3, 0x1.df207dc5c34c6p-58},
    {0x1.a6c90d44b704ep-3, 0x1.203f213ce9578p-58},
    {0x1.a93ed3c8ad9e3p-3, 0x1.bcafa9de97203p-57},
    {0x1.ae2ca6f672bd4p-3, 0x1.ab5ca9eaa088ap-57},
    {0x1.b0a4b48fc1b46p-3, 0x1.a5478cce26344p-58},
    {0x1.b31d8575bce3dp-3, -0x1.6353ab386a94dp-57},
    {0x1.b811730b823d2p-3, 0x1.a0ee735d9f0ecp-60},
    {0x1.ba8c90ae4ad19p-3, 0x1.cfe88865b42bdp-57},
    {0x1.bd087383bd8adp-3, 0x1.dd355f6a516d7p-60},
    {0x1.c2028ab17f9b4p-3, 0x1.f11aa3853a5f1p-57},
    {0x1.c480c0005ccd1p-3, 0x1.29abc89ceca68p-57},
    {0x1.c6ffbc6f00f71p-3, -0x1.8e58b2c57a4a5p-57},
    {0x1.c97f8079d44ecp-3, 0x1.61a8c6e6c4ee7p-57},
    {0x1.ce816157f1988p-3, -0x1.5744132a297bp-58},
    {0x1.d1037f2655e7bp-3, 0x1.60629242471a2p-57},
    {0x1.d38666871f465p-3, -0x1.42ce24caf02f9p-58},
    {0x1.d88e93fb2f45p-3, 0x1.affb96815e081p-57},
    {0x1.db13db0d4894p-3, 0x1.aa11d49f96cb9p-58},
    {0x1.dd99edaf6d7e9p-3, -0x1.8cd38eadd64c6p-57},
    {0x1.e020cc6235ab5p-3, 0x1.fea48dd7b81d1p-58},
    {0x1.e530effe71012p-3, 0x1.2276041f43042p-59},
    {0x1.e7ba35eb77e2ap-3, 0x1.11dc86c9b7564p-59},
    {0x1.ea4449f04aaf5p-3, -0x1.d33919ab94074p-57},
    {0x1.eccf2c8fe920ap-3, -0x1.8e8f9d590196p-58},
    {0x1.ef5ade4dcffe6p-3, -0x1.08ab2ddc708ap-58},
    {0x1.f474b134df229p-3, -0x1.27c77ded76aadp-58},
    {0x1.f702d36777dfp-3, -0x1.9516673e2299cp-58},
    {0x1.f991c6cb3b379p-3, 0x1.f665066f980a2p-57},
    {0x1.fc218be620a5ep-3, -0x1.6e438c258187fp-58},
    {0x1.00a1c6adda473p-2, 0x1.8d688b9e17a8ap-56},
    {0x1.01eae5626c691p-2, -0x1.18290bd2932e2p-59},
    {0x1.03346e0106062p-2, 0x1.ff8a966395c73p-56},
    {0x1.047e60cde83b8p-2, -0x1.0779634061cbcp-56},
    {0x1.05c8be0d9635ap-2, 0x1.e38ef996b0c96p-58},
    {0x1.085eb8f8ae797p-2, 0x1.513f45fe7a977p-56},
    {0x1.09aa572e6c6d4p-2, 0x1.43c2e68684d53p-57},
    {0x1.0af660eb9e279p-2, -0x1.e056b93fb20cfp-60},
    {0x1.0c42d676162e3p-2, 0x1.162c79d5d11eep-58},
    {0x1.0d8fb813eb1efp-2, -0x1.cdde2b0172bd5p-56},
    {0x1.102ac0a35cc1cp-2, 0x1.088080a5e68b4p-59},
    {0x1.1178e8227e47cp-2, -0x1.0e63a5f01c691p-57},
    {0x1.12c77cd00713bp-2, 0x1.4a4508fbcba26p-57},
    {0x1.14167ef367783p-2, 0x1.e0936abd4fa6ep-62},
    {0x1.1565eed455fc3p-2, 0x1.4a4092a8bb5ep-58},
    {0x1.16b5ccbacfb73p-2, 0x1.66fbd28b40935p-56},
    {0x1.1956d3b9bc2fap-2, 0x1.7b9d68d50a15dp-56},
    {0x1.1aa7fd638d33fp-2, -0x1.f569e908600b2p-57},
    {0x1.1bf99635a6b95p-2, -0x1.12aeb84249223p-57},
    {0x1.1d4b9e796c245p-2, 0x1.333e2172b6715p-56},
    {0x1.1e9e1678899f4p-2, 0x1.512c3749a1e4ep-56},
    {0x1.1ff0fe7cf47a7p-2, 0x1.5b513ff0c145p-56},
    {0x1.214456d0eb8d4p-2, 0x1.f7ae91aeba60ap-57},
    {0x1.23ec5991eba49p-2, 0x1.bb75d1addf87p-60},
    {0x1.25410494e56c7p-2, 0x1.7ac0ef77f252ap-56},
    {0x1.269621134db92p-2, 0x1.e0efadd9db02bp-56},
    {0x1.27ebaf58d8c9dp-2, -0x1.8800b4bda6c97p-57},
    {0x1.2941afb186b7cp-2, -0x1.856e61c51574p-57},
    {0x1.2a982269a3dbfp-2, 0x1.38d546bd18905p-56},
    {0x1.2bef07cdc9354p-2, -0x1.82dad7fd86088p-56},
    {0x1.2d46602adcceep-2, 0x1.7911955f3520fp-56},
    {0x1.2ff66b04ea9d4p-2, 0x1.2dabe191d1c94p-56},
    {0x1.314f1e1d35ce4p-2, -0x1.3d69909e5c3dcp-56},
    {0x1.32a84565120a8p-2, -0x1.b5fb19428a75ep-57},
    {0x1.3401e12aecba1p-2, -0x1.cd55b8a4746cp-58},
    {0x1.355bf1bd82c8bp-2, -0x1.9b8964f0e80dcp-57},
    {0x1.36b6776be1117p-2, -0x1.324f0e883858ep-58},
    {0x1.3811728564cb2p-2, -0x1.e493a0702b236p-57},
    {0x1.396ce359bbf54p-2, -0x1.ce2b31b31e8bp-58},
    {0x1.3ac8ca38e5c5fp-2, -0x1.f7de015f253eep-56},
    {0x1.3c25277333184p-2, -0x1.2ad27e50a8ec6p-56},
    {0x1.3d81fb5946dbap-2, 0x1.c1eab1642e36dp-56},
    {0x1.3edf463c1683ep-2, 0x1.83d680d3c1084p-56},
    {0x1.419b423d5e8c7p-2, 0x1.0dbb243827392p-57},
    {0x1.42f9f3ff62642p-2, -0x1.bbf082ccabbaep-56},
    {0x1.44591e0539f49p-2, -0x1.2b125247b0fa5p-56},
    {0x1.45b8c0a17df13p-2, 0x1.dbe305eaf5a2p-56},
    {0x1.4718dc271c41bp-2, 0x1.8fb4c14c56eefp-60},
    {0x1.487970e95877p-2, 0x1.b8465cf25f4c6p-56},
    {0x1.49da7f3bcc41fp-2, -0x1.9964a168ccacap-57},
    {0x1.4b3c077267e9ap-2, 0x1.2e5fbeb518508p-56},
    {0x1.4c9e09e172c3cp-2, -0x1.123615b147a5dp-58},
    {0x1.4e0086dd8bacap-2, 0x1.6d5e1bb877c2ep-56},
    {0x1.4f637ebba981p-2, -0x1.58cb3124b9245p-56},
    {0x1.50c6f1d11b97cp-2, 0x1.94817d83d3edp-56},
    {0x1.522ae0738a3d8p-2, -0x1.8f7e9b38a6979p-57},
    {0x1.538f4af8f72fep-2, -0x1.fe8dd55c19315p-56},
    {0x1.54f431b7be1a9p-2, -0x1.aacfdbbdab914p-56},
    {0x1.565995069514cp-2, 0x1.7f4aeb71dce6p-56},
    {0x1.57bf753c8d1fbp-2, -0x1.0908d15f88b63p-57},
    {0x1.5925d2b112a59p-2, 0x1.4b0b52198dbd9p-61},
    {0x1.5a8cadbbedfa1p-2, -0x1.e6c2bdfb3e037p-58},
    {0x1.5bf406b543db2p-2, -0x1.1f5b44c0df7e7p-60},
    {0x1.5d5bddf595f3p-2, -0x1.6541148cbb8a2p-56},
    {0x1.5ec433d5c35aep-2, -0x1.cbdbac5d0228ep-57},
    {0x1.602d08af091ecp-2, -0x1.6e8920c09b73fp-58},
    {0x1.61965cdb02c1fp-2, -0x1.5ac078911cb74p-58},
    {0x1.630030b3aac49p-2, 0x1.dc18ce51fff99p-57},
};

/* The following is a degree-6 polynomial generated by Sollya over
   [-0.00202941894531250,0.00212097167968735],
   with absolute error < 2^-70.278.
   The polynomial is P[0]*x + P[1]*x^2 + ... + P[5]*x^6. */
static const double P[] = {0x1p0,                 /* degree 1 */
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

/* Given 1 <= x < 2, where x = v.f:
 * if x < sqrt(2): put in h+l a double-double approximation of log(x)
   and leaves e unchanged
 * if x > sqrt(2): put in h+l a double-double approximation of log(x/2),
   and increases e by 1
*/
static void
cr_log_fast (double *h, double *l, int *e, d64u64 v)
{
  uint64_t m = 0x10000000000000 + (v.u & 0xfffffffffffff);
  /* x = m/2^52 */
  //  if (v.f == TRACEM) printf ("m=%lu\n", m);
  /* if x > sqrt(2), we divide it by 2 to avoid cancellation */
  int c = m >= 0x16a09e667f3bcd;
  *e += c;
  // if (v.f == TRACEM) printf ("c=%d\n", c);
  static const double cy[] = {1.0, 0.5};
  static const uint64_t cm[] = {43, 44};

  int i = m >> cm[c];
  // if (v.f == TRACEM) printf ("i=%d\n", i);
  double y = v.f * cy[c];
  // if (v.f == TRACEM) printf ("y=%la\n", y);
#define OFFSET 362
  double r = (_INVERSE - OFFSET)[i];
  double l1 = (_LOG_INV - OFFSET)[i][0];
  double l2 = (_LOG_INV - OFFSET)[i][1];
  // if (v.f == TRACEM) printf ("r=%la\n", r);
  double z = __builtin_fma (r, y, -1.0); /* exact */
  // if (v.f == TRACEM) printf ("z=%la\n", z);
  /* evaluate P(z) */
  double ph, pl, z2 = z * z;
  double p56 = __builtin_fma (P[6], z, P[5]);
  double p34 = __builtin_fma (P[4], z, P[3]);
  ph = __builtin_fma (p56, z2, p34);
  double p12 = __builtin_fma (P[2], z, P[1]);
  ph = __builtin_fma (ph, z2, p12);
  // if (v.f == TRACEM) printf ("ph=%la z=%la P[0]=%la\n", ph, z, P[0]);
  ph *= z2;
  // if (v.f == TRACEM) printf ("ph=%la\n", ph);
  /* add z since P[0]=1 */
  fast_two_sum (&ph, &pl, z, ph);
  // if (v.f == TRACEM) printf ("log(1+z): h=%la l=%la\n", ph, pl);
  /* add l1 + l2 */
  fast_two_sum (h, l, l1, ph);
  *l += pl + l2;
  // if (v.f == TRACEM) printf ("log(y): h=%la l=%la\n", *h, *l);
}

static inline void dint_fromd (dint64_t *a, double b);
static void log_2 (dint64_t *r, dint64_t *x);
static inline double dint_tod (dint64_t *a);

/* accurate path, using Tom Hubrecht's code below */
static double
cr_log_accurate (double x)
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
  // if (x == TRACE) printf ("X=(%lu,%ld,%lu,%lu)\n", X.sgn, X.ex, X.hi, X.lo);
  log_2 (&Y, &X);
  // if (x == TRACE) printf ("Y=(%lu,%ld,%lu,%lu)\n", Y.sgn, Y.ex, Y.hi, Y.lo);
  return dint_tod (&Y);
}

double
cr_log (double x)
{
  /* now x > 0 */
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
  /* normalize v in [1,2) */
  v.u = (0x3fful << 52) | (v.u & 0xfffffffffffff);
  /* now x = m*2^e with 1 <= m < 2 (m = v.f) */
  double h, l;
  // if (x == TRACE) printf ("x=%la e=%d m=%la\n", x, e, v.f);
  cr_log_fast (&h, &l, &e, v);
  /* err=0x1.21p-69 + ... fails for x=0x1.830124938cfeap-85 (rndz) */
  static double err = 0x1.22p-69 + 0x1.04p-85;
  /* 0x1.04p-85 is the maximal error for the addition of e*log(2) below */

  /* Add e*log(2) to (h,l), where -1074 <= e <= 1023, thus e has at most
     11 bits. We store log2_h on 42 bits, so that e*log2_h is exact. */
  static double log2_h = 0x1.62e42fefa38p-1, log2_l = 0x1.ef35793c7673p-45;
  /* |log(2) - (h+l)| < 2^-102.01 */
  double hh = e * log2_h; /* exact */
  double ll = __builtin_fma (e, log2_l, l);
  /* we have |l| < 2^-50 (from the analysis of cr_log_fast)
     and |e*log2_l| <= 1074*0x1.ef35793c7673p-45
     thus |ll| < 2^-33.9 and err(ll) <= ulp(2^-33.9) = 2^-86 */
  fast_two_sum (&h, &l, hh, h); /* rounding error bounded by 2^-104*|hh| < 2^-94.45 */
  l += ll; /* |l| < 2^-50 and |ll| < 2^-33.9 thus |l+ll| < 2^-33.8
              and the rounding error is less than ulp(2^-33.8) = 2^-86 */
  /* Additional rounding error:
     - e times the approximation error on log2_h+log2_l: 1074*2^-102.01 < 2^-91.94
     - error on ll < 2^-86
     - error in the fast_two_sum < 2^-94.45
     - error in l + ll < 2^-86
     Total < 2^-84.98 < 1.04e-85 */

  // if (x == TRACE) printf ("h=%la l=%la err=%la\n", h, l, err);
  double left = h + (l - err), right = h + (l + err);
  // if (x == TRACE) printf ("left=%la right=%la\n", left, right);
  if (left == right)
  {
    // if (x == TRACE) printf ("fast path succeeded\n");
    return left;
  }
  // if (x == TRACE) printf ("fast path failed\n");
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
