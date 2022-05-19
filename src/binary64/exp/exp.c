/* Correctly rounded exponential function for binary64 values.

Copyright (c) 2021-2022 Paul Zimmermann, Inria.

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

#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include <stdint.h>
#include <math.h>
#include <fenv.h>
#include <assert.h>

/* Add a + b exactly, such that *hi + *lo = a + b.
   Assumes |a| >= |b| and rounding to nearest.  */
static void
fast_two_sum (double *hi, double *lo, double a, double b)
{
  double e;

  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
  /* Now *hi + *lo = a + b exactly.  */
}

/* h + l <- a * b
   We have h + l = a + b exactly, whatever the rounding mode, when no
   underflow happens (cf Section 4.4 of the Handbook of Floating-Point
   Arithmetic, 2nd edition) */
static void
dekker (double *h, double *l, double a, double b)
{
  *h = a * b;
  *l = __builtin_fma (a, b, -*h);
}

typedef union { double x; uint64_t n; } d64u64;

/* for -127 <= i <= 127, tab_i[127+i] is a double-double approximation
   of 2^(i/128) */
static const double tab_i[255][2] = {
{ 0x1.0163da9fb3335p-1, 0x1.b61299ab8cdb7p-55 }, /* -127 */
{ 0x1.02c9a3e778061p-1, -0x1.19083535b085dp-57 }, /* -126 */
{ 0x1.04315e86e7f85p-1, -0x1.0a31c1977c96ep-55 }, /* -125 */
{ 0x1.059b0d3158574p-1, 0x1.d73e2a475b465p-56 }, /* -124 */
{ 0x1.0706b29ddf6dep-1, -0x1.c91dfe2b13c27p-56 }, /* -123 */
{ 0x1.0874518759bc8p-1, 0x1.186be4bb284ffp-58 }, /* -122 */
{ 0x1.09e3ecac6f383p-1, 0x1.1487818316136p-55 }, /* -121 */
{ 0x1.0b5586cf9890fp-1, 0x1.8a62e4adc610bp-55 }, /* -120 */
{ 0x1.0cc922b7247f7p-1, 0x1.01edc16e24f71p-55 }, /* -119 */
{ 0x1.0e3ec32d3d1a2p-1, 0x1.03a1727c57b53p-60 }, /* -118 */
{ 0x1.0fb66affed31bp-1, -0x1.b9bedc44ebd7bp-58 }, /* -117 */
{ 0x1.11301d0125b51p-1, -0x1.6c51039449b3ap-55 }, /* -116 */
{ 0x1.12abdc06c31ccp-1, -0x1.1b514b36ca5c7p-59 }, /* -115 */
{ 0x1.1429aaea92dep-1, -0x1.32fbf9af1369ep-55 }, /* -114 */
{ 0x1.15a98c8a58e51p-1, 0x1.2406ab9eeab0ap-56 }, /* -113 */
{ 0x1.172b83c7d517bp-1, -0x1.19041b9d78a76p-56 }, /* -112 */
{ 0x1.18af9388c8deap-1, -0x1.11023d1970f6cp-55 }, /* -111 */
{ 0x1.1a35beb6fcb75p-1, 0x1.e5b4c7b4968e4p-56 }, /* -110 */
{ 0x1.1bbe084045cd4p-1, -0x1.95386352ef607p-55 }, /* -109 */
{ 0x1.1d4873168b9aap-1, 0x1.e016e00a2643cp-55 }, /* -108 */
{ 0x1.1ed5022fcd91dp-1, -0x1.1df98027bb78cp-55 }, /* -107 */
{ 0x1.2063b88628cd6p-1, 0x1.dc775814a8495p-56 }, /* -106 */
{ 0x1.21f49917ddc96p-1, 0x1.2a97e9494a5eep-56 }, /* -105 */
{ 0x1.2387a6e756238p-1, 0x1.9b07eb6c70573p-55 }, /* -104 */
{ 0x1.251ce4fb2a63fp-1, 0x1.ac155bef4f4a4p-56 }, /* -103 */
{ 0x1.26b4565e27cddp-1, 0x1.2bd339940e9d9p-56 }, /* -102 */
{ 0x1.284dfe1f56381p-1, -0x1.a4c3a8c3f0d7ep-55 }, /* -101 */
{ 0x1.29e9df51fdee1p-1, 0x1.612e8afad1255p-56 }, /* -100 */
{ 0x1.2b87fd0dad99p-1, -0x1.10adcd6381aa4p-60 }, /* -99 */
{ 0x1.2d285a6e4030bp-1, 0x1.0024754db41d5p-55 }, /* -98 */
{ 0x1.2ecafa93e2f56p-1, 0x1.1ca0f45d52383p-57 }, /* -97 */
{ 0x1.306fe0a31b715p-1, 0x1.6f46ad23182e4p-56 }, /* -96 */
{ 0x1.32170fc4cd831p-1, 0x1.a9ce78e18047cp-56 }, /* -95 */
{ 0x1.33c08b26416ffp-1, 0x1.32721843659a6p-55 }, /* -94 */
{ 0x1.356c55f929ff1p-1, -0x1.b5cee5c4e4628p-56 }, /* -93 */
{ 0x1.371a7373aa9cbp-1, -0x1.63aeabf42eae2p-55 }, /* -92 */
{ 0x1.38cae6d05d866p-1, -0x1.e958d3c9904bdp-55 }, /* -91 */
{ 0x1.3a7db34e59ff7p-1, -0x1.5e436d661f5e3p-57 }, /* -90 */
{ 0x1.3c32dc313a8e5p-1, -0x1.efff8375d29c3p-55 }, /* -89 */
{ 0x1.3dea64c123422p-1, 0x1.ada0911f09ebcp-56 }, /* -88 */
{ 0x1.3fa4504ac801cp-1, -0x1.7d023f956f9f3p-55 }, /* -87 */
{ 0x1.4160a21f72e2ap-1, -0x1.ef3691c309278p-59 }, /* -86 */
{ 0x1.431f5d950a897p-1, -0x1.1c7dde35f7999p-56 }, /* -85 */
{ 0x1.44e086061892dp-1, 0x1.89b7a04ef80dp-60 }, /* -84 */
{ 0x1.46a41ed1d0057p-1, 0x1.c944bd1648a76p-55 }, /* -83 */
{ 0x1.486a2b5c13cdp-1, 0x1.3c1a3b69062fp-57 }, /* -82 */
{ 0x1.4a32af0d7d3dep-1, 0x1.9cb62f3d1be56p-55 }, /* -81 */
{ 0x1.4bfdad5362a27p-1, 0x1.d4397afec42e2p-57 }, /* -80 */
{ 0x1.4dcb299fddd0dp-1, 0x1.8ecdbbc6a7833p-55 }, /* -79 */
{ 0x1.4f9b2769d2ca7p-1, -0x1.4b309d25957e3p-55 }, /* -78 */
{ 0x1.516daa2cf6642p-1, -0x1.f768569bd93efp-56 }, /* -77 */
{ 0x1.5342b569d4f82p-1, -0x1.07abe1db13cadp-56 }, /* -76 */
{ 0x1.551a4ca5d920fp-1, -0x1.d689cefede59bp-56 }, /* -75 */
{ 0x1.56f4736b527dap-1, 0x1.9bb2c011d93adp-55 }, /* -74 */
{ 0x1.58d12d497c7fdp-1, 0x1.295e15b9a1de8p-56 }, /* -73 */
{ 0x1.5ab07dd485429p-1, 0x1.6324c054647adp-55 }, /* -72 */
{ 0x1.5c9268a5946b7p-1, 0x1.c4b1b816986a2p-61 }, /* -71 */
{ 0x1.5e76f15ad2148p-1, 0x1.ba6f93080e65ep-55 }, /* -70 */
{ 0x1.605e1b976dc09p-1, -0x1.3e2429b56de47p-55 }, /* -69 */
{ 0x1.6247eb03a5585p-1, -0x1.383c17e40b497p-55 }, /* -68 */
{ 0x1.6434634ccc32p-1, -0x1.c483c759d8933p-56 }, /* -67 */
{ 0x1.6623882552225p-1, -0x1.bb60987591c34p-55 }, /* -66 */
{ 0x1.68155d44ca973p-1, 0x1.038ae44f73e65p-58 }, /* -65 */
{ 0x1.6a09e667f3bcdp-1, -0x1.bdd3413b26456p-55 }, /* -64 */
{ 0x1.6c012750bdabfp-1, -0x1.2895667ff0b0dp-57 }, /* -63 */
{ 0x1.6dfb23c651a2fp-1, -0x1.bbe3a683c88abp-58 }, /* -62 */
{ 0x1.6ff7df9519484p-1, -0x1.83c0f25860ef6p-56 }, /* -61 */
{ 0x1.71f75e8ec5f74p-1, -0x1.16e4786887a99p-56 }, /* -60 */
{ 0x1.73f9a48a58174p-1, -0x1.0a8d96c65d53cp-55 }, /* -59 */
{ 0x1.75feb564267c9p-1, -0x1.0245957316dd3p-55 }, /* -58 */
{ 0x1.780694fde5d3fp-1, 0x1.866b80a02162dp-55 }, /* -57 */
{ 0x1.7a11473eb0187p-1, -0x1.41577ee04992fp-56 }, /* -56 */
{ 0x1.7c1ed0130c132p-1, 0x1.f124cd1164dd6p-55 }, /* -55 */
{ 0x1.7e2f336cf4e62p-1, 0x1.05d02ba15797ep-57 }, /* -54 */
{ 0x1.80427543e1a12p-1, -0x1.27c86626d972bp-55 }, /* -53 */
{ 0x1.82589994cce13p-1, -0x1.d4c1dd41532d8p-55 }, /* -52 */
{ 0x1.8471a4623c7adp-1, -0x1.8d684a341cdfbp-56 }, /* -51 */
{ 0x1.868d99b4492edp-1, -0x1.fc6f89bd4f6bap-55 }, /* -50 */
{ 0x1.88ac7d98a6699p-1, 0x1.994c2f37cb53ap-55 }, /* -49 */
{ 0x1.8ace5422aa0dbp-1, 0x1.6e9f156864b27p-55 }, /* -48 */
{ 0x1.8cf3216b5448cp-1, -0x1.0d55e32e9e3aap-57 }, /* -47 */
{ 0x1.8f1ae99157736p-1, 0x1.5cc13a2e3976cp-56 }, /* -46 */
{ 0x1.9145b0b91ffc6p-1, -0x1.dd6792e582524p-55 }, /* -45 */
{ 0x1.93737b0cdc5e5p-1, -0x1.75fc781b57ebcp-58 }, /* -44 */
{ 0x1.95a44cbc8520fp-1, -0x1.64b7c96a5f039p-57 }, /* -43 */
{ 0x1.97d829fde4e5p-1, -0x1.d185b7c1b85d1p-55 }, /* -42 */
{ 0x1.9a0f170ca07bap-1, -0x1.173bd91cee632p-55 }, /* -41 */
{ 0x1.9c49182a3f09p-1, 0x1.c7c46b071f2bep-57 }, /* -40 */
{ 0x1.9e86319e32323p-1, 0x1.824ca78e64c6ep-57 }, /* -39 */
{ 0x1.a0c667b5de565p-1, -0x1.359495d1cd533p-55 }, /* -38 */
{ 0x1.a309bec4a2d33p-1, 0x1.6305c7ddc36abp-55 }, /* -37 */
{ 0x1.a5503b23e255dp-1, -0x1.d2f6edb8d41e1p-55 }, /* -36 */
{ 0x1.a799e1330b358p-1, 0x1.bcb7ecac563c7p-55 }, /* -35 */
{ 0x1.a9e6b5579fdbfp-1, 0x1.0fac90ef7fd31p-55 }, /* -34 */
{ 0x1.ac36bbfd3f37ap-1, -0x1.f9234cae76cdp-56 }, /* -33 */
{ 0x1.ae89f995ad3adp-1, 0x1.7a1cd345dcc81p-55 }, /* -32 */
{ 0x1.b0e07298db666p-1, -0x1.bdef54c80e425p-55 }, /* -31 */
{ 0x1.b33a2b84f15fbp-1, -0x1.2805e3084d708p-58 }, /* -30 */
{ 0x1.b59728de5593ap-1, -0x1.c71dfbbba6de3p-55 }, /* -29 */
{ 0x1.b7f76f2fb5e47p-1, -0x1.5584f7e54ac3bp-57 }, /* -28 */
{ 0x1.ba5b030a1064ap-1, -0x1.efcd30e54292ep-55 }, /* -27 */
{ 0x1.bcc1e904bc1d2p-1, 0x1.23dd07a2d9e84p-56 }, /* -26 */
{ 0x1.bf2c25bd71e09p-1, -0x1.efdca3f6b9c73p-55 }, /* -25 */
{ 0x1.c199bdd85529cp-1, 0x1.11065895048ddp-56 }, /* -24 */
{ 0x1.c40ab5fffd07ap-1, 0x1.b4537e083c60ap-55 }, /* -23 */
{ 0x1.c67f12e57d14bp-1, 0x1.2884dff483cadp-55 }, /* -22 */
{ 0x1.c8f6d9406e7b5p-1, 0x1.1acbc48805c44p-57 }, /* -21 */
{ 0x1.cb720dcef9069p-1, 0x1.503cbd1e949dbp-57 }, /* -20 */
{ 0x1.cdf0b555dc3fap-1, -0x1.dd83b53829d72p-56 }, /* -19 */
{ 0x1.d072d4a07897cp-1, -0x1.cbc3743797a9cp-55 }, /* -18 */
{ 0x1.d2f87080d89f2p-1, -0x1.d487b719d8578p-55 }, /* -17 */
{ 0x1.d5818dcfba487p-1, 0x1.2ed02d75b3707p-56 }, /* -16 */
{ 0x1.d80e316c98398p-1, -0x1.11ec18beddfe8p-55 }, /* -15 */
{ 0x1.da9e603db3285p-1, 0x1.c2300696db532p-55 }, /* -14 */
{ 0x1.dd321f301b46p-1, 0x1.2da5778f018c3p-55 }, /* -13 */
{ 0x1.dfc97337b9b5fp-1, -0x1.1a5cd4f184b5cp-55 }, /* -12 */
{ 0x1.e264614f5a129p-1, -0x1.7b627817a1496p-55 }, /* -11 */
{ 0x1.e502ee78b3ff6p-1, 0x1.39e8980a9cc8fp-56 }, /* -10 */
{ 0x1.e7a51fbc74c83p-1, 0x1.2d522ca0c8de2p-55 }, /* -9 */
{ 0x1.ea4afa2a490dap-1, -0x1.e9c23179c2893p-55 }, /* -8 */
{ 0x1.ecf482d8e67f1p-1, -0x1.c93f3b411ad8cp-55 }, /* -7 */
{ 0x1.efa1bee615a27p-1, 0x1.dc7f486a4b6bp-55 }, /* -6 */
{ 0x1.f252b376bba97p-1, 0x1.3a1a5bf0d8e43p-55 }, /* -5 */
{ 0x1.f50765b6e454p-1, 0x1.9d3e12dd8a18bp-55 }, /* -4 */
{ 0x1.f7bfdad9cbe14p-1, -0x1.dbb12d006350ap-55 }, /* -3 */
{ 0x1.fa7c1819e90d8p-1, 0x1.74853f3a5931ep-56 }, /* -2 */
{ 0x1.fd3c22b8f71f1p-1, 0x1.2eb74966579e7p-58 }, /* -1 */
{ 0x1p+0, 0x0p+0 }, /* 0 */
{ 0x1.0163da9fb3335p+0, 0x1.b61299ab8cdb7p-54 }, /* 1 */
{ 0x1.02c9a3e778061p+0, -0x1.19083535b085dp-56 }, /* 2 */
{ 0x1.04315e86e7f85p+0, -0x1.0a31c1977c96ep-54 }, /* 3 */
{ 0x1.059b0d3158574p+0, 0x1.d73e2a475b465p-55 }, /* 4 */
{ 0x1.0706b29ddf6dep+0, -0x1.c91dfe2b13c27p-55 }, /* 5 */
{ 0x1.0874518759bc8p+0, 0x1.186be4bb284ffp-57 }, /* 6 */
{ 0x1.09e3ecac6f383p+0, 0x1.1487818316136p-54 }, /* 7 */
{ 0x1.0b5586cf9890fp+0, 0x1.8a62e4adc610bp-54 }, /* 8 */
{ 0x1.0cc922b7247f7p+0, 0x1.01edc16e24f71p-54 }, /* 9 */
{ 0x1.0e3ec32d3d1a2p+0, 0x1.03a1727c57b53p-59 }, /* 10 */
{ 0x1.0fb66affed31bp+0, -0x1.b9bedc44ebd7bp-57 }, /* 11 */
{ 0x1.11301d0125b51p+0, -0x1.6c51039449b3ap-54 }, /* 12 */
{ 0x1.12abdc06c31ccp+0, -0x1.1b514b36ca5c7p-58 }, /* 13 */
{ 0x1.1429aaea92dep+0, -0x1.32fbf9af1369ep-54 }, /* 14 */
{ 0x1.15a98c8a58e51p+0, 0x1.2406ab9eeab0ap-55 }, /* 15 */
{ 0x1.172b83c7d517bp+0, -0x1.19041b9d78a76p-55 }, /* 16 */
{ 0x1.18af9388c8deap+0, -0x1.11023d1970f6cp-54 }, /* 17 */
{ 0x1.1a35beb6fcb75p+0, 0x1.e5b4c7b4968e4p-55 }, /* 18 */
{ 0x1.1bbe084045cd4p+0, -0x1.95386352ef607p-54 }, /* 19 */
{ 0x1.1d4873168b9aap+0, 0x1.e016e00a2643cp-54 }, /* 20 */
{ 0x1.1ed5022fcd91dp+0, -0x1.1df98027bb78cp-54 }, /* 21 */
{ 0x1.2063b88628cd6p+0, 0x1.dc775814a8495p-55 }, /* 22 */
{ 0x1.21f49917ddc96p+0, 0x1.2a97e9494a5eep-55 }, /* 23 */
{ 0x1.2387a6e756238p+0, 0x1.9b07eb6c70573p-54 }, /* 24 */
{ 0x1.251ce4fb2a63fp+0, 0x1.ac155bef4f4a4p-55 }, /* 25 */
{ 0x1.26b4565e27cddp+0, 0x1.2bd339940e9d9p-55 }, /* 26 */
{ 0x1.284dfe1f56381p+0, -0x1.a4c3a8c3f0d7ep-54 }, /* 27 */
{ 0x1.29e9df51fdee1p+0, 0x1.612e8afad1255p-55 }, /* 28 */
{ 0x1.2b87fd0dad99p+0, -0x1.10adcd6381aa4p-59 }, /* 29 */
{ 0x1.2d285a6e4030bp+0, 0x1.0024754db41d5p-54 }, /* 30 */
{ 0x1.2ecafa93e2f56p+0, 0x1.1ca0f45d52383p-56 }, /* 31 */
{ 0x1.306fe0a31b715p+0, 0x1.6f46ad23182e4p-55 }, /* 32 */
{ 0x1.32170fc4cd831p+0, 0x1.a9ce78e18047cp-55 }, /* 33 */
{ 0x1.33c08b26416ffp+0, 0x1.32721843659a6p-54 }, /* 34 */
{ 0x1.356c55f929ff1p+0, -0x1.b5cee5c4e4628p-55 }, /* 35 */
{ 0x1.371a7373aa9cbp+0, -0x1.63aeabf42eae2p-54 }, /* 36 */
{ 0x1.38cae6d05d866p+0, -0x1.e958d3c9904bdp-54 }, /* 37 */
{ 0x1.3a7db34e59ff7p+0, -0x1.5e436d661f5e3p-56 }, /* 38 */
{ 0x1.3c32dc313a8e5p+0, -0x1.efff8375d29c3p-54 }, /* 39 */
{ 0x1.3dea64c123422p+0, 0x1.ada0911f09ebcp-55 }, /* 40 */
{ 0x1.3fa4504ac801cp+0, -0x1.7d023f956f9f3p-54 }, /* 41 */
{ 0x1.4160a21f72e2ap+0, -0x1.ef3691c309278p-58 }, /* 42 */
{ 0x1.431f5d950a897p+0, -0x1.1c7dde35f7999p-55 }, /* 43 */
{ 0x1.44e086061892dp+0, 0x1.89b7a04ef80dp-59 }, /* 44 */
{ 0x1.46a41ed1d0057p+0, 0x1.c944bd1648a76p-54 }, /* 45 */
{ 0x1.486a2b5c13cdp+0, 0x1.3c1a3b69062fp-56 }, /* 46 */
{ 0x1.4a32af0d7d3dep+0, 0x1.9cb62f3d1be56p-54 }, /* 47 */
{ 0x1.4bfdad5362a27p+0, 0x1.d4397afec42e2p-56 }, /* 48 */
{ 0x1.4dcb299fddd0dp+0, 0x1.8ecdbbc6a7833p-54 }, /* 49 */
{ 0x1.4f9b2769d2ca7p+0, -0x1.4b309d25957e3p-54 }, /* 50 */
{ 0x1.516daa2cf6642p+0, -0x1.f768569bd93efp-55 }, /* 51 */
{ 0x1.5342b569d4f82p+0, -0x1.07abe1db13cadp-55 }, /* 52 */
{ 0x1.551a4ca5d920fp+0, -0x1.d689cefede59bp-55 }, /* 53 */
{ 0x1.56f4736b527dap+0, 0x1.9bb2c011d93adp-54 }, /* 54 */
{ 0x1.58d12d497c7fdp+0, 0x1.295e15b9a1de8p-55 }, /* 55 */
{ 0x1.5ab07dd485429p+0, 0x1.6324c054647adp-54 }, /* 56 */
{ 0x1.5c9268a5946b7p+0, 0x1.c4b1b816986a2p-60 }, /* 57 */
{ 0x1.5e76f15ad2148p+0, 0x1.ba6f93080e65ep-54 }, /* 58 */
{ 0x1.605e1b976dc09p+0, -0x1.3e2429b56de47p-54 }, /* 59 */
{ 0x1.6247eb03a5585p+0, -0x1.383c17e40b497p-54 }, /* 60 */
{ 0x1.6434634ccc32p+0, -0x1.c483c759d8933p-55 }, /* 61 */
{ 0x1.6623882552225p+0, -0x1.bb60987591c34p-54 }, /* 62 */
{ 0x1.68155d44ca973p+0, 0x1.038ae44f73e65p-57 }, /* 63 */
{ 0x1.6a09e667f3bcdp+0, -0x1.bdd3413b26456p-54 }, /* 64 */
{ 0x1.6c012750bdabfp+0, -0x1.2895667ff0b0dp-56 }, /* 65 */
{ 0x1.6dfb23c651a2fp+0, -0x1.bbe3a683c88abp-57 }, /* 66 */
{ 0x1.6ff7df9519484p+0, -0x1.83c0f25860ef6p-55 }, /* 67 */
{ 0x1.71f75e8ec5f74p+0, -0x1.16e4786887a99p-55 }, /* 68 */
{ 0x1.73f9a48a58174p+0, -0x1.0a8d96c65d53cp-54 }, /* 69 */
{ 0x1.75feb564267c9p+0, -0x1.0245957316dd3p-54 }, /* 70 */
{ 0x1.780694fde5d3fp+0, 0x1.866b80a02162dp-54 }, /* 71 */
{ 0x1.7a11473eb0187p+0, -0x1.41577ee04992fp-55 }, /* 72 */
{ 0x1.7c1ed0130c132p+0, 0x1.f124cd1164dd6p-54 }, /* 73 */
{ 0x1.7e2f336cf4e62p+0, 0x1.05d02ba15797ep-56 }, /* 74 */
{ 0x1.80427543e1a12p+0, -0x1.27c86626d972bp-54 }, /* 75 */
{ 0x1.82589994cce13p+0, -0x1.d4c1dd41532d8p-54 }, /* 76 */
{ 0x1.8471a4623c7adp+0, -0x1.8d684a341cdfbp-55 }, /* 77 */
{ 0x1.868d99b4492edp+0, -0x1.fc6f89bd4f6bap-54 }, /* 78 */
{ 0x1.88ac7d98a6699p+0, 0x1.994c2f37cb53ap-54 }, /* 79 */
{ 0x1.8ace5422aa0dbp+0, 0x1.6e9f156864b27p-54 }, /* 80 */
{ 0x1.8cf3216b5448cp+0, -0x1.0d55e32e9e3aap-56 }, /* 81 */
{ 0x1.8f1ae99157736p+0, 0x1.5cc13a2e3976cp-55 }, /* 82 */
{ 0x1.9145b0b91ffc6p+0, -0x1.dd6792e582524p-54 }, /* 83 */
{ 0x1.93737b0cdc5e5p+0, -0x1.75fc781b57ebcp-57 }, /* 84 */
{ 0x1.95a44cbc8520fp+0, -0x1.64b7c96a5f039p-56 }, /* 85 */
{ 0x1.97d829fde4e5p+0, -0x1.d185b7c1b85d1p-54 }, /* 86 */
{ 0x1.9a0f170ca07bap+0, -0x1.173bd91cee632p-54 }, /* 87 */
{ 0x1.9c49182a3f09p+0, 0x1.c7c46b071f2bep-56 }, /* 88 */
{ 0x1.9e86319e32323p+0, 0x1.824ca78e64c6ep-56 }, /* 89 */
{ 0x1.a0c667b5de565p+0, -0x1.359495d1cd533p-54 }, /* 90 */
{ 0x1.a309bec4a2d33p+0, 0x1.6305c7ddc36abp-54 }, /* 91 */
{ 0x1.a5503b23e255dp+0, -0x1.d2f6edb8d41e1p-54 }, /* 92 */
{ 0x1.a799e1330b358p+0, 0x1.bcb7ecac563c7p-54 }, /* 93 */
{ 0x1.a9e6b5579fdbfp+0, 0x1.0fac90ef7fd31p-54 }, /* 94 */
{ 0x1.ac36bbfd3f37ap+0, -0x1.f9234cae76cdp-55 }, /* 95 */
{ 0x1.ae89f995ad3adp+0, 0x1.7a1cd345dcc81p-54 }, /* 96 */
{ 0x1.b0e07298db666p+0, -0x1.bdef54c80e425p-54 }, /* 97 */
{ 0x1.b33a2b84f15fbp+0, -0x1.2805e3084d708p-57 }, /* 98 */
{ 0x1.b59728de5593ap+0, -0x1.c71dfbbba6de3p-54 }, /* 99 */
{ 0x1.b7f76f2fb5e47p+0, -0x1.5584f7e54ac3bp-56 }, /* 100 */
{ 0x1.ba5b030a1064ap+0, -0x1.efcd30e54292ep-54 }, /* 101 */
{ 0x1.bcc1e904bc1d2p+0, 0x1.23dd07a2d9e84p-55 }, /* 102 */
{ 0x1.bf2c25bd71e09p+0, -0x1.efdca3f6b9c73p-54 }, /* 103 */
{ 0x1.c199bdd85529cp+0, 0x1.11065895048ddp-55 }, /* 104 */
{ 0x1.c40ab5fffd07ap+0, 0x1.b4537e083c60ap-54 }, /* 105 */
{ 0x1.c67f12e57d14bp+0, 0x1.2884dff483cadp-54 }, /* 106 */
{ 0x1.c8f6d9406e7b5p+0, 0x1.1acbc48805c44p-56 }, /* 107 */
{ 0x1.cb720dcef9069p+0, 0x1.503cbd1e949dbp-56 }, /* 108 */
{ 0x1.cdf0b555dc3fap+0, -0x1.dd83b53829d72p-55 }, /* 109 */
{ 0x1.d072d4a07897cp+0, -0x1.cbc3743797a9cp-54 }, /* 110 */
{ 0x1.d2f87080d89f2p+0, -0x1.d487b719d8578p-54 }, /* 111 */
{ 0x1.d5818dcfba487p+0, 0x1.2ed02d75b3707p-55 }, /* 112 */
{ 0x1.d80e316c98398p+0, -0x1.11ec18beddfe8p-54 }, /* 113 */
{ 0x1.da9e603db3285p+0, 0x1.c2300696db532p-54 }, /* 114 */
{ 0x1.dd321f301b46p+0, 0x1.2da5778f018c3p-54 }, /* 115 */
{ 0x1.dfc97337b9b5fp+0, -0x1.1a5cd4f184b5cp-54 }, /* 116 */
{ 0x1.e264614f5a129p+0, -0x1.7b627817a1496p-54 }, /* 117 */
{ 0x1.e502ee78b3ff6p+0, 0x1.39e8980a9cc8fp-55 }, /* 118 */
{ 0x1.e7a51fbc74c83p+0, 0x1.2d522ca0c8de2p-54 }, /* 119 */
{ 0x1.ea4afa2a490dap+0, -0x1.e9c23179c2893p-54 }, /* 120 */
{ 0x1.ecf482d8e67f1p+0, -0x1.c93f3b411ad8cp-54 }, /* 121 */
{ 0x1.efa1bee615a27p+0, 0x1.dc7f486a4b6bp-54 }, /* 122 */
{ 0x1.f252b376bba97p+0, 0x1.3a1a5bf0d8e43p-54 }, /* 123 */
{ 0x1.f50765b6e454p+0, 0x1.9d3e12dd8a18bp-54 }, /* 124 */
{ 0x1.f7bfdad9cbe14p+0, -0x1.dbb12d006350ap-54 }, /* 125 */
{ 0x1.fa7c1819e90d8p+0, 0x1.74853f3a5931ep-55 }, /* 126 */
{ 0x1.fd3c22b8f71f1p+0, 0x1.2eb74966579e7p-57 }, /* 127 */
};

static const double xmax = 0x1p1023;

#define MAYBE_UNUSED __attribute__ ((unused))

/* compute a double-double approximation of
   exp(x) ~ 2^e * 2^(i/128) * 2^h * 2^l
   where -127 <= i <= 127, |h| < 1/128, and |l| < 2^-42.
   We use a degree-9 polynomial to approximate 2^h on [-1/128,1/128]
   (cf exp2_acc.sollya) with 105.765 bits of relative precision.
   Coefficients of degree 0-4 are double-double, 5-9 are doubles. */
static double
cr_exp_accurate (double x, int e, int i)
{
  double h, l, yh, yl;

  /* we recompute h, l such that x/log(2) = e + i/128 + h + l,
     to get more accuracy on l (in the fast path we extract the high 8 bits
     of h to go into i, thus we lose 8 bits on l) */
  static const double log2_h = 0x1.71547652b82fep+0;
  static const double log2_m = 0x1.777d0ffda0ep-56;
  static const double log2_l = -0x1.b8b05dc52a215p-101;
  double l1;
  double ah, al, bh, bl;
  dekker (&ah, &al, x, log2_h);
  dekker (&bh, &bl, x, log2_m);
  ah -= e + i / 128.0;
  /* h + l = ah + (al + bh) + (bl + c) */
  fast_two_sum (&h, &l, ah, al);
  fast_two_sum (&h, &l1, h, bh);
  /* h + (l+l1) = ah + (al + bh) */
  l += l1 + (bl + x * log2_l);

  /* now x/log(2) ~ e + i/128 + h + l */

  static const double p[14] = { 0x1p0,            /* p[0]: degree 0 */
    0x1.62e42fefa39efp-1, 0x1.abc9e3b397ebp-56,   /* p[1,2]: degree 1 */
    0x1.ebfbdff82c58fp-3, -0x1.5e43a5429b326p-57, /* p[3,4]: degree 2 */
    0x1.c6b08d704a0cp-5, -0x1.d331600cee073p-59,  /* p[5,6]: degree 3 */
    0x1.3b2ab6fba4e77p-7, 0x1.4fb30e5c2c8bcp-62,  /* p[7,8]: degree 4 */
    0x1.5d87fe78a6731p-10,                        /* p[9]: degree 5 */
    0x1.430912f86bfb8p-13,                        /* p[10]: degree 6 */
    0x1.ffcbfc58b51c9p-17,                        /* p[11]: degree 7 */
    0x1.62c034be4ffd9p-20,                        /* p[12]: degree 8 */
    0x1.b523023e3d552p-24 };                      /* p[13]: degree 9 */
  double t, u;
  yh = p[12] + h * p[13];
  yh = p[11] + h * yh;
  yh = p[10] + h * yh;
  /* ulp(p6*h^6) <= 2^-107, thus we can compute h*yh in double,
     but ulp(p5*h^5) <= 2^-97, so we need a double-double to store
     p5 + h*yh */
  fast_two_sum (&yh, &yl, p[9], h * yh);
  /* multiply (yh,yl) by h and add p4=(p[7],p[8]) */
  dekker (&t, &u, yh, h);
  u += yl * h;
  fast_two_sum (&yh, &yl, p[7], t);
  yl += u + p[8];
  /* multiply (yh,yl) by h and add p3=(p[5],p[6]) */
  dekker (&t, &u, yh, h);
  u += yl * h;
  fast_two_sum (&yh, &yl, p[5], t);
  yl += u + p[6];
  /* multiply (yh,yl) by h and add p2=(p[3],p[4]) */
  dekker (&t, &u, yh, h);
  u += yl * h;
  fast_two_sum (&yh, &yl, p[3], t);
  yl += u + p[4];
  /* multiply (yh,yl) by h and add p1=(p[1],p[2]) */
  dekker (&t, &u, yh, h);
  u += yl * h;
  fast_two_sum (&yh, &yl, p[1], t);
  yl += u + p[2];
  /* multiply (yh,yl) by h and add p0=1 */
  dekker (&t, &u, yh, h);
  u += yl * h;
  fast_two_sum (&yh, &yl, p[0], t);
  yl += u;
  /* yh+yl is an approximation of 2^h */

  /* multiply by 2^l */
  static const double l2 = 0x1.62e42fefa39efp-1;
  t = l2 * l * yh;
  fast_two_sum (&yh, &u, yh, t);
  u += yl;

  /* multiply by 2^(i/128) */
  t = yh * tab_i[127+i][1];
  dekker (&yh, &yl, yh, tab_i[127+i][0]);
  yl += t + u * tab_i[127+i][0];
  d64u64 v;
  v.x = yh + yl;
  unsigned int f = v.n >> 52;
  f += e;
  if (__builtin_expect((f - 1) > 0x7ff, 0)) /* overflow or subnormal */
  {
    if (e > 0) /* overflow */
      return xmax + xmax;
    /* Same method as in the fast path. */
    double magic = ldexp (1.0, -1022 - e);
    fast_two_sum (&t, &u, magic, yh);
    t += u + yl; /* here comes the rounding */
    t -= magic;
    return ldexp (t, e);
  }

  d64u64 w;
  w.x = x;
  switch (w.n & 0x3f) {
    case 0:
      if (x == 0x1p-53)
        return 0x1.0000000000001p+0 - 0x1.fffffffffffffp-54;
      if (x == -0x1.60296a66b43p+8)
        return 0x1.ea71d85cee02p-509 - 0x1.de4ec450e71c7p-614;
      break;
    case 37:
      if (x == 0x1.00091a4a0dae5p+2)
        return 0x1.b50726fd1ccb3p+5 - 0x1.ffffffffffffep-49;
      if (x == -0x1.ff171507f8ba5p-30)
        return 0x1.fffffff007475p-1 + 0x1.fffffffffffffp-55;
      if (x == 0x1.d7a7d893609e5p-26)
        return 0x1.00000075e9f64p+0 + 0x1.63d01ed72ff7bp-109;
      if (x == 0x1.d7a7d893609e5p-26)
        return 0x1.00000075e9f64p+0 + 0x1.63d01ed72ff7bp-109;
      break;
    case 43:
      if (x == 0x1.005ae04256babp-1)
        return 0x1.a65d89abf3d1fp+0 - 0x1.fffffffffffffp-54;
      break;
    case 41:
      if (x == -0x1.02393d5976769p+1)
        return 0x1.1064b2c103ddbp-3 - 0x1.fffffffffffffp-57;
      break;
    case 29:
      if (x == -0x1.150dd1bf6ae1dp+3)
        return 0x1.6c5d3ae7c88f7p-13 - 0x1.ffffffffffffep-67;
      if (x == 0x1.7a60ee15e3e9dp+6)
        return 0x1.62e4dc3bbf53fp+136 + 0x1.ae7c8eddb6bcbp+30;
      if (x == 0x1.7a60ee15e3e9dp+6)
        return 0x1.62e4dc3bbf53fp+136 + 0x1.ae7c8eddb6bcbp+30;
      break;
    case 20:
      if (x == 0x1.273c188aa7b14p+2)
        return 0x1.93295a96ec6ebp+6 - 0x1.fffffffffffffp-48;
      if (x == 0x1.369bc71a44b14p-22)
        return 0x1.000004da6f283p+0 - 0x1.f9bc07bd14aeap-107;
      if (x == 0x1.458f7365fd894p-8)
        return 0x1.01465ece08736p+0 - 0x1.092f785a56fd9p-106;
      break;
    case 39:
      if (x == 0x1.282aae3169ee7p-20)
        return 0x1.00001282ab8e7p+0 - 0x1.fffffffffffffp-54;
      if (x == -0x1.7fb235d76cce7p-8)
        return 0x1.fd02d98c24bbbp-1 - 0x1.fffffffffffffp-55;
      if (x == 0x1.ba07d73250de7p-14)
        return 0x1.0006e83736f8dp+0 - 0x1.fffffffffffffp-54;
      if (x == 0x1.79c2a7777bee7p-10)
        return 0x1.005e8217799bbp+0 + 0x1.500dede5d5458p-106;
      if (x == 0x1.1d5c2daebe367p+4)
        return 0x1.a8c02e974c315p+25 - 0x1.de0fc9395bbd4p-83;
      if (x == 0x1.79c2a7777bee7p-10)
        return 0x1.005e8217799bbp+0 + 0x1.500dede5d5458p-106;
      break;
    case 7:
      if (x == 0x1.413f33fe6ce07p-31)
        return 0x1.00000002827e7p+0 - 0x1.fffffffffffffp-54;
      break;
    case 12:
      if (x == -0x1.54511e930898cp-7)
        return 0x1.fab5c6e464e0dp-1 + 0x1.ffffffffffffdp-55;
      if (x == 0x1.7d7fc2e4f5fccp-3)
        return 0x1.346b07b6dd7f3p+0 - 0x1.fffffffffffffp-54;
      if (x == -0x1.057b641debc0cp-15)
        return 0x1.fffbea169bd8dp-1 + 0x1.da08d6f93c81ep-108;
      if (x == -0x1.8f80e06f3a04cp+4)
        return 0x1.f80aafa92b498p-37 - 0x1.191595b7e88edp-143;
      break;
    case 28:
      if (x == -0x1.59f038076039cp+6)
        return 0x1.2c0fa76a0e15fp-125 + 0x1.fffffffffffffp-179;
      if (x == 0x1.8cd99fffb319cp-34)
        return 0x1.0000000063367p+0 - 0x1.fffffffffffffp-54;
      if (x == -0x1.85e60704a3a9cp-30)
        return 0x1.fffffff3d0cfdp-1 - 0x1.fffffffffffffp-55;
      if (x == 0x1.51d0f4f0a901cp+5)
        return 0x1.e4a01c9ddbc87p+60 + 0x1.676fbe202d089p-44;
      if (x == 0x1.51d0f4f0a901cp+5)
        return 0x1.e4a01c9ddbc87p+60 + 0x1.676fbe202d089p-44;
      break;
    case 34:
      if (x == 0x1.6b7fef325ada2p+2)
        return 0x1.24db524842d71p+8 - 0x1.ffffffffffffep-46;
      if (x == 0x1.1b231819a1ce2p-4)
        return 0x1.12527095bbd08p+0 + 0x1.55e8d3c024021p-106;
      if (x == -0x1.2a9cad9998262p+0)
        return 0x1.3ef1e9b3a81c8p-2 - 0x1.e757fe830d60ep-109;
      if (x == 0x1.32d3b310ba562p-14)
        return 0x1.0004cb5a4a4b8p+0 - 0x1.6a1e33d1252d9p-108;
      if (x == 0x1.1b231819a1ce2p-4)
        return 0x1.12527095bbd08p+0 + 0x1.55e8d3c024021p-106;
      break;
    case 54:
      if (x == -0x1.85068c07fbbf6p-1)
        return 0x1.defa8f4a8af21p-2 - 0x1.ffffffffffffep-56;
      break;
    case 57:
      if (x == -0x1.a4187f2ca71f9p-6)
        return 0x1.f309f46111221p-1 - 0x1.ffffffffffffep-55;
      if (x == -0x1.290ea09e36479p-3)
        return 0x1.baded30cbf1c4p-1 - 0x1.09bec3f4113eep-111;
      if (x == 0x1.d8d934d46a2b9p-7)
        return 0x1.03b88d9bb35e2p+0 - 0x1.061ccfb697d4ap-104;
      break;
    case 1:
      if (x == -0x1.ac026dccaa781p+5)
        return 0x1.c219cd8029075p-78 - 0x1.ffffffffffffdp-132;
      if (x == 0x1.ffff7fffe0001p-36)
        return 0x1.000000001ffffp+0 + 0x1.fffffffffffffp-54;
      if (x == -0x1.cc37ef7de7501p+0)
        return 0x1.534d4de870713p-3 + 0x1.e37dbca61347dp-110;
      if (x == -0x1.cc37ef7de7501p+0)
        return 0x1.534d4de870713p-3 + 0x1.e37dbca61347dp-110;
      break;
    case 59:
      if (x == 0x1.aca7ae8da5a7bp+0)
        return 0x1.557d4acd7e557p+2 - 0x1.fffffffffffffp-52;
      if (x == 0x1.be2caeebfc83bp-4)
        return 0x1.1d761d8637563p+0 + 0x1.a3cd02c39fb3ep-106;
      if (x == 0x1.16e7724ac74bbp-19)
        return 0x1.000022dcf0a91p+0 - 0x1.7a2300c13e9cfp-107;
      if (x == 0x1.be2caeebfc83bp-4)
        return 0x1.1d761d8637563p+0 + 0x1.a3cd02c39fb3ep-106;
      break;
    case 48:
      if (x == 0x1.accfbe46b4efp-1)
        return 0x1.27c2e4bc1ee7p+1 + 0x1.fffffffffffffp-53;
      break;
    case 30:
      if (x == 0x1.bcab27d05abdep-2)
        return 0x1.8b367381d82f5p+0 - 0x1.ffffffffffffep-54;
      if (x == 0x1.9c026d62f631ep-11)
        return 0x1.0033857c34f23p+0 + 0x1.1be674ce71eccp-107;
      if (x == 0x1.9c026d62f631ep-11)
        return 0x1.0033857c34f23p+0 + 0x1.1be674ce71eccp-107;
      if (x == -0x1.cc83748b7669ep-2)
        return 0x1.468e956d45383p-1 - 0x1.05665296469fdp-106;
      if (x == 0x1.9c026d62f631ep-11)
        return 0x1.0033857c34f23p+0 + 0x1.1be674ce71eccp-107;
      break;
    case 26:
      if (x == -0x1.c794ddcbd661ap+3)
        return 0x1.6040718b6b7cdp-21 - 0x1.fffffffffffffp-75;
      if (x == -0x1.daf693d64fadap-4)
        return 0x1.c7f14af0a08ebp-1 - 0x1.2aa7e4ef70195p-109;
      break;
    case 19:
      if (x == -0x1.d1d85a7f253d3p-15)
        return 0x1.fff8b8abd4c21p-1 + 0x1.fffffffffffffp-55;
      if (x == -0x1.8770955e29c93p+4)
        return 0x1.a12b15f39cf17p-36 + 0x1.556e099f7c80bp-140;
      if (x == -0x1.8770955e29c93p+4)
        return 0x1.a12b15f39cf17p-36 + 0x1.556e099f7c80bp-140;
      break;
    case 21:
      if (x == -0x1.d3f3799439415p-3)
        return 0x1.976a4c9985f5bp-1 + 0x1.fffffffffffffp-55;
      break;
    case 47:
      if (x == 0x1.e5ba92aa9b1efp-20)
        return 0x1.00001e5baaf77p+0 + 0x1.fffffffffffffp-54;
      if (x == -0x1.4bd46601ae1efp-31)
        return 0x1.fffffffad0ae7p-1 - 0x1.fffffffffffffp-55;
      if (x == -0x1.0a54d87783d6fp+0)
        return 0x1.69cef05657108p-2 + 0x1.d81f352752164p-108;
      if (x == 0x1.e07e71bfcf06fp+5)
        return 0x1.91ec4412c344fp+86 + 0x1.09d2b56d79a72p-23;
      if (x == -0x1.cddf723d3e52fp-3)
        return 0x1.98a04e0833091p-1 - 0x1.4c519851f4cf7p-106;
      if (x == 0x1.e07e71bfcf06fp+5)
        return 0x1.91ec4412c344fp+86 + 0x1.09d2b56d79a72p-23;
      if (x == -0x1.0a54d87783d6fp+0)
        return 0x1.69cef05657108p-2 + 0x1.d81f352752164p-108;
      if (x == 0x1.e07e71bfcf06fp+5)
        return 0x1.91ec4412c344fp+86 + 0x1.09d2b56d79a72p-23;
      break;
    case 55:
      if (x == 0x1.f8e165f8388f7p-30)
        return 0x1.00000007e3859p+0 + 0x1.fffffffffffffp-54;
      break;
    case 53:
      if (x == -0x1.0401ae48409b5p-28)
        return 0x1.ffffffdf7fca3p-1 + 0x1.fffffffffffffp-55;
      if (x == -0x1.8aeb636f3ce35p-3)
        return 0x1.a634ae87df6aep-1 + 0x1.82b6b66e03876p-110;
      break;
    case 3:
      if (x == 0x1.08f51434652c3p+4)
        return 0x1.daac459b157e5p+23 + 0x1.c6823badae774p-84;
      if (x == 0x1.08f51434652c3p+4)
        return 0x1.daac459b157e5p+23 + 0x1.c6823badae774p-84;
      break;
    case 62:
      if (x == -0x1.1ff9b8e8b38bep-7)
        return 0x1.fb85251a3f26fp-1 + 0x1.08e2a7c5ac444p-109;
      if (x == 0x1.333a83013057ep+2)
        return 0x1.e642354c34a34p+6 + 0x1.3df4a4fb0b639p-99;
      if (x == -0x1.1ff9b8e8b38bep-7)
        return 0x1.fb85251a3f26fp-1 + 0x1.08e2a7c5ac444p-109;
      if (x == 0x1.333a83013057ep+2)
        return 0x1.e642354c34a34p+6 + 0x1.3df4a4fb0b639p-99;
      break;
    case 25:
      if (x == -0x1.381126525f9d9p-7)
        return 0x1.fb25a83c4532p-1 + 0x1.34b866e39d79dp-106;
      if (x == 0x1.d707029bb59d9p-2)
        return 0x1.958497f7b353fp+0 + 0x1.e27333f64e0b9p-106;
      if (x == -0x1.381126525f9d9p-7)
        return 0x1.fb25a83c4532p-1 + 0x1.34b866e39d79dp-106;
      if (x == 0x1.d707029bb59d9p-2)
        return 0x1.958497f7b353fp+0 + 0x1.e27333f64e0b9p-106;
      break;
    case 16:
      if (x == 0x1.39fc4d3bb711p+5)
        return 0x1.8a4e90733b95ep+56 + 0x1.6fc2e81137fa9p-49;
      if (x == 0x1.39fc4d3bb711p+5)
        return 0x1.8a4e90733b95ep+56 + 0x1.6fc2e81137fa9p-49;
      break;
    case 27:
      if (x == 0x1.46370d915991bp-4)
        return 0x1.1538ea18a4585p+0 + 0x1.cf02a409a72c9p-107;
      if (x == 0x1.46370d915991bp-4)
        return 0x1.1538ea18a4585p+0 + 0x1.cf02a409a72c9p-107;
      break;
    case 58:
      if (x == 0x1.54cd1fea7663ap+7)
        return 0x1.c90810d354618p+245 + 0x1.2925a9627fb2cp+136;
      if (x == 0x1.54cd1fea7663ap+7)
        return 0x1.c90810d354618p+245 + 0x1.2925a9627fb2cp+136;
      break;
    case 33:
      if (x == 0x1.62f71c4656b61p-1)
        return 0x1.000976581ce4ep+1 + 0x1.01dc6b104a893p-105;
      if (x == 0x1.62f71c4656b61p-1)
        return 0x1.000976581ce4ep+1 + 0x1.01dc6b104a893p-105;
      break;
    case 4:
      if (x == 0x1.7486e56dbca44p-14)
        return 0x1.0005d22c869a6p+0 + 0x1.878d58ee35bc1p-107;
      if (x == -0x1.4156584bcd084p+7)
        return 0x1.26e9c4d32796p-232 + 0x1.566f8d0f3440dp-339;
      if (x == 0x1.7486e56dbca44p-14)
        return 0x1.0005d22c869a6p+0 + 0x1.878d58ee35bc1p-107;
      break;
    case 36:
      if (x == -0x1.78ba0840d79e4p+4)
        return 0x1.05960be55d55ap-34 + 0x1.2b040838e056bp-139;
      if (x == 0x1.c7206c1b753e4p+8)
        return 0x1.8670de0b68cadp+656 - 0x1.7599cebd802f7p+549;
      break;
    case 52:
      if (x == 0x1.83d4bcdebb3f4p+2)
        return 0x1.ac50b409c8aeep+8 + 0x1.16719fcede453p-103;
      if (x == 0x1.b9d483946ef74p+2)
        return 0x1.f1ecb20dbd6d9p+9 + 0x1.90a0f3bbcc087p-94;
      if (x == 0x1.83d4bcdebb3f4p+2)
        return 0x1.ac50b409c8aeep+8 + 0x1.16719fcede453p-103;
      if (x == 0x1.b9d483946ef74p+2)
        return 0x1.f1ecb20dbd6d9p+9 + 0x1.90a0f3bbcc087p-94;
      break;
    case 63:
      if (x == 0x1.97efcad0968ffp-4)
        return 0x1.1acf134ae11cfp+0 + 0x1.58240d09ab4abp-105;
      if (x == 0x1.d77fd13d27fffp-11)
        return 0x1.003af6c37c1d3p+0 + 0x1.288ddd3991ec4p-109;
      if (x == 0x1.97efcad0968ffp-4)
        return 0x1.1acf134ae11cfp+0 + 0x1.58240d09ab4abp-105;
      if (x == 0x1.d77fd13d27fffp-11)
        return 0x1.003af6c37c1d3p+0 + 0x1.288ddd3991ec4p-109;
      break;
    case 14:
      if (x == -0x1.a01c846cb4a4ep-15)
        return 0x1.fff97f987fb48p-1 + 0x1.28458fcaa6d88p-110;
      if (x == 0x1.c1b4603b2f18ep-8)
        return 0x1.01c3404506a2cp+0 + 0x1.1f8cef4615c73p-105;
      if (x == -0x1.e8bdbfcd9144ep+3)
        return 0x1.f3e558cf4de54p-23 + 0x1.f5808cf1d9ee5p-131;
      if (x == 0x1.46d06e5d90f4ep-22)
        return 0x1.0000051b41c68p+0 - 0x1.90412d7b31c2cp-107;
      if (x == -0x1.a01c846cb4a4ep-15)
        return 0x1.fff97f987fb48p-1 + 0x1.28458fcaa6d88p-110;
      if (x == 0x1.c1b4603b2f18ep-8)
        return 0x1.01c3404506a2cp+0 + 0x1.1f8cef4615c73p-105;
      if (x == -0x1.e8bdbfcd9144ep+3)
        return 0x1.f3e558cf4de54p-23 + 0x1.f5808cf1d9ee5p-131;
      break;
    case 31:
      if (x == -0x1.a2fefefd580dfp-13)
        return 0x1.ffe5d0bb7eabfp-1 + 0x1.81905a57d6a44p-112;
      if (x == -0x1.a2fefefd580dfp-13)
        return 0x1.ffe5d0bb7eabfp-1 + 0x1.81905a57d6a44p-112;
      break;
    case 45:
      if (x == -0x1.b604e4b77d96dp-5)
        return 0x1.e557c34fef2b1p-1 + 0x1.5822731ec52fbp-106;
      break;
    case 11:
      if (x == -0x1.cecc4ad8d358bp-8)
        return 0x1.fc65aa1908a66p-1 + 0x1.631ff977c5109p-107;
      if (x == 0x1.26ee1a46d8c8bp+9)
        return 0x1.fbe20477df4a7p+850 - 0x1.556f0ed19479ep+746;
      if (x == -0x1.cecc4ad8d358bp-8)
        return 0x1.fc65aa1908a66p-1 + 0x1.631ff977c5109p-107;
      break;
    case 42:
      if (x == 0x1.d6336a88077aap+0)
        return 0x1.91a8dff540ff7p+2 + 0x1.78f1982b593afp-105;
      if (x == 0x1.d6336a88077aap+0)
        return 0x1.91a8dff540ff7p+2 + 0x1.78f1982b593afp-105;
      break;
    case 50:
      if (x == 0x1.f6e4c3ced7c72p-3)
        return 0x1.47408cb9583cep+0 + 0x1.644b7f5399dfp-107;
      if (x == 0x1.f6e4c3ced7c72p-3)
        return 0x1.47408cb9583cep+0 + 0x1.644b7f5399dfp-107;
      break;
    case 5:
      if (x == 0x1.4ff7529737745p-7)
        return 0x1.02a3637d2021fp+0 - 0x1.8ce83a1239239p-105;
      break;
    case 23:
      if (x == 0x1.6474c604cc0d7p+6)
        return 0x1.7a8f65ad009bdp+128 - 0x1.0b611158ec877p+21;
      break;
    case 46:
      if (x == 0x1.968d2004d83eep-6)
        return 0x1.066e8c9ca7589p+0 - 0x1.fc9cb62e1cf95p-105;
      break;
    case 17:
      if (x == -0x1.a2772020006d1p-12)
        return 0x1.ffcbb3c7edf17p-1 - 0x1.41d718c5b866ep-109;
      break;
    case 2:
      if (x == 0x1.dd2ae021a5cc2p-25)
        return 0x1.000000ee95708p+0 - 0x1.517fb684e390ap-106;
      break;
    case 8:
      if (x == 0x1.eb9914d4ac1c8p+8)
        return 0x1.2b67eff65dce8p+709 - 0x1.01c4a555ef227p+604;
      break;
    case 35:
      if (x == -0x1.f39c6d7449a23p+2)
        return 0x1.aae341b425b8bp-12 - 0x1.428a476d4fcabp-117;
      break;
    case 10:
      if (x == 0x1.1a0408712e00ap-2)
        return 0x1.512b3126454f3p+0 + 0x1.758d621c3b9dep-106;
      break;
  };

  /* special code for tiny x */
  if (__builtin_fabs (x) < 0x1p-27)
  {
    /* example: x=0x1.051699fdeb71dp-30: 49 identical bits after round bit */
    /* -0x1.cfd2c0003485cp-35: 53 */
    /* -0x1.ff45860ff45a8p-29: 44 */
    h = x * 0x1.5555555555555p-3; /* x/6 */
    fast_two_sum (&yh, &yl, 0.5, h); /* yh+l ~ 1/2 + x/6 */
    dekker (&h, &l, yh, x);
    l += yl * x;                     /* h+l ~ x/2 + x^2/6 */
    fast_two_sum (&yh, &yl, 1.0, h);
    yl += l;                         /* yh+yl ~ 1 + x/2+x^2/6 */
    dekker (&h, &l, yh, x);
    l += yl * x;                     /* h+l ~ x + x^2/2 + x^3/6 */
    fast_two_sum (&yh, &yl, 1.0, h);
    yl += l;                         /* yh+yl ~ 1 + x + x^2/2 + x^3/6 */
    return yh + yl;
  }

  return ldexp (v.x, e);
}

double
cr_exp (double x)
{
  d64u64 v;
  v.x = x;
  int e = ((v.n >> 52) & 0x7ff) - 0x3ff;

  if (e >= 9) /* potential underflow or overflow */
  {
    if (isnan (x))
      return x + x; /* always return qNaN, even for sNaN input */

    if (x >= 0x1.62e42fefa39fp+9) /* exp(x) > 2^1024*(1-2^-54) */
      return xmax + xmax;
    else if (x <= -0x1.74910d52d3052p+9) /* exp(x) < 2^-1075 */
    {
      static const double xmin = 0x1p-1074;
      return xmin / 2.0;
    }
    /* otherwise go through the main path */
  }

  /* now |x| < 746 */

  /* first multiply x by a double-double approximation of 1/log(2) */
  static const double log2_h = 0x1.71547652b82fep+0;
  static const double log2_l = 0x1.777d0ffda0ep-56;
  /* |1/log(2) - log2_h - log2_l| < 2^-100.2 */
  double h, l, l0;
  dekker (&h, &l0, x, log2_h);
  /* h + l0 = x * log2_h exactly */
  l = l0 + x * log2_l;
  /* |x * log2_l| < 2^-45 thus the error on x * log2_h (even without fma)
     is bounded by 2^-98.
     |h| < 2^11 thus |l0| < 2^-42 thus the rounding error on l is bounded
     by 2^-95.
     The total rounding error on l is bounded by 2^-98+2^-95 < 2^-94.83.
     |x/log(2) - h - l| < |x| * 2^-100.2 + 2^-94.83 < 2^-90.57. */

  /* now x/log(2) ~ h + l thus exp(x) ~ 2^h * 2^l where |l| < 2^-42 */
  double t = __builtin_trunc (128.0 * h), u;
  h = h - t / 128.0; /* exact */
  e = t;
  int i = e % 128;
  e = (e - i) >> 7;
  /* exp(x) ~ 2^e * 2^(i/128) * 2^h * 2^l where |h| < 1/128 and |l| < 2^-42,
     where -127 <= i <= 127 */

  /* p[i] are the coefficients of a degree-6 polynomial approximating 2^x
     over [-1/128,1/128], with double coefficients, except p[1] which is
     double-double (here we give only the upper term, the lower is p1l.
     The relative accuracy is 70.862 bits. */
  static const double p[7] = {
  0x1p0, 0x1.62e42fefa39efp-1, 0x1.ebfbdff82c58fp-3,
    0x1.c6b08d70484c1p-5, 0x1.3b2ab6fb663a2p-7, 0x1.5d881a764d899p-10,
    0x1.430bba9c70dddp-13 };
  static const double p1l = 0x1.b2ca0bb577094p-56;
  double hh = h * h;
  /* |h| < 2^-7 thus hh < 2^-14, and the error on hh is bounded by 2^-67 */
  double yl = p[5] + h * p[6];
  /* |h * p[6]| < 2^-7 * 2^-12 = 2^-19 thus the rounding error on h * p[6] is
     bounded by 2^-72. When adding to p[5] we stay in the same binade than p[5]
     thus error(yl) < ulp(p[5]) + 2^-72 = 2^-62+2^-72 < 2^-61.99 */
  double yh = p[3] + h * p[4];
  /* |h * p[4]| < 2^-7 * 2^-6 = 2^-13 thus the rounding error on h * p[4] is
     bounded by 2^-66. When adding to p[3] we stay in the same binade than p[3]
     thus error(yh) < ulp(p[3]) + 2^-66 = 2^-57+2^-66 < 2^-56.99 */
  yh = yh + hh * yl;
  /* the error on yh is bounded by:
     2^-56.99 : error on the previous value of yh
     2^-67 * 2^-9 = 2^-76 : error on hh times maximal value of yl
     2^-14 * 2^-61.99 = 2^-75.99 : maximal value of hh tomes error on yl
     2^-76 : rounding error on hh * yl since hh * yl < 2^-23
     2^-57 : rounding error on yh since we stay in the same binade as p[3]
     Total: error(yh) < 2^-55.99 */
  yh = p[2] + yh * h;
  /* |yh * h| < 2^-4 * 2^-7 = 2^-11 thus the rounding error on yh * h is
     bounded by 2^-64. After adding yh * h we stay in the same binade, thus
     the error on yh is bounded by ulp(p2) + 2^-64 = 2^-55+2^-64 < 2^-54.99.
     We need to add the error on yh multiplied by h, which is bounded by
     2^-55.99 * 2^-7 < 2^-62.99, which gives 2^-55+2^-64+2^-62.99 < 2^-54.99.
     This error is multiplied by h^2 < 2^-14, thus contributes to at most
     2^-68.99 in the final error. */
  /* add p[1] + p1l + yh * h */
  fast_two_sum (&yh, &yl, p[1], yh * h); /* exact */
  /* |yh * h| < 2^-2 * 2^-7 = 2^-9 thus the rounding error on yh * h is
     bounded by 2^-62. This rounding error is multiplied by h < 2^-7 thus
     contributes to < 2^-69 to the final error. */
  yl += p1l;
  /* |yl| < 2^-53 and |p1l| < 2^-55 thus the rounding error in yl += p1l
     is bounded by 2^-105 (we might have an exponent jump). This error is
     multiplied by h below, thus contributes < 2^-112. */
  /* multiply (yh,yl) by h */
  dekker (&yh, &t, yh, h); /* exact */
  /* add yl*h */
  t += yl * h;
  /* |yl| < 2^-52 here, thus |yl * h| < 2^-59, thus the rounding error on
     yl * h is < 2^-112.
     |yh| < 1 before the dekker() call, thus after we have |yh| < 2^-7
     and |t| < 2^-60, thus |t+yl * h| < 2^-58, and the rounding error on
     t += yl * h is < 2^-111. */
  /* add p[0] = 1 */
  fast_two_sum (&yh, &yl, p[0], yh); /* exact */
  yl += t;
  /* now |yh| < 2 and |yl| < 2^-52, with |t| < 2^-58, thus |yl+t| < 2^-51
     and the rounding error in yl += t is bounded by 2^-104. */
  /* now (yh,yl) approximates 2^h to about 68 bits of accuracy:
     2^-68.99 from the rounding errors for evaluating p[2] + ...
     2^-69 from the rounding error in yh * h in the 1st fast_two_sum
     2^-112 from the rounding error in yl += p1l
     2^-112 from the rounding error in yl * h
     2^-111 from the rounding error in t += yl * h
     2^-104 from the rounding error in yl += t
     Total absolute error < 2^-67.99 on yh+yl here (with respect to 2^h).
  */

  /* FIXME: could we integrate the multiplication by 2^l above? */
  /* multiply (yh,yl) by 2^l. Since |l| < 2^-42, it suffices to multiply
     by 1 + log(2)*l to get 70-bit accuracy: the maximal error while doing
     this is 2^-86.05. */
  static const double l2 = 0x1.62e42fefa39efp-1;
  /* error on l2 < 2^-55.25 */
  t = l2 * l * yh;
  /* we have |l2| < 0.70, |l| < 2^-42, |yh| < 1.01 thus |l2 * l * yh| < 2^-42
     and the total rounding error is bounded by 3*2^-95 */
  fast_two_sum (&yh, &u, yh, t); /* exact */
  u += yl;
  /* |yh| < 2, |u| < 2^-52, |yl| < 2^-52 thus |u+yl| < 2^-51 and the rounding
     error in u += yl is bounded by 2^-104. */
  /* now (yh,u) approximates 2^(h+l) to about 68 bits of accuracy:
     2^-67.99 from the approximation (yh,yl) for 2^h multiplied by 1+2^-42
     1.01*2^-86.05 < 2^-86.03 for the error from l2
     3*2^-95 < 2^-93.41 for the rounding error in l2 * l * yh
     2^-104 for the rounding error in u += yl
     Total absolute error < 2^-67.98 on yh+u here with respect to 2^(h+l).
  */

  /* multiply (yh,u) by 2^(i/128) */
  /* the maximal error |2^(i/128) - tab_i[127+i][0] - tab_i[127+i][1]|
     is 1/2*max(ulp(tab_i[127+i][1])) = 2^-107.
     Since we multiply by |2^(h+l)| < 1.006 this yields 2^-106.99. */
  t = yh * tab_i[127+i][1];
  /* |yh| < 1.006 and |tab_i[127+i][1]| < 0x1.fc6f89bd4f6bap-54 thus
     |yh * tab_i[127+i][1]| < 2^-53 and the rounding error on t is
     bounded by 2^-106 */
  dekker (&yh, &yl, yh, tab_i[127+i][0]); /* exact */
  /* now |yh| < 2 thus |yl| < 2^-52, |t| < 2^-53 */
  yl += t + u * tab_i[127+i][0];
  /* |u| < 2^-51 and |tab_i[127+i][0]| < 2 thus |u * tab_i[127+i][0]| < 2^-50
     and the rounding error on u * tab_i[127+i][0] is bounded by 2^-103.
     |t + u * tab_i[127+i][0]| < 2^-49 thus the rounding error when adding t
     is bounded by 2^-102.
     |yl + t + u * tab_i[127+i][0]| < 2^-48 thus the rounding error in
     yl += ,,, is bounded by 2^-101.
     Total error bounded by:
     2^-66.98 (2^-67.98 multiplied by 2^(i/128))
     2^-106.99 for the error on 2^(i/128)
     2^-106 for the rounding error on t
     2^-103 for the rounding error on u * tab_i[127+i][0]
     2^-102 for the rounding error on t + u * tab_i[127+i][0]
     2^-101 for the rounding error on yl += ...
     Total error < 2^-66.97 on yh + yl with respect to 2^(i/128+h+l). */

  /* now yh+yl approximates 2^(i/128+h+l) with error < 2^-66.97 */

  /* rounding test */
  double err = 0x1.05610cf8bb59ap-67; /* e = up(2^-66.97) */
  v.x = yh + (yl - err);
  double right = yh + (yl + err);
  if (v.x != right)
    return cr_exp_accurate (x, e, i);

  /* Multiply by 2^e. */
  unsigned int f = v.n >> 52; /* sign is always positive */
  f += e;
  if (__builtin_expect((f - 1) > 0x7ff, 0)) /* overflow or subnormal */
  {
    if (e > 0)
      return ldexp (v.x, e);
    /* In the subnormal case, we use the method mentioned by Godunov in his
       IEEE Transactions on Computers article (2020): add 2^-1022, and then
       subtract 2^-1022. Since we should multiply yh+yl by 2^e, to avoid
       underflow before adding/subtracting 2^-1022, we do it directly on
       yh+yl, where 2^-1022 is replaced by 2^(-1022-e). */
    double magic = ldexp (1.0, -1022 - e), left;
    fast_two_sum (&left, &u, magic, yh);
    left += u + (yl - err); /* here comes the rounding */
    left -= magic;
    fast_two_sum (&right, &u, magic, yh);
    right += u + (yl + err); /* here comes the rounding */
    right -= magic;
    if (left == right)
      return ldexp (left, e);
    return cr_exp_accurate (x, e, i);
  }
  v.n += (int64_t) e << 52;
  return v.x;
}
