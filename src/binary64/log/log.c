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

/* for 181 <= i <= 362, r[i] = _INVERSE[i-181] is a 9-bit approximation of
   1/x[i], where i*2^-8 <= x[i] < (i+1)*2^-8.
   More precisely r[i] is a 9-bit value such that r[i]*y-1 is representable
   exactly on 53 bits for for any y, i*2^-8 <= y < (i+1)*2^-8.
   Moreover |r[i]*y-1| <= 0.0040283203125. */
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
    0x1.05p+0, 0x1.04p+0, 0x1.03p+0, 0x1.02p+0, 0x1.01p+0, 0x1.ffp-1, 0x1.fdp-1,
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
   approximation of -log(_INVERSE[i-181]) */
static const double _LOG_INV[182][2] = {
    {-0x1.5ff3070a793d4p-2, 0x1.bc60efafc6f6ep-57},
    {-0x1.5a42ab0f4cfe2p-2, 0x1.8ebcb7dee9a3dp-56},
    {-0x1.548a2c3add263p-2, 0x1.819cf7e308ddbp-57},
    {-0x1.4ec973260026ap-2, 0x1.42a87d977dc5ep-56},
    {-0x1.49006804009d1p-2, 0x1.9ffc341f177dcp-57},
    {-0x1.432ef2a04e814p-2, 0x1.29931715ac903p-56},
    {-0x1.404308686a7e4p-2, 0x1.0bcfb6082ce6dp-56},
    {-0x1.3a64c556945eap-2, 0x1.c68651945f97cp-57},
    {-0x1.347dd9a987d55p-2, 0x1.4dd4c580919f8p-57},
    {-0x1.2e8e2bae11d31p-2, 0x1.8f4cdb95ebdf9p-56},
    {-0x1.2895a13de86a3p-2, -0x1.7ad24c13f040ep-56},
    {-0x1.22941fbcf7966p-2, 0x1.76f5eb09628afp-56},
    {-0x1.1f8ff9e48a2f3p-2, 0x1.c9fdf9a0c4b07p-56},
    {-0x1.1980d2dd4236fp-2, -0x1.9d3d1b0e4d147p-56},
    {-0x1.136870293a8bp-2, -0x1.7b66298edd24ap-56},
    {-0x1.1058bf9ae4ad5p-2, -0x1.89fa0ab4cb31dp-58},
    {-0x1.0a324e27390e3p-2, -0x1.7dcfde8061c03p-56},
    {-0x1.0402594b4d041p-2, 0x1.28ec217a5022dp-57},
    {-0x1.fb9186d5e3e2bp-3, 0x1.caaae64f21acbp-57},
    {-0x1.f550a564b7b37p-3, -0x1.c5f6dfd018c37p-61},
    {-0x1.e8c0252aa5a6p-3, 0x1.6e03a39bfc89bp-59},
    {-0x1.e27076e2af2e6p-3, 0x1.61578001e0162p-59},
    {-0x1.d5c216b4fbb91p-3, -0x1.6e443597e4d4p-57},
    {-0x1.c8ff7c79a9a22p-3, 0x1.4f689f8434012p-57},
    {-0x1.c2968558c18c1p-3, 0x1.73dee38a3fb6bp-57},
    {-0x1.b5b519e8fb5a4p-3, -0x1.ba27fdc19e1ap-57},
    {-0x1.af3c94e80bff3p-3, 0x1.398cff3641985p-58},
    {-0x1.a23bc1fe2b563p-3, -0x1.93711b07a998cp-59},
    {-0x1.9bb362e7dfb83p-3, -0x1.575e31f003e0cp-57},
    {-0x1.8e928de886d41p-3, 0x1.569d851a5677p-57},
    {-0x1.87fa06520c911p-3, 0x1.bf7fdbfa08d9ap-57},
    {-0x1.7ab890210d909p-3, -0x1.be36b2d6a0608p-59},
    {-0x1.740f8f54037a5p-3, 0x1.b264062a84cdbp-58},
    {-0x1.6d60fe719d21dp-3, 0x1.caae268ecd179p-57},
    {-0x1.5ff3070a793d4p-3, 0x1.bc60efafc6f6ep-58},
    {-0x1.59338d9982086p-3, 0x1.65d22aa8ad7cfp-58},
    {-0x1.4ba36f39a55e5p-3, -0x1.68981bcc36756p-57},
    {-0x1.44d2b6ccb7d1ep-3, -0x1.9f4f6543e1f88p-57},
    {-0x1.3dfc2b0ecc62ap-3, 0x1.ab3a8e7d81017p-58},
    {-0x1.303d718e47fd3p-3, 0x1.6b9c7d96091fap-63},
    {-0x1.29552f81ff523p-3, -0x1.301771c407dbfp-57},
    {-0x1.2266f190a5acbp-3, -0x1.f547bf1809e88p-57},
    {-0x1.14785846742acp-3, -0x1.a28813e3a7f07p-57},
    {-0x1.0d77e7cd08e59p-3, -0x1.9a5dc5e9030acp-57},
    {-0x1.0671512ca596ep-3, -0x1.50c647eb86499p-58},
    {-0x1.f0a30c01162a6p-4, -0x1.85f325c5bbacdp-58},
    {-0x1.e27076e2af2e6p-4, 0x1.61578001e0162p-60},
    {-0x1.d4313d66cb35dp-4, -0x1.790dd951d90fap-58},
    {-0x1.c5e548f5bc743p-4, -0x1.5d617ef8161b1p-60},
    {-0x1.a926d3a4ad563p-4, -0x1.942f48aa70ea9p-58},
    {-0x1.9ab42462033adp-4, 0x1.2099e1c184e8ep-59},
    {-0x1.8c345d6319b21p-4, 0x1.4a697ab3424a9p-61},
    {-0x1.7da766d7b12cdp-4, 0x1.eeedfcdd94131p-58},
    {-0x1.60658a93750c4p-4, 0x1.388458ec21b6ap-58},
    {-0x1.51b073f06183fp-4, -0x1.a49e39a1a8be4p-58},
    {-0x1.42edcbea646fp-4, -0x1.ddd4f935996c9p-59},
    {-0x1.341d7961bd1d1p-4, 0x1.b599f227becbbp-58},
    {-0x1.253f62f0a1417p-4, 0x1.c125963fc4cfdp-62},
    {-0x1.16536eea37ae1p-4, 0x1.79da3e8c22cdap-60},
    {-0x1.f0a30c01162a6p-5, -0x1.85f325c5bbacdp-59},
    {-0x1.d276b8adb0b52p-5, -0x1.1e3c53257fd47p-61},
    {-0x1.b42dd711971bfp-5, 0x1.eb9759c130499p-60},
    {-0x1.95c830ec8e3ebp-5, -0x1.f5a0e80520bf2p-59},
    {-0x1.77458f632dcfcp-5, -0x1.18d3ca87b9296p-59},
    {-0x1.58a5bafc8e4d5p-5, 0x1.ce55c2b4e2b72p-59},
    {-0x1.39e87b9febd6p-5, 0x1.5bfa937f551bbp-59},
    {-0x1.1b0d98923d98p-5, 0x1.e9ae889bac481p-60},
    {-0x1.f829b0e7833p-6, -0x1.33e3f04f1ef23p-60},
    {-0x1.b9fc027af9198p-6, 0x1.0ae69229dc868p-64},
    {-0x1.7b91b07d5b11bp-6, 0x1.5b602ace3a51p-60},
    {-0x1.3cea44346a575p-6, 0x1.0cb5a902b3a1cp-62},
    {-0x1.fc0a8b0fc03e4p-7, 0x1.83092c59642a1p-62},
    {-0x1.7dc475f810a77p-7, 0x1.16d7687d3df21p-62},
    {-0x1.fe02a6b106789p-8, 0x1.e44b7e3711ebfp-67},
    {-0x1.ff00aa2b10bcp-9, -0x1.2821ad5a6d353p-63},
    {0x1.0040155d5889ep-9, -0x1.8f98e1113f403p-65},
    {0x1.8121214586b54p-8, 0x1.c14b9f9377a1dp-65},
    {0x1.41929f96832fp-7, -0x1.c5517f64bc223p-61},
    {0x1.c317384c75f06p-7, 0x1.806208c04c22p-61},
    {0x1.228fb1fea2e28p-6, -0x1.cd7b66e01c26dp-61},
    {0x1.63d6178690bd6p-6, -0x1.8ed4d357c9c97p-64},
    {0x1.a55f548c5c43fp-6, 0x1.ec1a5f86d41f9p-62},
    {0x1.e72bf2813ce51p-6, 0x1.75b44595cab18p-60},
    {0x1.0415d89e74444p-5, 0x1.c05cf1d753622p-59},
    {0x1.252f32f8d183fp-5, -0x1.947f792615916p-59},
    {0x1.466aed42de3eap-5, -0x1.cdd6f7f4a137ep-59},
    {0x1.67c94f2d4bb58p-5, 0x1.0413e6505e603p-59},
    {0x1.894aa149fb343p-5, 0x1.a8be97660a23dp-60},
    {0x1.aaef2d0fb10fcp-5, 0x1.a353bb42e0addp-61},
    {0x1.bbcebfc68f42p-5, 0x1.e5cf3a0f56f72p-60},
    {0x1.dda8adc67ee4ep-5, 0x1.4e6c986f44c55p-59},
    {0x1.ffa6911ab9301p-5, -0x1.cd9f1f95c2eedp-59},
    {0x1.10e45b3cae831p-4, -0x1.a4a128d192686p-58},
    {0x1.2207b5c78549ep-4, -0x1.cc0fbce104eaap-58},
    {0x1.2aa04a44717a5p-4, -0x1.d15d38d2fa3f7p-58},
    {0x1.3bdf5a7d1ee64p-4, 0x1.7a976d3b5b45fp-59},
    {0x1.4d3115d207eacp-4, 0x1.769f42c7842ccp-58},
    {0x1.55e10050e0384p-4, -0x1.45f9d61c68c1bp-58},
    {0x1.674f089365a7ap-4, -0x1.9acd8b33f8fdcp-58},
    {0x1.78d02263d82d3p-4, 0x1.abca5b4fdb88p-58},
    {0x1.8197e2f40e3fp-4, 0x1.b9f2dffbeed43p-60},
    {0x1.9335e5d594989p-4, -0x1.478a85704ccb7p-58},
    {0x1.a4e7640b1bc38p-4, -0x1.5b5ca203e4259p-58},
    {0x1.adc77ee5aea8cp-4, 0x1.37d8f39bee659p-58},
    {0x1.bf968769fca11p-4, -0x1.cdc9f6f5f38c7p-59},
    {0x1.d179788219364p-4, 0x1.9daf7df76ad2ap-59},
    {0x1.da727638446a2p-4, 0x1.401fa71733019p-58},
    {0x1.ec739830a112p-4, -0x1.a2bf991780d3fp-59},
    {0x1.f57bc7d9005dbp-4, -0x1.9361574fb24e2p-58},
    {0x1.03cdc0a51ec0dp-3, 0x1.39e2d3f8b7d1p-57},
    {0x1.08598b59e3a07p-3, -0x1.dd7009902bf32p-57},
    {0x1.1178e8227e47cp-3, -0x1.0e63a5f01c691p-58},
    {0x1.160c8024b27b1p-3, -0x1.2d56ff61c2bfbp-57},
    {0x1.1f3b925f25d41p-3, 0x1.62c9ef939ac5dp-59},
    {0x1.23d712a49c202p-3, -0x1.6e38161051d69p-57},
    {0x1.2d1610c86813ap-3, -0x1.499a3f25af95fp-58},
    {0x1.31b994d3a4f85p-3, -0x1.c4716bdfc0cc9p-58},
    {0x1.3b08b6757f2a9p-3, 0x1.70d6cdf05266cp-60},
    {0x1.3fb45a59928ccp-3, -0x1.d87e6a354d056p-57},
    {0x1.4913d8333b561p-3, -0x1.0d5604930f135p-58},
    {0x1.4dc7b897bc1c8p-3, -0x1.927d47803c5f4p-57},
    {0x1.5737cc9018cddp-3, 0x1.4f4d710fec38ep-57},
    {0x1.5bf406b543db2p-3, -0x1.1f5b44c0df7e7p-61},
    {0x1.6574ebe8c133ap-3, -0x1.d34f0f4621bedp-60},
    {0x1.6a399dabbd383p-3, 0x1.96332bd4b341fp-57},
    {0x1.6f0128b756abcp-3, -0x1.8de59c21e166cp-57},
    {0x1.7898d85444c73p-3, 0x1.ef8f6ebcfb201p-58},
    {0x1.7d6903caf5adp-3, -0x1.ac5f0c075b847p-59},
    {0x1.871213750e994p-3, 0x1.d685f35eea2ap-57},
    {0x1.8beafeb38fe8cp-3, 0x1.55aa8b6997a4p-58},
    {0x1.90c6db9fcbcd9p-3, 0x1.054473941ad99p-57},
    {0x1.9a8778debaa38p-3, 0x1.f47dfd871f87fp-57},
    {0x1.9f6c407089664p-3, 0x1.35a19605e67efp-59},
    {0x1.a454082e6ab05p-3, 0x1.df207dc5c34c6p-58},
    {0x1.ae2ca6f672bd4p-3, 0x1.ab5ca9eaa088ap-57},
    {0x1.b31d8575bce3dp-3, -0x1.6353ab386a94dp-57},
    {0x1.b811730b823d2p-3, 0x1.a0ee735d9f0ecp-60},
    {0x1.bd087383bd8adp-3, 0x1.dd355f6a516d7p-60},
    {0x1.c6ffbc6f00f71p-3, -0x1.8e58b2c57a4a5p-57},
    {0x1.cc000c9db3c52p-3, 0x1.53d154280394fp-57},
    {0x1.d1037f2655e7bp-3, 0x1.60629242471a2p-57},
    {0x1.db13db0d4894p-3, 0x1.aa11d49f96cb9p-58},
    {0x1.e020cc6235ab5p-3, 0x1.fea48dd7b81d1p-58},
    {0x1.e530effe71012p-3, 0x1.2276041f43042p-59},
    {0x1.ea4449f04aaf5p-3, -0x1.d33919ab94074p-57},
    {0x1.f474b134df229p-3, -0x1.27c77ded76aadp-58},
    {0x1.f991c6cb3b379p-3, 0x1.f665066f980a2p-57},
    {0x1.feb2233ea07cdp-3, 0x1.8de00938b4c4p-61},
    {0x1.01eae5626c691p-2, -0x1.18290bd2932e2p-59},
    {0x1.047e60cde83b8p-2, -0x1.0779634061cbcp-56},
    {0x1.09aa572e6c6d4p-2, 0x1.43c2e68684d53p-57},
    {0x1.0c42d676162e3p-2, 0x1.162c79d5d11eep-58},
    {0x1.0edd060b78081p-2, -0x1.92b49ef282b09p-57},
    {0x1.1178e8227e47cp-2, -0x1.0e63a5f01c691p-57},
    {0x1.14167ef367783p-2, 0x1.e0936abd4fa6ep-62},
    {0x1.16b5ccbacfb73p-2, 0x1.66fbd28b40935p-56},
    {0x1.1bf99635a6b95p-2, -0x1.12aeb84249223p-57},
    {0x1.1e9e1678899f4p-2, 0x1.512c3749a1e4ep-56},
    {0x1.214456d0eb8d4p-2, 0x1.f7ae91aeba60ap-57},
    {0x1.23ec5991eba49p-2, 0x1.bb75d1addf87p-60},
    {0x1.269621134db92p-2, 0x1.e0efadd9db02bp-56},
    {0x1.2941afb186b7cp-2, -0x1.856e61c51574p-57},
    {0x1.2bef07cdc9354p-2, -0x1.82dad7fd86088p-56},
    {0x1.314f1e1d35ce4p-2, -0x1.3d69909e5c3dcp-56},
    {0x1.3401e12aecba1p-2, -0x1.cd55b8a4746cp-58},
    {0x1.36b6776be1117p-2, -0x1.324f0e883858ep-58},
    {0x1.396ce359bbf54p-2, -0x1.ce2b31b31e8bp-58},
    {0x1.3c25277333184p-2, -0x1.2ad27e50a8ec6p-56},
    {0x1.3edf463c1683ep-2, 0x1.83d680d3c1084p-56},
    {0x1.419b423d5e8c7p-2, 0x1.0dbb243827392p-57},
    {0x1.44591e0539f49p-2, -0x1.2b125247b0fa5p-56},
    {0x1.4718dc271c41bp-2, 0x1.8fb4c14c56eefp-60},
    {0x1.49da7f3bcc41fp-2, -0x1.9964a168ccacap-57},
    {0x1.4c9e09e172c3cp-2, -0x1.123615b147a5dp-58},
    {0x1.4f637ebba981p-2, -0x1.58cb3124b9245p-56},
    {0x1.522ae0738a3d8p-2, -0x1.8f7e9b38a6979p-57},
    {0x1.54f431b7be1a9p-2, -0x1.aacfdbbdab914p-56},
    {0x1.57bf753c8d1fbp-2, -0x1.0908d15f88b63p-57},
    {0x1.5a8cadbbedfa1p-2, -0x1.e6c2bdfb3e037p-58},
    {0x1.5d5bddf595f3p-2, -0x1.6541148cbb8a2p-56},
    {0x1.602d08af091ecp-2, -0x1.6e8920c09b73fp-58},
    {0x1.630030b3aac49p-2, 0x1.dc18ce51fff99p-57},
};

/* The following is a degree-7 polynomial generated by Sollya over
   [-0.0040283203125,0.0040283203125],
   with relative error < 2^-73.148.
   The polynomial is P[0]*x + P[1]*x^2 + ... + P[6]*x^7. */
static const double P[] = {0x1p0,                 /* degree 1 */
                           -0x1.0000000000001p-1, /* degree 2 */
                           0x1.5555555555557p-2,  /* degree 3 */
                           -0x1.fffffffea76acp-3, /* degree 4 */
                           0x1.9999999870db3p-3,  /* degree 5 */
                           -0x1.55576f0ef6485p-3, /* degree 6 */
                           0x1.2494212200e1bp-3,  /* degree 7 */
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
  static const uint64_t cm[] = {44, 45};

  int i = m >> cm[c];
  // if (v.f == TRACEM) printf ("i=%d\n", i);
  double y = v.f * cy[c];
  // if (v.f == TRACEM) printf ("y=%la\n", y);
#define OFFSET 181
  double r = (_INVERSE - OFFSET)[i];
  double l1 = (_LOG_INV - OFFSET)[i][0];
  double l2 = (_LOG_INV - OFFSET)[i][1];
  // if (v.f == TRACEM) printf ("r=%la\n", r);
  double z = __builtin_fma (r, y, -1.0); /* exact */
  // if (v.f == TRACEM) printf ("z=%la\n", z);
  /* evaluate P(z) */
  double ph, pl, z2 = z * z;
  double p56 = __builtin_fma (P[6], z, P[5]);
  ph = __builtin_fma (P[7], z2, p56);
  double p34 = __builtin_fma (P[4], z, P[3]);
  ph = __builtin_fma (ph, z2, p34);
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
  if (x <= 0.0)
  {
    /* f(x<0) is NaN, f(+/-0) is -Inf and raises DivByZero */
    if (x < 0)
      return 0.0 / 0.0;
    else
      return 1.0 / -0.0;
  }
  /* now x > 0 */
  d64u64 v = {.f = x};
  int e = (v.u >> 52) - 0x3ff, bias = 0;
  if (e == 0x400) /* +Inf or NaN */
    return x;
  /* now 0 < x < +Inf */
  if (e == -0x3ff)
  {
    v.f *= 0x1p52;
    bias = 52;
    e = (v.u >> 52) - 0x3ff;
  }
  v.u -= (int64_t) e << 52;
  e -= bias;
  /* now x = m*2^e with 1 <= m < 2 (m = v.f) */
  double h, l;
  // if (x == TRACE) printf ("x=%la e=%d m=%la\n", x, e, v.f);
  cr_log_fast (&h, &l, &e, v);
  /* err=0x1.7fp-68 + ... fails for x=0x1.8e0c521132157p-639 (rndz) */
  static double err = 0x1.80p-68 + 0x1.04p-85;
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
