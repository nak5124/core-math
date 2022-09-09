/* Correctly rounded logarithm of binary64 values.

Copyright (c) 2022 Paul Zimmermann, Inria.

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

#define TRACE 0x1.19847c76d3bf3p-3
#define TRACEM 0x1.19847c76d3bf3p+0

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

/* for 0 <= i < 64, red1[i] is such that max(abs(m*red1[i]-2^58))
   is minimized for m[i] <= m < m[i+1], with m[i] = 2^52+i*2^46 */
static int
red1[64] = {64, 63, 62, 61, 60, 59, 58, 57, 56, 56, 55, 54, 54, 53, 52, 52, 51, 50, 50, 49, 48, 48, 47, 47, 46, 46, 45, 45, 44, 44, 43, 43, 42, 42, 42, 41, 41, 40, 40, 40, 39, 39, 38, 38, 38, 37, 37, 37, 36, 36, 36, 35, 35, 35, 35, 34, 34, 34, 33, 33, 33, 33, 32, 32};

/* for 0 <= i < 64, logred1[i] = {h,l} is an approximation of log(red1[i]/64),
   with h an integer multiple of 2^-52, and l an integer multiple of 2^-104,
   with |log(red1[i]/64) - (h+l)| <= 2^-105 */
static double
logred1[64][2] = {{0x0p3, 0x0p3}, {-0x1.020565893584p-6, -0x1.d27c8e8416e7p-56}, {-0x1.0415d89e7444p-5, -0x1.1c05cf1d7536p-55}, {-0x1.894aa149fb34p-5, -0x1.9a8be97660a2p-56}, {-0x1.08598b59e3ap-4, -0x1.a228ff66fd40cp-54}, {-0x1.4d3115d207ebp-4, 0x1.d12c17a70f7a8p-55}, {-0x1.9335e5d59499p-4, 0x1.d478a85704cccp-54}, {-0x1.da727638446ap-4, -0x1.2803f4e2e66p-55}, {-0x1.1178e8227e478p-3, -0x1.ef19c5a0fe398p-54}, {-0x1.1178e8227e478p-3, -0x1.ef19c5a0fe398p-54}, {-0x1.365fcb0159018p-3, 0x1.d057dcb48d768p-55}, {-0x1.5bf406b543dbp-3, -0x1.fb8292ecfc82p-55}, {-0x1.5bf406b543dbp-3, -0x1.fb8292ecfc82p-55}, {-0x1.823c16551a3cp-3, -0x1.bb734c63d062p-55}, {-0x1.a93ed3c8ad9ep-3, -0x1.b795f53bd2e4p-54}, {-0x1.a93ed3c8ad9ep-3, -0x1.b795f53bd2e4p-54}, {-0x1.d1037f2655e78p-3, -0x1.ac0c524848e34p-54}, {-0x1.f991c6cb3b378p-3, -0x1.7d99419be6028p-55}, {-0x1.f991c6cb3b378p-3, -0x1.7d99419be6028p-55}, {-0x1.1178e8227e47cp-2, 0x1.0e63a5f01c6ap-57}, {-0x1.269621134db94p-2, 0x1.87c41489893f4p-54}, {-0x1.269621134db94p-2, 0x1.87c41489893f4p-54}, {-0x1.3c25277333184p-2, 0x1.2ad27e50a8ecp-56}, {-0x1.3c25277333184p-2, 0x1.2ad27e50a8ecp-56}, {-0x1.522ae0738a3d8p-2, 0x1.8f7e9b38a698p-57}, {-0x1.522ae0738a3d8p-2, 0x1.8f7e9b38a698p-57}, {-0x1.68ac83e9c6a14p-2, -0x1.a64eadd74018p-58}, {-0x1.68ac83e9c6a14p-2, -0x1.a64eadd74018p-58}, {-0x1.7fafa3bd8151cp-2, 0x1.219024acd3b8p-58}, {-0x1.7fafa3bd8151cp-2, 0x1.219024acd3b8p-58}, {-0x1.973a3431356acp-2, -0x1.cec5afd260f9p-54}, {-0x1.973a3431356acp-2, -0x1.cec5afd260f9p-54}, {-0x1.af5295248cddp-2, -0x1.9d56c45dd3e8p-56}, {-0x1.af5295248cddp-2, -0x1.9d56c45dd3e8p-56}, {-0x1.af5295248cddp-2, -0x1.9d56c45dd3e8p-56}, {-0x1.c7ff9c74554c8p-2, -0x1.2447d5b6ca368p-54}, {-0x1.c7ff9c74554c8p-2, -0x1.2447d5b6ca368p-54}, {-0x1.e148a1a2726ccp-2, -0x1.94df8cdd6c81p-54}, {-0x1.e148a1a2726ccp-2, -0x1.94df8cdd6c81p-54}, {-0x1.e148a1a2726ccp-2, -0x1.94df8cdd6c81p-54}, {-0x1.fb358af7a4884p-2, -0x1.7e8f05924d26p-57}, {-0x1.fb358af7a4884p-2, -0x1.7e8f05924d26p-57}, {-0x1.0ae76e2d054fap-1, -0x1.0d710fcfc4e1p-55}, {-0x1.0ae76e2d054fap-1, -0x1.0d710fcfc4e1p-55}, {-0x1.0ae76e2d054fap-1, -0x1.0d710fcfc4e1p-55}, {-0x1.188ee40f23ca6p-1, -0x1.89df1568ca0bp-55}, {-0x1.188ee40f23ca6p-1, -0x1.89df1568ca0bp-55}, {-0x1.188ee40f23ca6p-1, -0x1.89df1568ca0bp-55}, {-0x1.269621134db92p-1, -0x1.e0efadd9db028p-55}, {-0x1.269621134db92p-1, -0x1.e0efadd9db028p-55}, {-0x1.269621134db92p-1, -0x1.e0efadd9db028p-55}, {-0x1.35028ad9d8c86p-1, 0x1.f01ab6065516p-56}, {-0x1.35028ad9d8c86p-1, 0x1.f01ab6065516p-56}, {-0x1.35028ad9d8c86p-1, 0x1.f01ab6065516p-56}, {-0x1.35028ad9d8c86p-1, 0x1.f01ab6065516p-56}, {-0x1.43d9ff2f923c4p-1, -0x1.9ec2dfbeb8238p-54}, {-0x1.43d9ff2f923c4p-1, -0x1.9ec2dfbeb8238p-54}, {-0x1.43d9ff2f923c4p-1, -0x1.9ec2dfbeb8238p-54}, {-0x1.5322e26867858p-1, 0x1.99dd16d4567acp-54}, {-0x1.5322e26867858p-1, 0x1.99dd16d4567acp-54}, {-0x1.5322e26867858p-1, 0x1.99dd16d4567acp-54}, {-0x1.5322e26867858p-1, 0x1.99dd16d4567acp-54}, {-0x1.62e42fefa39fp-1, 0x1.950d871319ffp-54}, {-0x1.62e42fefa39fp-1, 0x1.950d871319ffp-54}};

/* for 0 <= j < 140, red2[j] is such that max(abs(m*red2[j]-2^75))
   is minimized for m[j] <= m < m[j+1], with m[j] = (4025+j)*2^46 */
static int
red2[140] = {133368, 133334, 133301, 133268, 133235, 133202, 133169, 133136, 133103, 133070, 133037, 133004, 132971, 132938, 132905, 132872, 132840, 132807, 132774, 132741, 132708, 132675, 132643, 132610, 132577, 132544, 132512, 132479, 132446, 132414, 132381, 132348, 132316, 132283, 132251, 132218, 132185, 132153, 132120, 132088, 132055, 132023, 131990, 131958, 131926, 131893, 131861, 131828, 131796, 131764, 131731, 131699, 131667, 131634, 131602, 131570, 131538, 131505, 131473, 131441, 131409, 131377, 131345, 131312, 131280, 131248, 131216, 131184, 131152, 131120, 131088, 131056, 131024, 130992, 130960, 130928, 130896, 130864, 130832, 130801, 130769, 130737, 130705, 130673, 130641, 130610, 130578, 130546, 130514, 130483, 130451, 130419, 130388, 130356, 130324, 130293, 130261, 130229, 130198, 130166, 130135, 130103, 130072, 130040, 130009, 129977, 129946, 129914, 129883, 129851, 129820, 129789, 129757, 129726, 129695, 129663, 129632, 129601, 129569, 129538, 129507, 129476, 129444, 129413, 129382, 129351, 129320, 129289, 129257, 129226, 129195, 129164, 129133, 129102, 129071, 129040, 129009, 128978, 128947, 128916};

/* for 0 <= j < 140, logred2[j] = {h,l} is an approximation of log(red2[j]/2^17),
   with h an integer multiple of 2^-52, and l an integer multiple of 2^-104,
   with |log(red2[j]/2^17) - (h+l)| <= 2^-105 */
static double
logred2[140][2] = {{0x1.1c83e8e4fffcp-6, 0x1.4a649b2148174p-54}, {0x1.185681008e7p-6, 0x1.fe9244a51b118p-54}, {0x1.14484a56728cp-6, -0x1.dab74b011438cp-54}, {0x1.1039d1de2dcp-6, 0x1.94d9b74a862p-57}, {0x1.0c2b178f6804p-6, 0x1.b6aac03b232p-54}, {0x1.081c1b61c7b4p-6, 0x1.e277680a193ccp-54}, {0x1.040cdd4cf194p-6, 0x1.20c57c5f1055p-54}, {0x1.fffaba9111ap-7, 0x1.5dc3c54caa378p-55}, {0x1.f7db36985ep-7, -0x1.caba49985fd3p-54}, {0x1.efbb2e9f0838p-7, -0x1.3fcdcc681b57p-55}, {0x1.e79aa2944d18p-7, 0x1.77952b71f1378p-54}, {0x1.df799267664p-7, 0x1.f08a2a17ad28p-55}, {0x1.d757fe078a1p-7, 0x1.2b26263c73294p-54}, {0x1.cf35e563ebcp-7, -0x1.09644373368cp-56}, {0x1.c713486bbb48p-7, 0x1.66153805d395p-56}, {0x1.bef0270e2578p-7, 0x1.e45b1aeff8b4p-57}, {0x1.b70ba73ae378p-7, -0x1.593fa5a75fb78p-54}, {0x1.aee780e453dp-7, -0x1.688660f0e78fp-55}, {0x1.a6c2d5f654c8p-7, -0x1.7b04632f02e48p-54}, {0x1.9e9da660066p-7, 0x1.256868de0b6ep-56}, {0x1.9677f210857p-7, 0x1.baec64f23bd4p-58}, {0x1.8e51b8f6eb88p-7, 0x1.437f3543d1c68p-55}, {0x1.866a39055838p-7, -0x1.7948d3c578b8p-58}, {0x1.7e42fa2c323p-7, 0x1.92d07be612428p-54}, {0x1.761b3656b008p-7, -0x1.d61c0ef67323p-56}, {0x1.6df2ed73de7p-7, 0x1.750f514e3481cp-54}, {0x1.66096d772e18p-7, 0x1.ecb798a90d92p-57}, {0x1.5de01e504758p-7, -0x1.ae9cde8b3aef8p-54}, {0x1.55b649e9a7c8p-7, -0x1.d6e70cdbfe7e8p-55}, {0x1.4dcb4a353828p-7, 0x1.7bc43d68ff9acp-54}, {0x1.45a06f271b8p-7, 0x1.2412cbed5a974p-54}, {0x1.3d750ea6cp-7, 0x1.f0ddf7fe707ep-57}, {0x1.35888ea912d8p-7, 0x1.f8b95b2c2da08p-55}, {0x1.2d5c271d99c8p-7, -0x1.bb37052259d6cp-54}, {0x1.256ea7ed4d68p-7, -0x1.fb81d12261b48p-55}, {0x1.1d4139148338p-7, 0x1.d3a5f5433b348p-55}, {0x1.1513447531b8p-7, -0x1.e7f58d2e9355p-54}, {0x1.0d244408f1c8p-7, 0x1.5af8c2fcbcffp-54}, {0x1.04f547b8504p-7, 0x1.bb4184fedce3cp-54}, {0x1.fa0a8ef0532p-8, -0x1.3a8bb947269bp-55}, {0x1.e9aa86678b7p-8, -0x1.a522bcf330a6p-56}, {0x1.d9c885be7a6p-8, 0x1.f467983a7fb4p-55}, {0x1.c9666cc9173p-8, 0x1.4d76c2c65a374p-54}, {0x1.b9826b761bbp-8, 0x1.d3258947c1e54p-54}, {0x1.a99d6d90388p-8, -0x1.9001f4ca70d8p-58}, {0x1.99383f10963p-8, -0x1.25e34c2c6dde4p-54}, {0x1.89513fbfb1bp-8, -0x1.92c5de88ca6e8p-55}, {0x1.78e9ff86e56p-8, 0x1.b0bce63025d14p-54}, {0x1.6900fe495edp-8, 0x1.49caeadff3dcp-58}, {0x1.5916ffd9e35p-8, -0x1.e736a595b0c9p-55}, {0x1.48aca825a6dp-8, 0x1.94e2e7a1ace2p-55}, {0x1.38c0a707bbep-8, 0x1.b88c7146a87p-55}, {0x1.28d3a85858p-8, 0x1.32ca3894d5484p-54}, {0x1.186637fe16fp-8, -0x1.9db9f7984dbcp-57}, {0x1.087735de08ap-8, -0x1.b233fdbd4279p-55}, {0x1.f10e6b998e8p-9, 0x1.99c0e19fb5698p-55}, {0x1.d12c6f5565p-9, 0x1.4576a055dc558p-55}, {0x1.b0494eb95e6p-9, 0x1.778858e930118p-54}, {0x1.9063498deep-9, -0x1.63d8865506458p-54}, {0x1.707b4780d22p-9, 0x1.e509a49a3fcbp-56}, {0x1.509148529a8p-9, 0x1.aeafa608e18p-63}, {0x1.30a54bc3cacp-9, -0x1.8fdfd3560e8f8p-54}, {0x1.10b75194da8p-9, 0x1.e1b595450c988p-55}, {0x1.df8fa31ba88p-10, 0x1.0b73d948f8a58p-54}, {0x1.9fab96dbb3cp-10, 0x1.e0bbae6e15d08p-54}, {0x1.5fc38dd9c34p-10, 0x1.63765fe8bf4p-59}, {0x1.1fd78796664p-10, 0x1.85bfdedb4492p-55}, {0x1.bfcf0724298p-11, -0x1.62ea24b2d4c1p-55}, {0x1.3fe7029a5c8p-11, 0x1.46d6550e8172p-55}, {0x1.7fee011fecp-12, -0x1.f3da8c861458p-55}, {0x1.fff8002aa8p-14, 0x1.aab110e6678bp-54}, {-0x1.0004001556p-13, 0x1.553bbb110c7fp-56}, {-0x1.8012012014p-12, -0x1.0613acbcf76ep-54}, {-0x1.4019029af9p-11, 0x1.5c807d230a938p-54}, {-0x1.c0310726818p-11, -0x1.4f241e016a61p-54}, {-0x1.202887999a8p-10, 0x1.3a6c156950404p-54}, {-0x1.603c8de0e98p-10, 0x1.ebe193def6e74p-54}, {-0x1.a05496e9a6p-10, 0x1.add98ca8dc6p-60}, {-0x1.e070a33460cp-10, 0x1.75392455b7e98p-55}, {-0x1.0f47d198becp-9, 0x1.76bb696aa587p-54}, {-0x1.2f59cbaf1bep-9, -0x1.01ce9a670e088p-54}, {-0x1.4f6dc82595cp-9, 0x1.a611b7ad40f7p-55}, {-0x1.6f83c73ca46p-9, 0x1.9d25ca9757d2p-57}, {-0x1.8f9bc934cc2p-9, -0x1.951185341598cp-54}, {-0x1.afb5ce4e9dap-9, 0x1.bc29c49db91ep-56}, {-0x1.ced0eeb9136p-9, 0x1.ba47b5ae63158p-55}, {-0x1.eeeeeaba07ep-9, -0x1.382db953a78ep-54}, {-0x1.0787754e4d5p-8, 0x1.bc6646de7677p-56}, {-0x1.17987750c62p-8, 0x1.d7058b9158af8p-54}, {-0x1.2729e77a61cp-8, -0x1.1a03dcfde50a8p-55}, {-0x1.373ce5ed653p-8, 0x1.295852b0c06ccp-54}, {-0x1.4750e6d1da1p-8, -0x1.d9f78e1daa68p-55}, {-0x1.56e53e423e1p-8, -0x1.62bf1b529d6a8p-54}, {-0x1.66fb3c54f4bp-8, 0x1.4f628884ad56p-55}, {-0x1.77123d3994bp-8, -0x1.f0975c302825cp-54}, {-0x1.86a97d06321p-8, 0x1.1bb14bc402cp-58}, {-0x1.96c27bd7352p-8, 0x1.7da31f4ab1f08p-55}, {-0x1.a6dc7ddacfap-8, 0x1.3310362bb74cp-55}, {-0x1.b676a719b2cp-8, 0x1.3498c30071788p-55}, {-0x1.c692a7c83c2p-8, -0x1.a9415ef965004p-54}, {-0x1.d62ebffa027p-8, -0x1.fb5b68d37a42p-54}, {-0x1.e64cbfd2128p-8, 0x1.1138c113c0b44p-54}, {-0x1.f5eac77175cp-8, -0x1.ebad95e6b08p-63}, {-0x1.03056378e998p-7, 0x1.ad2d42334782p-54}, {-0x1.0ad55f3cdd68p-7, 0x1.5b129be549068p-55}, {-0x1.12e65f10adc8p-7, -0x1.a62370435a11p-54}, {-0x1.1ab753066e68p-7, 0x1.d46281a7af788p-54}, {-0x1.22c9532d7318p-7, -0x1.e96f3a849bf1p-55}, {-0x1.2a9b3f92a2p-7, 0x1.4e3d853bb2fb8p-54}, {-0x1.32ae404c8558p-7, 0x1.ecc3ae2d922c8p-54}, {-0x1.3a81255edbp-7, 0x1.7b32008ddee6p-54}, {-0x1.425484e53408p-7, -0x1.aa30274aca84p-56}, {-0x1.4a6904e8aafp-7, 0x1.a8a070b3b9584p-54}, {-0x1.523d5d786cd8p-7, -0x1.cb8395e889e4cp-54}, {-0x1.5a1230a9a468p-7, -0x1.d57212b481908p-54}, {-0x1.62283084eee8p-7, 0x1.d90f66bee01f8p-55}, {-0x1.69fdfd1c04d8p-7, 0x1.3f41c6cb876p-60}, {-0x1.71d444821bf8p-7, -0x1.54075429ae77cp-54}, {-0x1.79ebc4c3c958p-7, 0x1.9fdc4bd21b294p-54}, {-0x1.81c305ec67ep-7, 0x1.1bce74b43cac8p-55}, {-0x1.899ac211ac9p-7, 0x1.0bc0f8a954678p-54}, {-0x1.9172f942aa8p-7, -0x1.5dd03a2a1c744p-54}, {-0x1.998c798d4458p-7, 0x1.be533990c3a34p-54}, {-0x1.a165aafc51bp-7, 0x1.df3752c57cd5p-55}, {-0x1.a93f57a4df2p-7, -0x1.657db5a92422cp-54}, {-0x1.b1197f960b08p-7, 0x1.5a763218bb16cp-54}, {-0x1.b8f422def67p-7, -0x1.b926d2001b67p-56}, {-0x1.c0cf418ec548p-7, -0x1.dc44d2233dbd8p-55}, {-0x1.c8ebc1b3b688p-7, -0x1.593da636bbf4cp-54}, {-0x1.d0c7db5b0da8p-7, 0x1.426370b4d4eep-56}, {-0x1.d8a4709741f8p-7, 0x1.b199dab13c874p-54}, {-0x1.e08181778298p-7, -0x1.423fd87688478p-55}, {-0x1.e85f0e0b0188p-7, -0x1.732b221c01cfp-56}, {-0x1.f03d1660f38p-7, -0x1.d9045b3d7a158p-54}, {-0x1.f81b9a889018p-7, -0x1.a00bb40f2832p-56}, {-0x1.fffa9a9111a8p-7, 0x1.511e1daf00398p-54}, {-0x1.03ed0b44daacp-6, 0x1.01f6a89c181fcp-54}, {-0x1.07dd0740dd94p-6, 0x1.593c480e90c1p-56}, {-0x1.0bcd414432f4p-6, 0x1.3b04523ae5e18p-55}, {-0x1.0fbdb9567d98p-6, -0x1.8e3c214b11e84p-54}};

/* given 1 <= x < 2, put in h+l a double-double approximation of log(x),
   with absolute error less than 2^-62.93 (see details below), where x=v.f.
   We also have |l| < 2^-50. */
static void
cr_log_fast (double *h, double *l, d64u64 v)
{
  int64_t m = 0x10000000000000 + (v.u & 0xfffffffffffff);
  /* m/2^52 = x */
  int i = (v.u >> 46) & 0x3f; /* 0 <= i < 64 */
  // if (x == TRACEM) printf ("x=%la m=%ld i=%d\n", x, m, i);
  m = m * (int64_t) red1[i];
  /* now m/2^58 = x * red1[i]/2^6 */
  // if (x == TRACEM) printf	("red1[i]=%d m=%ld\n", red1[i], m);
  /* -0x11c00000000000 <= m - 2^58 <= 0x113fffffffffdd */
  int j = (m >> 46) - 4025; /* 0 <= j <= 139 */
  m = m * (int64_t) red2[j];
  // if (x == TRACEM) printf ("j=%d red2[j]=%d m=%ld\n", j, red2[j], m);
  /* -0x42a1000000000000 <= m - 2^75 <= 0x430bfffffffdf708
     (where the 2^75 disappeared mod 2^64) */
  /* now 1+m/2^75 = x * red1[i]/2^6 * red2[i]/2^17
     thus log(x) = -log(red1[i]/2^6) - log(red2[i]/2^17) + log(1+m/2^75) */
  double y = m * 0x1p-75; /* rounding error < 2^-65 since m has at most 63
                             bits, thus at most the low 10 bits are lost */
  /* We use the following degree-4 polynomial generated by Sollya over
     [-0.00012709, 0.00012789], with maximal absolute error < 2^-69.90:
     c[0]*y + c[1]*y^2 + c[2]*y^3 + c[3]*y^4 */
  // if (x == TRACEM) printf ("y=%la\n", y);
  double c[4] = {1.0, -0x1.0000000000033p-1, 0x1.5555558631a59p-2,
                 -0x1.ffffb56b03d95p-3};
  double yy = y * y; /* error < ulp(yy) <= 2^-78 */
  // if (x == TRACEM) printf ("yy=%la\n", yy);
  double c23 = c[2] + y * c[3]; /* err(y*c[3]) < 2^-67, c23 is in the same binade as c[2],
                                   thus err(c23) <= ulp(c[2]) + 2^-67 <= 2^-54 + 2^-67
                                   <= 2^-53.99 and |c23| < 2^-1.58 */
  // if (x == TRACEM) printf ("c23=%la\n", c23);
  double c01 = c[0] + y * c[1]; /* err(y*c[1]) < 2^-66, c01 might be in the [1,2) binade,
                                   thus err(c01) <= ulp(1) + 2^-66 <= 2^-52 + 2^-66
                                   <= 2^-51.99 */
  // if (x == TRACEM) printf ("c01=%la\n", c01);
  double p = c01 + yy * c23; /* |yy*c23| < 2^-27.44 thus err(yy*c23) <= 2^-80
                                p is at most in the [1,2) binade,
                                thus err(p) <= ulp(1) + 2^-80 <= 2^-52+2^-80 <= 2^-51.99 */
  p = y * p; /* now |p|<0.000128 thus rounding error bounded by ulp(0.000128) < 2^-65 */
  /* Now p approximates log(1+m/2^75), with total error bounded by:
c     - 2^-65 for the difference between P(y) and P(m/2^-75), since
       |y - m/2^75| < 2^-65 and |P(y)-P(y')| <= |y-y'| * P'(t) for t in (y,y')
     - 2^-69.90 for the Sollya approximation error
     - 2^-92.51 for the error on yy multiplied by max(c23)*ymax
     - 2^-92.78 for the error on c23 multiplied by y^3
     - 2^-64.92 for the error on c01 multiplied by y
     - 2^-64.92 for the error on p (before y*p) multiplied by y
     - 2^-65 for the final rounding error on y * p
     This gives a bound of 2^-62.94 for |p - log(1+m/2^75)| for the last value of m.
   */
  // if (x == TRACEM) printf ("p=%la\n", p);
  *h = - logred1[i][0] - logred2[j][0]; /* no rounding error since both logred1[i][0]
                                           and logred2[j][0] are integer multiples of
                                           2^-52, with multipliers at most 2^52 */
  // if (x == TRACEM) printf ("logred1[i][0]=%la\n", logred1[i][0]);
  // if (x == TRACEM) printf ("logred2[j][0]=%la\n", logred2[j][0]);
  *l = - logred1[i][1] - logred2[j][1]; /* idem: both logred1[i][1]
                                           and logred2[j][1] are integer multiples of
                                           2^-104, with multipliers at most 2^52 */

  double ll;
  fast_two_sum (h, &ll, *h, p);
  *l += ll;
  // if (x == TRACEM) printf ("h=%la l=%la\n", *h, *l);
  /* Total approximation error:
     a) error in log(red1[i]/64) < 2^-105
     b) error in log(red2[j]/2^17) < 2^-105
     c) error in p < 2^-62.94
     d) error in fast_two_sum < 2^-104 |h| <= 2^-104*log(2) <= 2^-104.52
     e) error in *l + ll: since |h| < log(2) < 1, we have |ll| < ulp(log(2)) <= 2^-53
        and |*l| <= 2*2^-104*2^52 = 2^-51, thus |*l + ll| <= 2^-50 and that error is
        less than 2^-103
     Total error < 2^-62.93 */
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
  /* now x = m*2^e with 1 <= m < 2 */
  double h, l;
  // if (x == TRACE) printf ("x=%la e=%d m=%la\n", x, e, m);
  cr_log_fast (&h, &l, v);
  double err = 0x1.0dp-63; /* 2^-62.93 < err (maximal error for cr_log_fast) */
  if (e != 0)
  {
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
       Total < 2^-84.98 */
    err += 0x1.04p-85; /* 2^-84.98 < 1.04e-85 */
  }
  // if (x == TRACE) printf ("h=%la l=%la err=%la\n", h, l, err);
  double left = h + (l - err), right = h + (l + err);
  if (left == right)
    return h + l;
  return 0;
}
