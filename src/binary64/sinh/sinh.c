/* Correctly rounded sinh for binary64 values.

Copyright (c) 2022-2023 INRIA.
Author: Paul Zimmermann.

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

/* References:
   [1] IA-64 and Elementary Functions, Peter Markstein,
       Hewlett-Packard Professional Books, 2000, Chapter 16. */

#define TRACE 0x1.003a9d3e03c2p-20

#include <stdio.h>
#include <stdint.h>

#if 0
sage: build_table_T()
static const double T[256][3] = {
   {0x0p+0, 0x0p+0, 0x1p+0}, /* i=0 */
   {0x1.633d9a9a4e353p+1, 0x1.ff678ce2fd4c2p+2, 0x1.01b26143a4463p+3}, /* i=1 */
   {0x1.633d9a9affe92p+2, 0x1.0165a66164185p+7, 0x1.0167a398018a8p+7}, /* i=2 */
   {0x1.0a6e33f544a6cp+3, 0x1.021ab2c424f31p+11, 0x1.021ab4bff8544p+11}, /* i=3 */
   {0x1.633d9a68a71d3p+3, 0x1.02cf3ed313846p+15, 0x1.02cf3ed50df58p+15}, /* i=4 */

   {0x1.633d9a9a6d16cp+4, 0x1.05a665980f2f6p+31, 0x1.05a665980f2f6p+31}, /* i=8 */
   {0x1.8fa54dedba469p+4, 0x1.065d6d3c0e58cp+35, 0x1.065d6d3c0e58cp+35}, /* i=9 */
   {0x1.bc0d014107efcp+4, 0x1.0714f4e89ef61p+39, 0x1.0714f4e89ef61p+39}, /* i=10 */
   {0x1.e874b49455ac3p+4, 0x1.07ccfcf74a6a3p+43, 0x1.07ccfcf74a6a3p+43}, /* i=11 */
   {0x1.0a6e33f3d1985p+5, 0x1.088585c1da689p+47, 0x1.088585c1da689p+47}, /* i=12 */
   {0x1.20a20d9d78bcbp+5, 0x1.093e8fa26c9fbp+51, 0x1.093e8fa26c9fbp+51}, /* i=13 */
   {0x1.36d5e7471f37p+5, 0x1.09f81af32ab2ap+55, 0x1.09f81af32ab2ap+55}, /* i=14 */
   {0x1.4d09c0f0c66ebp+5, 0x1.0ab2280ecec9bp+59, 0x1.0ab2280ecec9bp+59}, /* i=15 */
   {0x1.633d9a9a6cebdp+5, 0x1.0b6cb74ff2e05p+63, 0x1.0b6cb74ff2e05p+63}, /* i=16 */
   {0x1.79717444138c9p+5, 0x1.0c27c911be3c9p+67, 0x1.0c27c911be3c9p+67}, /* i=17 */
   {0x1.8fa54dedba627p+5, 0x1.0ce35daf7d453p+71, 0x1.0ce35daf7d453p+71}, /* i=18 */
   {0x1.a5d927976143cp+5, 0x1.0d9f7584b4788p+75, 0x1.0d9f7584b4788p+75}, /* i=19 */
   {0x1.bc0d0141080e3p+5, 0x1.0e5c10ed29475p+79, 0x1.0e5c10ed29475p+79}, /* i=20 */
   {0x1.d240daeaaed49p+5, 0x1.0f193044e8536p+83, 0x1.0f193044e8536p+83}, /* i=21 */
   {0x1.e874b49455ceap+5, 0x1.0fd6d3e8438ecp+87, 0x1.0fd6d3e8438ecp+87}, /* i=22 */
   {0x1.fea88e3dfc7bap+5, 0x1.1094fc33b5011p+91, 0x1.1094fc33b5011p+91}, /* i=23 */
   {0x1.0a6e33f3d1953p+6, 0x1.1153a98412f21p+95, 0x1.1153a98412f21p+95}, /* i=24 */
   {0x1.158820c8a4f9ep+6, 0x1.1212dc366d78bp+99, 0x1.1212dc366d78bp+99}, /* i=25 */
   {0x1.20a20d9d7865bp+6, 0x1.12d294a810ccbp+103, 0x1.12d294a810ccbp+103}, /* i=26 */
   {0x1.2bbbfa724bbdfp+6, 0x1.1392d33684e19p+107, 0x1.1392d33684e19p+107}, /* i=27 */
   {0x1.36d5e7471ec7cp+6, 0x1.1453983f8a598p+111, 0x1.1453983f8a598p+111}, /* i=28 */
   {0x1.41efd41bf2ba3p+6, 0x1.1514e4218764p+115, 0x1.1514e4218764p+115}, /* i=29 */
   {0x1.4d09c0f0c62dfp+6, 0x1.15d6b73a6f05dp+119, 0x1.15d6b73a6f05dp+119}, /* i=30 */
   {0x1.5823adc599accp+6, 0x1.169911e8fc6bdp+123, 0x1.169911e8fc6bdp+123}, /* i=31 */
   {0x1.633d9a9a6cedap+6, 0x1.175bf48bf3d6dp+127, 0x1.175bf48bf3d6dp+127}, /* i=32 */
   {0x1.6e57876f402ebp+6, 0x1.181f5f82809f5p+131, 0x1.181f5f82809f5p+131}, /* i=33 */
   {0x1.7971744413adcp+6, 0x1.18e3532c10be2p+135, 0x1.18e3532c10be2p+135}, /* i=34 */
   {0x1.848b6118e6d4cp+6, 0x1.19a7cfe81acap+139, 0x1.19a7cfe81acap+139}, /* i=35 */
   {0x1.8fa54dedba6b5p+6, 0x1.1a6cd616b83cap+143, 0x1.1a6cd616b83cap+143}, /* i=36 */
   {0x1.9abf3ac28de24p+6, 0x1.1b326617e76d8p+147, 0x1.1b326617e76d8p+147}, /* i=37 */
   {0x1.a5d9279761395p+6, 0x1.1bf8804c112bp+151, 0x1.1bf8804c112bp+151}, /* i=38 */
   {0x1.b0f3146c33f86p+6, 0x1.1cbf2513c05e9p+155, 0x1.1cbf2513c05e9p+155}, /* i=39 */
   {0x1.bc0d01410812bp+6, 0x1.1d8654d06fe7dp+159, 0x1.1d8654d06fe7dp+159}, /* i=40 */
   {0x1.c726ee15db3bap+6, 0x1.1e4e0fe2afd0cp+163, 0x1.1e4e0fe2afd0cp+163}, /* i=41 */
   {0x1.d240daeaaf18cp+6, 0x1.1f1656ac6d8dcp+167, 0x1.1f1656ac6d8dcp+167}, /* i=42 */
   {0x1.dd5ac7bf82218p+6, 0x1.1fdf298ef72c9p+171, 0x1.1fdf298ef72c9p+171}, /* i=43 */
   {0x1.e874b494558c3p+6, 0x1.20a888eca4077p+175, 0x1.20a888eca4077p+175}, /* i=44 */
   {0x1.f38ea16928e89p+6, 0x1.2172752799893p+179, 0x1.2172752799893p+179}, /* i=45 */
   {0x1.fea88e3dfc95cp+6, 0x1.223ceea27c49p+183, 0x1.223ceea27c49p+183}, /* i=46 */
   {0x1.04e13d8967e5fp+7, 0x1.2307f5bfe23fdp+187, 0x1.2307f5bfe23fdp+187}, /* i=47 */
   {0x1.0a6e33f3d1d4p+7, 0x1.23d38ae330196p+191, 0x1.23d38ae330196p+191}, /* i=48 */
   {0x1.0ffb2a5e3b3bap+7, 0x1.249fae6f42f72p+195, 0x1.249fae6f42f72p+195}, /* i=49 */
   {0x1.158820c8a5169p+7, 0x1.256c60c847dc2p+199, 0x1.256c60c847dc2p+199}, /* i=50 */
   {0x1.1b1517330ee38p+7, 0x1.2639a251d935cp+203, 0x1.2639a251d935cp+203}, /* i=51 */
   {0x1.20a20d9d784d7p+7, 0x1.2707736ff00b4p+207, 0x1.2707736ff00b4p+207}, /* i=52 */
   {0x1.262f0407e23c6p+7, 0x1.27d5d4878215dp+211, 0x1.27d5d4878215dp+211}, /* i=53 */
   {0x1.2bbbfa724be1p+7, 0x1.28a4c5fcce3fap+215, 0x1.28a4c5fcce3fap+215}, /* i=54 */
   {0x1.3148f0dcb5a82p+7, 0x1.29744835103f6p+219, 0x1.29744835103f6p+219}, /* i=55 */
   {0x1.36d5e7471f2bcp+7, 0x1.2a445b9550519p+223, 0x1.2a445b9550519p+223}, /* i=56 */
   {0x1.3c62ddb188ca2p+7, 0x1.2b1500834f862p+227, 0x1.2b1500834f862p+227}, /* i=57 */
   {0x1.41efd41bf22ap+7, 0x1.2be63764ab127p+231, 0x1.2be63764ab127p+231}, /* i=58 */
   {0x1.477cca865cafbp+7, 0x1.2cb800a04c859p+235, 0x1.2cb800a04c859p+235}, /* i=59 */
   {0x1.4d09c0f0c611p+7, 0x1.2d8a5c9b3c0dep+239, 0x1.2d8a5c9b3c0dep+239}, /* i=60 */
   {0x1.5296b75b2faadp+7, 0x1.2e5d4bbcede7ap+243, 0x1.2e5d4bbcede7ap+243}, /* i=61 */
   {0x1.5823adc599747p+7, 0x1.2f30ce6c4c9ffp+247, 0x1.2f30ce6c4c9ffp+247}, /* i=62 */
   {0x1.5db0a430042f6p+7, 0x1.3004e51102c5cp+251, 0x1.3004e51102c5cp+251}, /* i=63 */
   {0x1.633d9a9a6cfa3p+7, 0x1.30d99010da138p+255, 0x1.30d99010da138p+255}, /* i=64 */
   {0x1.68ca9104d6a26p+7, 0x1.31aecfd54485bp+259, 0x1.31aecfd54485bp+259}, /* i=65 */
   {0x1.6e57876f402eep+7, 0x1.3284a4c5beb68p+263, 0x1.3284a4c5beb68p+263}, /* i=66 */
   {0x1.73e47dd9a9916p+7, 0x1.335b0f4a9958p+267, 0x1.335b0f4a9958p+267}, /* i=67 */
   {0x1.797174441351p+7, 0x1.34320fccc7f6ap+271, 0x1.34320fccc7f6ap+271}, /* i=68 */
   {0x1.7efe6aae7d3f8p+7, 0x1.3509a6b51ab9fp+275, 0x1.3509a6b51ab9fp+275}, /* i=69 */
   {0x1.848b6118e70f7p+7, 0x1.35e1d46c9800bp+279, 0x1.35e1d46c9800bp+279}, /* i=70 */
   {0x1.8a18578350cacp+7, 0x1.36ba995cc446p+283, 0x1.36ba995cc446p+283}, /* i=71 */
   {0x1.8fa54dedbaa14p+7, 0x1.3793f5ef84edep+287, 0x1.3793f5ef84edep+287}, /* i=72 */
   {0x1.95324458241ddp+7, 0x1.386dea8ea5074p+291, 0x1.386dea8ea5074p+291}, /* i=73 */
   {0x1.9abf3ac28dbf3p+7, 0x1.394877a4ce455p+295, 0x1.394877a4ce455p+295}, /* i=74 */
   {0x1.a04c312cf7a4fp+7, 0x1.3a239d9cbb339p+299, 0x1.3a239d9cbb339p+299}, /* i=75 */
   {0x1.a5d9279761494p+7, 0x1.3aff5ce10b70ap+303, 0x1.3aff5ce10b70ap+303}, /* i=76 */
   {0x1.ab661e01cae72p+7, 0x1.3bdbb5dd1f773p+307, 0x1.3bdbb5dd1f773p+307}, /* i=77 */
   {0x1.b0f3146c34cf9p+7, 0x1.3cb8a8fcb0c6dp+311, 0x1.3cb8a8fcb0c6dp+311}, /* i=78 */
   {0x1.b6800ad69e9e1p+7, 0x1.3d9636ab540e4p+315, 0x1.3d9636ab540e4p+315}, /* i=79 */
   {0x1.bc0d014107b91p+7, 0x1.3e745f54c7bap+319, 0x1.3e745f54c7bap+319}, /* i=80 */
   {0x1.c199f7ab74e95p+7, 0x1.3f532368706d3p+323, 0x1.3f532368706d3p+323}, /* i=81 */
   {0x1.c726ee15db52bp+7, 0x1.4032834c3e623p+327, 0x1.4032834c3e623p+327}, /* i=82 */
   {0x1.ccb3e480451dcp+7, 0x1.41127f7389181p+331, 0x1.41127f7389181p+331}, /* i=83 */
   {0x1.d240daeaaea2cp+7, 0x1.41f31849566bfp+335, 0x1.41f31849566bfp+335}, /* i=84 */
   {0x1.d7cdd155188f4p+7, 0x1.42d44e3badbbfp+339, 0x1.42d44e3badbbfp+339}, /* i=85 */
   {0x1.dd5ac7bf823d6p+7, 0x1.43b621b80da76p+343, 0x1.43b621b80da76p+343}, /* i=86 */
   {0x1.e2e7be29ec09ap+7, 0x1.4498932ce457ap+347, 0x1.4498932ce457ap+347}, /* i=87 */
   {0x1.e874b49455d15p+7, 0x1.457ba3089cd5dp+351, 0x1.457ba3089cd5dp+351}, /* i=88 */
   {0x1.ee01aafec0feap+7, 0x1.465f51baec36dp+355, 0x1.465f51baec36dp+355}, /* i=89 */
   {0x1.f38ea16929805p+7, 0x1.47439fb05688ep+359, 0x1.47439fb05688ep+359}, /* i=90 */
   {0x1.f91b97d392d91p+7, 0x1.48288d5a81afp+363, 0x1.48288d5a81afp+363}, /* i=91 */
   {0x1.fea88e3dfc8e5p+7, 0x1.490e1b28d76a8p+367, 0x1.490e1b28d76a8p+367}, /* i=92 */
   {0x1.021ac25432eccp+8, 0x1.49f4498addcd8p+371, 0x1.49f4498addcd8p+371}, /* i=93 */
   {0x1.04e13d8967f41p+8, 0x1.4adb18f1ab05ap+375, 0x1.4adb18f1ab05ap+375}, /* i=94 */
   {0x1.07a7b8be9cbp+8, 0x1.4bc289cd0245ep+379, 0x1.4bc289cd0245ep+379}, /* i=95 */
   {0x1.0a6e33f3d1b4p+8, 0x1.4caa9c8e939f6p+383, 0x1.4caa9c8e939f6p+383}, /* i=96 */
   {0x1.0d34af290660bp+8, 0x1.4d9351a6d0729p+387, 0x1.4d9351a6d0729p+387}, /* i=97 */
   {0x1.0ffb2a5e3b79cp+8, 0x1.4e7ca98847006p+391, 0x1.4e7ca98847006p+391}, /* i=98 */
   {0x1.12c1a59370409p+8, 0x1.4f66a4a3dda27p+395, 0x1.4f66a4a3dda27p+395}, /* i=99 */
   {0x1.158820c8a5155p+8, 0x1.5051436c406cfp+399, 0x1.5051436c406cfp+399}, /* i=100 */
   {0x1.184e9bfdda00dp+8, 0x1.513c8653f999dp+403, 0x1.513c8653f999dp+403}, /* i=101 */
   {0x1.1b1517330ee38p+8, 0x1.52286dcdae45ep+407, 0x1.52286dcdae45ep+407}, /* i=102 */
   {0x1.1ddb92684382ep+8, 0x1.5314fa4c2fbafp+411, 0x1.5314fa4c2fbafp+411}, /* i=103 */
   {0x1.20a20d9d786a9p+8, 0x1.54022c43a6159p+415, 0x1.54022c43a6159p+415}, /* i=104 */
   {0x1.236888d2ad2e5p+8, 0x1.54f00427421a2p+419, 0x1.54f00427421a2p+419}, /* i=105 */
   {0x1.262f0407e22e3p+8, 0x1.55de826b94cb3p+423, 0x1.55de826b94cb3p+423}, /* i=106 */
   {0x1.28f57f3d17044p+8, 0x1.56cda784789fap+427, 0x1.56cda784789fap+427}, /* i=107 */
   {0x1.2bbbfa724bc43p+8, 0x1.57bd73e6bb4a7p+431, 0x1.57bd73e6bb4a7p+431}, /* i=108 */
   {0x1.2e8275a780aacp+8, 0x1.58ade807b3b41p+435, 0x1.58ade807b3b41p+435}, /* i=109 */
   {0x1.3148f0dcb4d6ep+8, 0x1.599f045b89a67p+439, 0x1.599f045b89a67p+439}, /* i=110 */
   {0x1.340f6c11ea758p+8, 0x1.5a90c95ad6554p+443, 0x1.5a90c95ad6554p+443}, /* i=111 */
   {0x1.36d5e7471f393p+8, 0x1.5b8337787951cp+447, 0x1.5b8337787951cp+447}, /* i=112 */
   {0x1.399c627c53d57p+8, 0x1.5c764f2bb2a5cp+451, 0x1.5c764f2bb2a5cp+451}, /* i=113 */
   {0x1.3c62ddb188ecbp+8, 0x1.5d6a10ec007b4p+455, 0x1.5d6a10ec007b4p+455}, /* i=114 */
   {0x1.3f2958e6be119p+8, 0x1.5e5e7d2fc261p+459, 0x1.5e5e7d2fc261p+459}, /* i=115 */
   {0x1.41efd41bf2a6ap+8, 0x1.5f53946d6784fp+463, 0x1.5f53946d6784fp+463}, /* i=116 */
   {0x1.44b64f51274e5p+8, 0x1.6049571d68584p+467, 0x1.6049571d68584p+467}, /* i=117 */
   {0x1.477cca865c082p+8, 0x1.613fc5b7b2902p+471, 0x1.613fc5b7b2902p+471}, /* i=118 */
   {0x1.4a4345bb91424p+8, 0x1.6236e0b520ec1p+475, 0x1.6236e0b520ec1p+475}, /* i=119 */
   {0x1.4d09c0f0c5fb6p+8, 0x1.632ea88ce6148p+479, 0x1.632ea88ce6148p+479}, /* i=120 */
   {0x1.4fd03c25fad6cp+8, 0x1.64271db8ce3fep+483, 0x1.64271db8ce3fep+483}, /* i=121 */
   {0x1.5296b75b2fff7p+8, 0x1.652040b255096p+487, 0x1.652040b255096p+487}, /* i=122 */
   {0x1.555d329064958p+8, 0x1.661a11f1d6083p+491, 0x1.661a11f1d6083p+491}, /* i=123 */
   {0x1.5823adc5994b4p+8, 0x1.671491f232ac3p+495, 0x1.671491f232ac3p+495}, /* i=124 */
   {0x1.5aea28facdb9fp+8, 0x1.680fc12d19165p+499, 0x1.680fc12d19165p+499}, /* i=125 */
   {0x1.5db0a430030bcp+8, 0x1.690ba01ec0cb5p+503, 0x1.690ba01ec0cb5p+503}, /* i=126 */
   {0x1.60771f6537f7fp+8, 0x1.6a082f40450afp+507, 0x1.6a082f40450afp+507}, /* i=127 */
   {0x1.633d9a9a6cdaap+8, 0x1.6b056f0d66d12p+511, 0x1.6b056f0d66d12p+511}, /* i=128 */
   {0x1.660415cfa19eep+8, 0x1.6c0360019df87p+515, 0x1.6c0360019df87p+515}, /* i=129 */
   {0x1.68ca9104d68dbp+8, 0x1.6d0202993e6b3p+519, 0x1.6d0202993e6b3p+519}, /* i=130 */
   {0x1.6b910c3a0b853p+8, 0x1.6e0157505b196p+523, 0x1.6e0157505b196p+523}, /* i=131 */
   {0x1.6e57876f4026ap+8, 0x1.6f015ea306132p+527, 0x1.6f015ea306132p+527}, /* i=132 */
   {0x1.711e02a475703p+8, 0x1.7002190f9cb1bp+531, 0x1.7002190f9cb1bp+531}, /* i=133 */
   {0x1.73e47dd9a9aa4p+8, 0x1.71038710ef972p+535, 0x1.71038710ef972p+535}, /* i=134 */
   {0x1.76aaf90ede9dbp+8, 0x1.7205a9272feb4p+539, 0x1.7205a9272feb4p+539}, /* i=135 */
   {0x1.7971744413ap+8, 0x1.73087fcf60a0dp+543, 0x1.73087fcf60a0dp+543}, /* i=136 */
   {0x1.7c37ef794859fp+8, 0x1.740c0b8753153p+547, 0x1.740c0b8753153p+547}, /* i=137 */
   {0x1.7efe6aae7d601p+8, 0x1.75104cce8744p+551, 0x1.75104cce8744p+551}, /* i=138 */
   {0x1.81c4e5e3b23cdp+8, 0x1.76154423533e2p+555, 0x1.76154423533e2p+555}, /* i=139 */
   {0x1.848b6118e6803p+8, 0x1.771af2046c331p+559, 0x1.771af2046c331p+559}, /* i=140 */
   {0x1.8751dc4e1bd19p+8, 0x1.782156f3f136dp+563, 0x1.782156f3f136dp+563}, /* i=141 */
   {0x1.8a18578350c2p+8, 0x1.7928736fd522ap+567, 0x1.7928736fd522ap+567}, /* i=142 */
   {0x1.8cded2b88594cp+8, 0x1.7a3047f8dfeaep+571, 0x1.7a3047f8dfeaep+571}, /* i=143 */
   {0x1.8fa54dedba894p+8, 0x1.7b38d5102ec3ap+575, 0x1.7b38d5102ec3ap+575}, /* i=144 */
   {0x1.926bc922edcf4p+8, 0x1.7c421b342a25dp+579, 0x1.7c421b342a25dp+579}, /* i=145 */
   {0x1.9532445823e63p+8, 0x1.7d4c1aecf3c9ap+583, 0x1.7d4c1aecf3c9ap+583}, /* i=146 */
   {0x1.97f8bf8d58e19p+8, 0x1.7e56d4b686c16p+587, 0x1.7e56d4b686c16p+587}, /* i=147 */
   {0x1.9abf3ac28e37cp+8, 0x1.7f6249153991ep+591, 0x1.7f6249153991ep+591}, /* i=148 */
   {0x1.9d85b5f7c2a2cp+8, 0x1.806e7889a82aap+595, 0x1.806e7889a82aap+595}, /* i=149 */
   {0x1.a04c312cf790cp+8, 0x1.817b6398d6b96p+599, 0x1.817b6398d6b96p+599}, /* i=150 */
   {0x1.a312ac622c53ap+8, 0x1.82890ac4fb296p+603, 0x1.82890ac4fb296p+603}, /* i=151 */
   {0x1.a5d927976167fp+8, 0x1.83976e92688ddp+607, 0x1.83976e92688ddp+607}, /* i=152 */
   {0x1.a89fa2cc95c86p+8, 0x1.84a68f8386431p+611, 0x1.84a68f8386431p+611}, /* i=153 */
   {0x1.ab661e01caef7p+8, 0x1.85b66e1ee19d1p+615, 0x1.85b66e1ee19d1p+615}, /* i=154 */
   {0x1.ae2c9936ffe77p+8, 0x1.86c70ae7b1f83p+619, 0x1.86c70ae7b1f83p+619}, /* i=155 */
   {0x1.b0f3146c349f8p+8, 0x1.87d86662e5a48p+623, 0x1.87d86662e5a48p+623}, /* i=156 */
   {0x1.b3b98fa169b81p+8, 0x1.88ea8116d7dd5p+627, 0x1.88ea8116d7dd5p+627}, /* i=157 */
   {0x1.b6800ad69e4b5p+8, 0x1.89fd5b87eb43ep+631, 0x1.89fd5b87eb43ep+631}, /* i=158 */
   {0x1.b946860bd32a8p+8, 0x1.8b10f63d81c21p+635, 0x1.8b10f63d81c21p+635}, /* i=159 */
   {0x1.bc0d014107f81p+8, 0x1.8c2551bd8a045p+639, 0x1.8c2551bd8a045p+639}, /* i=160 */
   {0x1.bed37c763ccdap+8, 0x1.8d3a6e8f0828dp+643, 0x1.8d3a6e8f0828dp+643}, /* i=161 */
   {0x1.c199f7ab71a38p+8, 0x1.8e504d392b3e9p+647, 0x1.8e504d392b3e9p+647}, /* i=162 */
   {0x1.c46072e0a6921p+8, 0x1.8f66ee43b2d1cp+651, 0x1.8f66ee43b2d1cp+651}, /* i=163 */
   {0x1.c726ee15db429p+8, 0x1.907e52360fdbp+655, 0x1.907e52360fdbp+655}, /* i=164 */
   {0x1.c9ed694b10581p+8, 0x1.9196799998d68p+659, 0x1.9196799998d68p+659}, /* i=165 */
   {0x1.ccb3e4804537ap+8, 0x1.92af64f61234fp+663, 0x1.92af64f61234fp+663}, /* i=166 */
   {0x1.cf7a5fb57a075p+8, 0x1.93c914d4cdcbfp+667, 0x1.93c914d4cdcbfp+667}, /* i=167 */
   {0x1.d240daeaaec93p+8, 0x1.94e389bf450d8p+671, 0x1.94e389bf450d8p+671}, /* i=168 */
   {0x1.d507561fe3d3ap+8, 0x1.95fec43fd7536p+675, 0x1.95fec43fd7536p+675}, /* i=169 */
   {0x1.d7cdd15518876p+8, 0x1.971ac4dfbf47p+679, 0x1.971ac4dfbf47p+679}, /* i=170 */
   {0x1.da944c8a4d2c4p+8, 0x1.98378c2a05f15p+683, 0x1.98378c2a05f15p+683}, /* i=171 */
   {0x1.dd5ac7bf8234fp+8, 0x1.99551aaa5b1bp+687, 0x1.99551aaa5b1bp+687}, /* i=172 */
   {0x1.e02142f4b6fa1p+8, 0x1.9a7370eb0d876p+691, 0x1.9a7370eb0d876p+691}, /* i=173 */
   {0x1.e2e7be29ebf9fp+8, 0x1.9b928f78a2b1p+695, 0x1.9b928f78a2b1p+695}, /* i=174 */
   {0x1.e5ae395f20e5p+8, 0x1.9cb276deb8243p+699, 0x1.9cb276deb8243p+699}, /* i=175 */
   {0x1.e874b49455b63p+8, 0x1.9dd327a9c4789p+703, 0x1.9dd327a9c4789p+703}, /* i=176 */
   {0x1.eb3b2fc98a87bp+8, 0x1.9ef4a266d31d9p+707, 0x1.9ef4a266d31d9p+707}, /* i=177 */
   {0x1.ee01aafebf599p+8, 0x1.a016e7a3280f1p+711, 0x1.a016e7a3280f1p+711}, /* i=178 */
   {0x1.f0c82633f4588p+8, 0x1.a139f7ecb2d56p+715, 0x1.a139f7ecb2d56p+715}, /* i=179 */
   {0x1.f38ea169292fbp+8, 0x1.a25dd3d0f338bp+719, 0x1.a25dd3d0f338bp+719}, /* i=180 */
   {0x1.f6551c9e5dd4ap+8, 0x1.a3827bde44a0ep+723, 0x1.a3827bde44a0ep+723}, /* i=181 */
   {0x1.f91b97d392a7cp+8, 0x1.a4a7f0a4159e3p+727, 0x1.a4a7f0a4159e3p+727}, /* i=182 */
   {0x1.fbe21308c821p+8, 0x1.a5ce32b2611adp+731, 0x1.a5ce32b2611adp+731}, /* i=183 */
   {0x1.fea88e3dfc7bp+8, 0x1.a6f54295d5209p+735, 0x1.a6f54295d5209p+735}, /* i=184 */
   {0x1.00b784b998d5ap+9, 0x1.a81d20e1a9ef2p+739, 0x1.a81d20e1a9ef2p+739}, /* i=185 */
   {0x1.021ac25433016p+9, 0x1.a945ce237878p+743, 0x1.a945ce237878p+743}, /* i=186 */
   {0x1.037dffeecd503p+9, 0x1.aa6f4aede2abfp+747, 0x1.aa6f4aede2abfp+747}, /* i=187 */
   {0x1.04e13d89680b3p+9, 0x1.ab9997d30a13dp+751, 0x1.ab9997d30a13dp+751}, /* i=188 */
   {0x1.06447b2402339p+9, 0x1.acc4b5612df02p+755, 0x1.acc4b5612df02p+755}, /* i=189 */
   {0x1.07a7b8be9ccc4p+9, 0x1.adf0a42da7a8dp+759, 0x1.adf0a42da7a8dp+759}, /* i=190 */
   {0x1.090af6593739fp+9, 0x1.af1d64c8cd02bp+763, 0x1.af1d64c8cd02bp+763}, /* i=191 */
   {0x1.0a6e33f3d1835p+9, 0x1.b04af7c57791p+767, 0x1.b04af7c57791p+767}, /* i=192 */
   {0x1.0bd1718e6c249p+9, 0x1.b1795db875b98p+771, 0x1.b1795db875b98p+771}, /* i=193 */
   {0x1.0d34af2906dacp+9, 0x1.b2a8973478677p+775, 0x1.b2a8973478677p+775}, /* i=194 */
   {0x1.0e97ecc3a1334p+9, 0x1.b3d8a4cbf2e7p+779, 0x1.b3d8a4cbf2e7p+779}, /* i=195 */
   {0x1.0ffb2a5e3b105p+9, 0x1.b5098712dd8ecp+783, 0x1.b5098712dd8ecp+783}, /* i=196 */
   {0x1.115e67f8d5b8p+9, 0x1.b63b3ea25990dp+787, 0x1.b63b3ea25990dp+787}, /* i=197 */
   {0x1.12c1a59370407p+9, 0x1.b76dcc0c7a184p+791, 0x1.b76dcc0c7a184p+791}, /* i=198 */
   {0x1.1424e32e0a74p+9, 0x1.b8a12fe61e6cep+795, 0x1.b8a12fe61e6cep+795}, /* i=199 */
   {0x1.158820c8a506ep+9, 0x1.b9d56ac7b38aap+799, 0x1.b9d56ac7b38aap+799}, /* i=200 */
   {0x1.16eb5e633f595p+9, 0x1.bb0a7d457d889p+803, 0x1.bb0a7d457d889p+803}, /* i=201 */
   {0x1.184e9bfdd9fdap+9, 0x1.bc4067f849efep+807, 0x1.bc4067f849efep+807}, /* i=202 */
   {0x1.19b1d998746d9p+9, 0x1.bd772b7584ebep+811, 0x1.bd772b7584ebep+811}, /* i=203 */
   {0x1.1b1517330ed5fp+9, 0x1.beaec855703b8p+815, 0x1.beaec855703b8p+815}, /* i=204 */
   {0x1.1c7854cda9b49p+9, 0x1.bfe73f31d4d6fp+819, 0x1.bfe73f31d4d6fp+819}, /* i=205 */
   {0x1.1ddb926843a93p+9, 0x1.c120909e5d77dp+823, 0x1.c120909e5d77dp+823}, /* i=206 */
   {0x1.1f3ed002ddf84p+9, 0x1.c25abd385edeap+827, 0x1.c25abd385edeap+827}, /* i=207 */
   {0x1.20a20d9d786e9p+9, 0x1.c395c59877fbep+831, 0x1.c395c59877fbep+831}, /* i=208 */
   {0x1.22054b3812df8p+9, 0x1.c4d1aa57c8412p+835, 0x1.c4d1aa57c8412p+835}, /* i=209 */
   {0x1.236888d2ad239p+9, 0x1.c60e6c0feb206p+839, 0x1.c60e6c0feb206p+839}, /* i=210 */
   {0x1.24cbc66d47af4p+9, 0x1.c74c0b5d113dcp+843, 0x1.c74c0b5d113dcp+843}, /* i=211 */
   {0x1.262f0407e228fp+9, 0x1.c88a88d8fcb32p+847, 0x1.c88a88d8fcb32p+847}, /* i=212 */
   {0x1.279241a27c875p+9, 0x1.c9c9e51ef7c7fp+851, 0x1.c9c9e51ef7c7fp+851}, /* i=213 */
   {0x1.28f57f3d16d14p+9, 0x1.cb0a20caf32bcp+855, 0x1.cb0a20caf32bcp+855}, /* i=214 */
   {0x1.2a58bcd7b130dp+9, 0x1.cc4b3c79cb5fep+859, 0x1.cc4b3c79cb5fep+859}, /* i=215 */
   {0x1.2bbbfa724bcb7p+9, 0x1.cd8d38c8bb2fap+863, 0x1.cd8d38c8bb2fap+863}, /* i=216 */
   {0x1.2d1f380ce6289p+9, 0x1.ced0165331358p+867, 0x1.ced0165331358p+867}, /* i=217 */
   {0x1.2e8275a780846p+9, 0x1.d013d5b79488dp+871, 0x1.d013d5b79488dp+871}, /* i=218 */
   {0x1.2fe5b3421ade1p+9, 0x1.d1587793df13cp+875, 0x1.d1587793df13cp+875}, /* i=219 */
   {0x1.3148f0dcb5a66p+9, 0x1.d29dfc881738bp+879, 0x1.d29dfc881738bp+879}, /* i=220 */
   {0x1.32ac2e774fdb3p+9, 0x1.d3e4652f6af66p+883, 0x1.d3e4652f6af66p+883}, /* i=221 */
   {0x1.340f6c11ea3c5p+9, 0x1.d52bb22bde51fp+887, 0x1.d52bb22bde51fp+887}, /* i=222 */
   {0x1.3572a9ac84f15p+9, 0x1.d673e41dbcd16p+891, 0x1.d673e41dbcd16p+891}, /* i=223 */
   {0x1.36d5e7471f262p+9, 0x1.d7bcfba223f3cp+895, 0x1.d7bcfba223f3cp+895}, /* i=224 */
   {0x1.383924e1b9befp+9, 0x1.d906f95cf3ae1p+899, 0x1.d906f95cf3ae1p+899}, /* i=225 */
   {0x1.399c627c54238p+9, 0x1.da51dded03bc8p+903, 0x1.da51dded03bc8p+903}, /* i=226 */
   {0x1.3affa016ee7fap+9, 0x1.db9da9f46edf7p+907, 0x1.db9da9f46edf7p+907}, /* i=227 */
   {0x1.3c62ddb188e1ep+9, 0x1.dcea5e1555958p+911, 0x1.dcea5e1555958p+911}, /* i=228 */
   {0x1.3dc61b4c236e7p+9, 0x1.de37faf29b034p+915, 0x1.de37faf29b034p+915}, /* i=229 */
   {0x1.3f2958e6bddc1p+9, 0x1.df86812dfb5e7p+919, 0x1.df86812dfb5e7p+919}, /* i=230 */
   {0x1.408c968158ccfp+9, 0x1.e0d5f16d162e8p+923, 0x1.e0d5f16d162e8p+923}, /* i=231 */
   {0x1.41efd41bf28afp+9, 0x1.e2264c4d2e935p+927, 0x1.e2264c4d2e935p+927}, /* i=232 */
   {0x1.435311b68d027p+9, 0x1.e3779279a2beap+931, 0x1.e3779279a2beap+931}, /* i=233 */
   {0x1.44b64f512772ap+9, 0x1.e4c9c49434e4fp+935, 0x1.e4c9c49434e4fp+935}, /* i=234 */
   {0x1.46198cebc1e11p+9, 0x1.e61ce34203caep+939, 0x1.e61ce34203caep+939}, /* i=235 */
   {0x1.477cca865c622p+9, 0x1.e770ef28da4bp+943, 0x1.e770ef28da4bp+943}, /* i=236 */
   {0x1.48e00820f6ba9p+9, 0x1.e8c5e8edc7cf4p+947, 0x1.e8c5e8edc7cf4p+947}, /* i=237 */
   {0x1.4a4345bb913a5p+9, 0x1.ea1bd138624fcp+951, 0x1.ea1bd138624fcp+951}, /* i=238 */
   {0x1.4ba683562b967p+9, 0x1.eb72a8ae64defp+955, 0x1.eb72a8ae64defp+955}, /* i=239 */
   {0x1.4d09c0f0c610fp+9, 0x1.ecca6ff81abc6p+959, 0x1.ecca6ff81abc6p+959}, /* i=240 */
   {0x1.4e6cfe8b60c89p+9, 0x1.ee2327bdbee11p+963, 0x1.ee2327bdbee11p+963}, /* i=241 */
   {0x1.4fd03c25faeb4p+9, 0x1.ef7cd0a45ebeep+967, 0x1.ef7cd0a45ebeep+967}, /* i=242 */
   {0x1.513379c09551ep+9, 0x1.f0d76b57edc45p+971, 0x1.f0d76b57edc45p+971}, /* i=243 */
   {0x1.5296b75b2fc51p+9, 0x1.f232f880bae65p+975, 0x1.f232f880bae65p+975}, /* i=244 */
   {0x1.53f9f4f5ca08p+9, 0x1.f38f78c773da9p+979, 0x1.f38f78c773da9p+979}, /* i=245 */
   {0x1.555d3290646a6p+9, 0x1.f4ececd7603b9p+983, 0x1.f4ececd7603b9p+983}, /* i=246 */
   {0x1.56c0702aff9dap+9, 0x1.f64b555dc1981p+987, 0x1.f64b555dc1981p+987}, /* i=247 */
   {0x1.5823adc59987ap+9, 0x1.f7aab2fd58bcep+991, 0x1.f7aab2fd58bcep+991}, /* i=248 */
   {0x1.5986eb6033de5p+9, 0x1.f90b0668522c8p+995, 0x1.f90b0668522c8p+995}, /* i=249 */
   {0x1.5aea28face48p+9, 0x1.fa6c50493dac3p+999, 0x1.fa6c50493dac3p+999}, /* i=250 */
   {0x1.5c4d669568b09p+9, 0x1.fbce914c33636p+1003, 0x1.fbce914c33636p+1003}, /* i=251 */
   {0x1.5db0a4300343bp+9, 0x1.fd31ca1ec10e9p+1007, 0x1.fd31ca1ec10e9p+1007}, /* i=252 */
   {0x1.5f13e1ca9e46ep+9, 0x1.fe95fb6f56358p+1011, 0x1.fe95fb6f56358p+1011}, /* i=253 */
   {0x1.60771f6538d28p+9, 0x1.fffb25e82a37ap+1015, 0x1.fffb25e82a37ap+1015}, /* i=254 */
   {0x1.61da5cffd2007p+9, 0x1.00b0a519f4fa6p+1020, 0x1.00b0a519f4fa6p+1020}, /* i=255 */
};
#endif

static const double U[256][3] = {
   {0x0p+0, 0x0p+0, 0x1p+0}, /* i=0 */
   {0x1.633d9a98c86abp-7, 0x1.633f62a20b333p-7, 0x1.0003d9ea4bffdp+0}, /* i=1 */
   {0x1.633d9a9a4ead4p-6, 0x1.6344bac795757p-6, 0x1.000f67c6da1acp+0}, /* i=2 */
   {0x1.0a6e33f18c63cp-5, 0x1.0a7a3a551b2d4p-5, 0x1.0022a9eea8778p+0}, /* i=3 */
   {0x1.633d9a9a45eabp-5, 0x1.635a1bd31b5efp-5, 0x1.003da0f60f265p+0}, /* i=4 */
   {0x1.bc0d01405faf3p-5, 0x1.bc44ae645943cp-5, 0x1.00604dacbe961p+0}, /* i=5 */
   {0x1.0a6e33f3e7daap-4, 0x1.0a9e4f7654b9bp-4, 0x1.008ab11dccb5fp+0}, /* i=6 */
   {0x1.36d5e747f76a8p-4, 0x1.37224d591c0bdp-4, 0x1.00bccc8fb84c9p+0}, /* i=7 */
   {0x1.633d9a99ec38bp-4, 0x1.63afa7b941a63p-4, 0x1.00f6a184710c1p+0}, /* i=8 */
   {0x1.8fa54ded0b802p-4, 0x1.9047b5c35b96fp-4, 0x1.013831b96d97p+0}, /* i=9 */
   {0x1.bc0d0141af0fdp-4, 0x1.bcebcef38197ap-4, 0x1.01817f27aedfcp+0}, /* i=10 */
   {0x1.e874b494f349ep-4, 0x1.e99d4b1f4c7fdp-4, 0x1.01d28c03cc329p+0}, /* i=11 */
   {0x1.0a6e33f3bab52p-3, 0x1.0b2ec1437b0f6p-3, 0x1.022b5abe0f712p+0}, /* i=12 */
   {0x1.20a20d9d4d628p-3, 0x1.2196e6ee7d6b3p-3, 0x1.028bee0284c2p+0}, /* i=13 */
   {0x1.36d5e7474d5b8p-3, 0x1.3807c3273f40fp-3, 0x1.02f448b90b0d4p+0}, /* i=14 */
   {0x1.4d09c0f0433f8p-3, 0x1.4e8202c5d715ap-3, 0x1.03646e0564306p+0}, /* i=15 */
   {0x1.633d9a9a1f225p-3, 0x1.650652ee1d78cp-3, 0x1.03dc6147662e8p+0}, /* i=16 */
   {0x1.797174465dd59p-3, 0x1.7b956110ff177p-3, 0x1.045c261b044fp+0}, /* i=17 */
   {0x1.8fa54ded0ff92p-3, 0x1.922fdae88f34dp-3, 0x1.04e3c05831bf6p+0}, /* i=18 */
   {0x1.a5d92798d7cfcp-3, 0x1.a8d66e99be79p-3, 0x1.05733413a8554p+0}, /* i=19 */
   {0x1.bc0d01411268p-3, 0x1.bf89ca9310f7dp-3, 0x1.060a859e20405p+0}, /* i=20 */
   {0x1.d240dae9f4f42p-3, 0x1.d64a9db298aeep-3, 0x1.06a9b9853b1a8p+0}, /* i=21 */
   {0x1.e874b49405c15p-3, 0x1.ed19973a5f9ddp-3, 0x1.0750d4933e734p+0}, /* i=22 */
   {0x1.fea88e3bbc287p-3, 0x1.01fbb36a9fa96p-2, 0x1.07ffdbcf2f233p+0}, /* i=23 */
   {0x1.0a6e33f34a7bap-2, 0x1.0d725e54d112bp-2, 0x1.08b6d47d5ee4ep+0}, /* i=24 */
   {0x1.158820c8b6d2fp-2, 0x1.18f124a65256ap-2, 0x1.0975c41f0a294p+0}, /* i=25 */
   {0x1.20a20d9d347fp-2, 0x1.24785ee8c80eap-2, 0x1.0a3cb072cd70dp+0}, /* i=26 */
   {0x1.2bbbfa721d88bp-2, 0x1.300865e95162fp-2, 0x1.0b0b9f74f3217p+0}, /* i=27 */
   {0x1.36d5e747d2ee2p-2, 0x1.3ba192b7d4e57p-2, 0x1.0be2975f69a3cp+0}, /* i=28 */
   {0x1.41efd421d6446p-2, 0x1.47443eadf8acfp-2, 0x1.0cc19eaa43b3bp+0}, /* i=29 */
   {0x1.4d09c0f04f9fcp-2, 0x1.52f0c35a33914p-2, 0x1.0da8bc0a1977cp+0}, /* i=30 */
   {0x1.5823adc5b71fcp-2, 0x1.5ea77abaded9bp-2, 0x1.0e97f6748d816p+0}, /* i=31 */
   {0x1.633d9a9a861e2p-2, 0x1.6a68bf01925f5p-2, 0x1.0f8f551ba9fe7p+0}, /* i=32 */
   {0x1.6e57876f71e8ep-2, 0x1.7634eab9d1244p-2, 0x1.108edf70d95d1p+0}, /* i=33 */
   {0x1.7971744450862p-2, 0x1.820c58c22f55fp-2, 0x1.11969d24622c8p+0}, /* i=34 */
   {0x1.848b6119879d3p-2, 0x1.8def64509375ep-2, 0x1.12a69625c2c2dp+0}, /* i=35 */
   {0x1.8fa54dede59d9p-2, 0x1.99de68f2a691p-2, 0x1.13bed2a3b9539p+0}, /* i=36 */
   {0x1.9abf3ac2b3ce9p-2, 0x1.a5d9c294f3a96p-2, 0x1.14df5b0ce8585p+0}, /* i=37 */
   {0x1.a5d92796b6e87p-2, 0x1.b1e1cd804c30ep-2, 0x1.1608380f9945ap+0}, /* i=38 */
   {0x1.b0f3146bb22e8p-2, 0x1.bdf6e662836bap-2, 0x1.1739729a8f951p+0}, /* i=39 */
   {0x1.bc0d014105485p-2, 0x1.ca196a4b6466dp-2, 0x1.187313dcbebffp+0}, /* i=40 */
   {0x1.c726ee15e8f4fp-2, 0x1.d649b6b1e68efp-2, 0x1.19b52545c9b2bp+0}, /* i=41 */
   {0x1.d240daeab2992p-2, 0x1.e2882978621cdp-2, 0x1.1affb08670ff9p+0}, /* i=42 */
   {0x1.dd5ac7c01995p-2, 0x1.eed520ee9bp-2, 0x1.1c52bf90cb347p+0}, /* i=43 */
   {0x1.e874b49449e79p-2, 0x1.fb30fbd155be6p-2, 0x1.1dae5c9836286p+0}, /* i=44 */
   {0x1.f38ea168ef806p-2, 0x1.03ce0ca9f3723p-1, 0x1.1f12921260006p+0}, /* i=45 */
   {0x1.fea88e3d9f9cep-2, 0x1.0a0b6c8e6e28dp-1, 0x1.207f6ab6eb7a4p+0}, /* i=46 */
   {0x1.04e13d895812ep-1, 0x1.1050cda60393bp-1, 0x1.21f4f1801b11ap+0}, /* i=47 */
   {0x1.0a6e33f42f5f8p-1, 0x1.169e603e37552p-1, 0x1.237331ab0afdp+0}, /* i=48 */
   {0x1.0ffb2a5e4c8c1p-1, 0x1.1cf454e257075p-1, 0x1.24fa36b7b8652p+0}, /* i=49 */
   {0x1.158820c91514dp-1, 0x1.2352dc61081b4p-1, 0x1.268a0c6a56215p+0}, /* i=50 */
   {0x1.1b15173257ebp-1, 0x1.29ba27c6dbfa7p-1, 0x1.2822bec9fa9dfp+0}, /* i=51 */
   {0x1.20a20d9d79ecdp-1, 0x1.302a686a58bf5p-1, 0x1.29c45a2399591p+0}, /* i=52 */
   {0x1.262f040562b2dp-1, 0x1.36a3cfdd2c5d3p-1, 0x1.2b6eeb06498bap+0}, /* i=53 */
   {0x1.2bbbfa6ca4527p-1, 0x1.3d2690006480ep-1, 0x1.2d227e4863137p+0}, /* i=54 */
   {0x1.3148f0dc3360ap-1, 0x1.43b2db04ab666p-1, 0x1.2edf2107b7e66p+0}, /* i=55 */
   {0x1.36d5e746d8077p-1, 0x1.4a48e34a0780cp-1, 0x1.30a4e0a0f7ba6p+0}, /* i=56 */
   {0x1.3c62ddb2d88fdp-1, 0x1.50e8db927b60bp-1, 0x1.3273cabd15bebp+0}, /* i=57 */
   {0x1.41efd41bcc36p-1, 0x1.5792f6dfe4274p-1, 0x1.344bed4832e5ap+0}, /* i=58 */
   {0x1.477cca862041ap-1, 0x1.5e47688cc6ff4p-1, 0x1.362d56785d9b9p+0}, /* i=59 */
   {0x1.4d09c0ef32d9dp-1, 0x1.6506643acecb6p-1, 0x1.381814c8be821p+0}, /* i=60 */
   {0x1.5296b75ae00ap-1, 0x1.6bd01de4b6363p-1, 0x1.3a0c36fe9d0f6p+0}, /* i=61 */
   {0x1.5823adc4f6df3p-1, 0x1.72a4c9ce2392ap-1, 0x1.3c09cc24c4e9dp+0}, /* i=62 */
   {0x1.5db0a4313aae3p-1, 0x1.79849c98ce061p-1, 0x1.3e10e3919b977p+0}, /* i=63 */
   {0x1.633d9a9aa2021p-1, 0x1.806fcb316fd34p-1, 0x1.40218ce190b08p+0}, /* i=64 */
   {0x1.68ca910445a8p-1, 0x1.87668ae653b8p-1, 0x1.423bd7fdc2f6dp+0}, /* i=65 */
   {0x1.6e57876ceec45p-1, 0x1.8e6911597fb6fp-1, 0x1.445fd517e9febp+0}, /* i=66 */
   {0x1.73e47dd890f7p-1, 0x1.9577948e3b3c5p-1, 0x1.468d94ae6ded2p+0}, /* i=67 */
   {0x1.79717444796f5p-1, 0x1.9c924adbbb041p-1, 0x1.48c527885172ep+0}, /* i=68 */
   {0x1.7efe6aaeaef84p-1, 0x1.a3b96af803a88p-1, 0x1.4b069eb87df6p+0}, /* i=69 */
   {0x1.848b611810801p-1, 0x1.aaed2bfc58c9ap-1, 0x1.4d520b9f2f7f4p+0}, /* i=70 */
   {0x1.8a18578296a7fp-1, 0x1.b22dc564bc7bep-1, 0x1.4fa77fe9d6cdp+0}, /* i=71 */
   {0x1.8fa54deda7092p-1, 0x1.b97b6f0ccdcf1p-1, 0x1.52070d9216085p+0}, /* i=72 */
   {0x1.95324457957fep-1, 0x1.c0d66133592c9p-1, 0x1.5470c6dedcb1fp+0}, /* i=73 */
   {0x1.9abf3ac15ae4cp-1, 0x1.c83ed480f5b54p-1, 0x1.56e4be66920d8p+0}, /* i=74 */
   {0x1.a04c312cee8f1p-1, 0x1.cfb502079fe9ap-1, 0x1.5963070efd69p+0}, /* i=75 */
   {0x1.a5d9279724603p-1, 0x1.d739233c287aap-1, 0x1.5bebb40b0ea73p+0}, /* i=76 */
   {0x1.ab661e01b6bfcp-1, 0x1.decb72056c3f5p-1, 0x1.5e7ed8dff377p+0}, /* i=77 */
   {0x1.b0f3146c697c4p-1, 0x1.e66c28b4d93f4p-1, 0x1.611c896299448p+0}, /* i=78 */
   {0x1.b6800ad6a87c1p-1, 0x1.ee1b820a56d4fp-1, 0x1.63c4d9b8fb767p+0}, /* i=79 */
   {0x1.bc0d0140fe8efp-1, 0x1.f5d9b93812557p-1, 0x1.6677de5b72f7ep+0}, /* i=80 */
   {0x1.c199f7ab79a4ap-1, 0x1.fda709e211184p-1, 0x1.6935ac1491c72p+0}, /* i=81 */
   {0x1.c726ee156f7a3p-1, 0x1.02c1d80fd077ap+0, 0x1.6bfe5801a068ap+0}, /* i=82 */
   {0x1.ccb3e4808cdb9p-1, 0x1.06b7f440bf9a1p+0, 0x1.6ed1f794cb557p+0}, /* i=83 */
   {0x1.d240daeac2fbap-1, 0x1.0ab5f80583b64p+0, 0x1.71b0a092b58b9p+0}, /* i=84 */
   {0x1.d7cdd154f7505p-1, 0x1.0ebc021ec7a18p+0, 0x1.749a691684294p+0}, /* i=85 */
   {0x1.dd5ac7bf96959p-1, 0x1.12ca318ab2388p+0, 0x1.778f6790cfa16p+0}, /* i=86 */
   {0x1.e2e7be29aa4f7p-1, 0x1.16e0a585245d2p+0, 0x1.7a8fb2c7cc2acp+0}, /* i=87 */
   {0x1.e874b4942365fp-1, 0x1.1aff7d8b131cp+0, 0x1.7d9b61d9bcb1fp+0}, /* i=88 */
   {0x1.ee01aaff2af5bp-1, 0x1.1f26d959845bep+0, 0x1.80b28c3c371bep+0}, /* i=89 */
   {0x1.f38ea16a05998p-1, 0x1.2356d8ee67ab7p+0, 0x1.83d549bcc09d1p+0}, /* i=90 */
   {0x1.f91b97d3862d1p-1, 0x1.278f9c89dd44ep+0, 0x1.8703b281bdeebp+0}, /* i=91 */
   {0x1.fea88e3d9f79p-1, 0x1.2bd144b1ea44ap+0, 0x1.8a3ddf0d3f846p+0}, /* i=92 */
   {0x1.021ac2543dc5bp+0, 0x1.301bf22fbcc3fp+0, 0x1.8d83e83af532dp+0}, /* i=93 */
   {0x1.04e13d88f963fp+0, 0x1.346fc6102ed01p+0, 0x1.90d5e7408833ep+0}, /* i=94 */
   {0x1.07a7b8bd5b1b3p+0, 0x1.38cce1a8c5752p+0, 0x1.9433f5b16bf22p+0}, /* i=95 */
   {0x1.0a6e33f3fd01cp+0, 0x1.3d33669a1cd23p+0, 0x1.979e2d80ca2bcp+0}, /* i=96 */
   {0x1.0d34af2904b1ep+0, 0x1.41a376c31539ap+0, 0x1.9b14a8f78b123p+0}, /* i=97 */
   {0x1.0ffb2a5e48ec5p+0, 0x1.461d3455143f7p+0, 0x1.9e9782c409bbep+0}, /* i=98 */
   {0x1.12c1a5931363dp+0, 0x1.4aa0c1c7f8a1bp+0, 0x1.a226d5f0bb389p+0}, /* i=99 */
   {0x1.158820c8f4c5p+0, 0x1.4f2e41e2eb096p+0, 0x1.a5c2bdeb17d8cp+0}, /* i=100 */
   {0x1.184e9bfd34f8ap+0, 0x1.53c5d7b2abc2ep+0, 0x1.a96b567becf85p+0}, /* i=101 */
   {0x1.1b1517334db7cp+0, 0x1.5867a69be2f0ap+0, 0x1.ad20bbd5e9772p+0}, /* i=102 */
   {0x1.1ddb9268056dp+0, 0x1.5d13d246d5f83p+0, 0x1.b0e30a85710d8p+0}, /* i=103 */
   {0x1.20a20d9d504bep+0, 0x1.61ca7eb44013ap+0, 0x1.b4b25f8146993p+0}, /* i=104 */
   {0x1.236888d351b2dp+0, 0x1.668bd032c0575p+0, 0x1.b88ed822190b2p+0}, /* i=105 */
   {0x1.262f04082ee1ep+0, 0x1.6b57eb5f8031ep+0, 0x1.bc789222fa6c8p+0}, /* i=106 */
   {0x1.28f57f3d53f74p+0, 0x1.702ef5306120fp+0, 0x1.c06faba9a8ec3p+0}, /* i=107 */
   {0x1.2bbbfa717d242p+0, 0x1.751112eacfb55p+0, 0x1.c474433f153fp+0}, /* i=108 */
   {0x1.2e8275a25eb1dp+0, 0x1.79fe6a27a0ed9p+0, 0x1.c88677d27e5efp+0}, /* i=109 */
   {0x1.3148f0dc8f6efp+0, 0x1.7ef720f0c1916p+0, 0x1.cca668d214e81p+0}, /* i=110 */
   {0x1.340f6c11f0d43p+0, 0x1.83fb5d780a8e2p+0, 0x1.d0d435ee57893p+0}, /* i=111 */
   {0x1.36d5e7468fa39p+0, 0x1.890b466873bddp+0, 0x1.d50fff5d503fcp+0}, /* i=112 */
   {0x1.399c627c0c6f6p+0, 0x1.8e2702c2a84adp+0, 0x1.d959e5bd4f3dcp+0}, /* i=113 */
   {0x1.3c62ddb05aa94p+0, 0x1.934eb9dba9be8p+0, 0x1.ddb20a13afe29p+0}, /* i=114 */
   {0x1.3f2958e7eda0bp+0, 0x1.98829370d8c9ap+0, 0x1.e2188dddb5348p+0}, /* i=115 */
   {0x1.41efd41c728c5p+0, 0x1.9dc2b788e26c8p+0, 0x1.e68d92f64a42p+0}, /* i=116 */
   {0x1.44b64f50e5ba8p+0, 0x1.a30f4e9a89f8ep+0, 0x1.eb113bb6cbc75p+0}, /* i=117 */
   {0x1.477cca86c4b06p+0, 0x1.a8688179c6c07p+0, 0x1.efa3aae71d101p+0}, /* i=118 */
   {0x1.4a4345bd16336p+0, 0x1.adce7956f9036p+0, 0x1.f44503bcec968p+0}, /* i=119 */
   {0x1.4d09c0f1bc76dp+0, 0x1.b3415fc2896a2p+0, 0x1.f8f569debfdd2p+0}, /* i=120 */
   {0x1.4fd03c258d664p+0, 0x1.b8c15eb634c7ap+0, 0x1.fdb5016bf30fbp+0}, /* i=121 */
   {0x1.5296b75b0f115p+0, 0x1.be4ea0940abc5p+0, 0x1.0141f77df5f77p+1}, /* i=122 */
   {0x1.555d3290420a4p+0, 0x1.c3e9501b447dap+0, 0x1.03b12bca368eap+1}, /* i=123 */
   {0x1.5823adc56efe6p+0, 0x1.c99198770f432p+0, 0x1.0628305b2edbbp+1}, /* i=124 */
   {0x1.5aea28f9db4e1p+0, 0x1.cf47a53940a8fp+0, 0x1.08a7182c99fc7p+1}, /* i=125 */
   {0x1.5db0a42e7458ep+0, 0x1.d50ba2611f2e3p+0, 0x1.0b2df678704dfp+1}, /* i=126 */
   {0x1.60771f642df45p+0, 0x1.daddbc596f107p+0, 0x1.0dbcdeb61023p+1}, /* i=127 */
   {0x1.633d9a99998b7p+0, 0x1.e0be1ff4a926ap+0, 0x1.1053e4989116ap+1}, /* i=128 */
   {0x1.660415d021e29p+0, 0x1.e6acfa79542cdp+0, 0x1.12f31c143238bp+1}, /* i=129 */
   {0x1.68ca910458747p+0, 0x1.ecaa7992fa3bap+0, 0x1.159a9957b5bfbp+1}, /* i=130 */
   {0x1.6b910c3a2b6f6p+0, 0x1.f2b6cb6dc49ebp+0, 0x1.184a70d8988ep+1}, /* i=131 */
   {0x1.6e57876d7cb6p+0, 0x1.f8d21e9692bb5p+0, 0x1.1b02b744ed615p+1}, /* i=132 */
   {0x1.711e02a57a2efp+0, 0x1.fefca226b6188p+0, 0x1.1dc38196cfd66p+1}, /* i=133 */
   {0x1.73e47dda2e13ap+0, 0x1.029b42c57e28ap+1, 0x1.208ce4faffc89p+1}, /* i=134 */
   {0x1.76aaf90ead5b7p+0, 0x1.05bffc5f87d82p+1, 0x1.235ef6eb5dd67p+1}, /* i=135 */
   {0x1.79717445d679cp+0, 0x1.08ec961b1fb43p+1, 0x1.2639cd2497d3ep+1}, /* i=136 */
   {0x1.7c37ef793a5fep+0, 0x1.0c212863eb8aap+1, 0x1.291d7d9d486a7p+1}, /* i=137 */
   {0x1.7efe6aaecdfaep+0, 0x1.0f5dcbf131853p+1, 0x1.2c0a1e9d09a77p+1}, /* i=138 */
   {0x1.81c4e5e2eb71dp+0, 0x1.12a299ad42ea7p+1, 0x1.2effc6a6503cp+1}, /* i=139 */
   {0x1.848b6118a797p+0, 0x1.15efaac93475ep+1, 0x1.31fe8c88381bdp+1}, /* i=140 */
   {0x1.8751dc4d9cd23p+0, 0x1.194518af48916p+1, 0x1.35068752397e6p+1}, /* i=141 */
   {0x1.8a185782e6e46p+0, 0x1.1ca2fd0e53f17p+1, 0x1.3817ce5e7c62cp+1}, /* i=142 */
   {0x1.8cded2b8bfc87p+0, 0x1.200971d552106p+1, 0x1.3b32794dd7846p+1}, /* i=143 */
   {0x1.8fa54dedd33e1p+0, 0x1.2378913349073p+1, 0x1.3e56a007b3385p+1}, /* i=144 */
   {0x1.926bc922a355cp+0, 0x1.26f0759c3d796p+1, 0x1.41845abe8f5a1p+1}, /* i=145 */
   {0x1.95324459ce4d8p+0, 0x1.2a7139ca73fbap+1, 0x1.44bbc1f12ed23p+1}, /* i=146 */
   {0x1.97f8bf8d8935ep+0, 0x1.2dfaf8b1e6a39p+1, 0x1.47fcee5f108f9p+1}, /* i=147 */
   {0x1.9abf3ac273558p+0, 0x1.318dcd98e7fc4p+1, 0x1.4b47f91f11465p+1}, /* i=148 */
   {0x1.9d85b5f6cd85cp+0, 0x1.3529d403b0235p+1, 0x1.4e9cfb8ca17f3p+1}, /* i=149 */
   {0x1.a04c312cb394ep+0, 0x1.38cf27c246d92p+1, 0x1.51fc0f5494be3p+1}, /* i=150 */
   {0x1.a312ac61c284p+0, 0x1.3c7de4e677d8ap+1, 0x1.55654e6bdaa4bp+1}, /* i=151 */
   {0x1.a5d927967f125p+0, 0x1.403627ce5678ep+1, 0x1.58d8d31935e16p+1}, /* i=152 */
   {0x1.a89fa2cd06b3p+0, 0x1.43f80d2382078p+1, 0x1.5c56b7f494425p+1}, /* i=153 */
   {0x1.ab661e016e378p+0, 0x1.47c3b1d1a2063p+1, 0x1.5fdf17de350a4p+1}, /* i=154 */
   {0x1.ae2c9936a51bap+0, 0x1.4b993318bf12fp+1, 0x1.63720e0fb6d8fp+1}, /* i=155 */
   {0x1.b0f3146b13af2p+0, 0x1.4f78ae7ea778fp+1, 0x1.670fb60e82129p+1}, /* i=156 */
   {0x1.b3b98fa1dcdfcp+0, 0x1.536241dc9ec12p+1, 0x1.6ab82bb88b28cp+1}, /* i=157 */
   {0x1.b6800ad6924d5p+0, 0x1.57560b4f0096ap+1, 0x1.6e6b8b3508169p+1}, /* i=158 */
   {0x1.b946860a8a5bcp+0, 0x1.5b542948e52abp+1, 0x1.7229f106cba3dp+1}, /* i=159 */
   {0x1.bc0d01386c475p+0, 0x1.5f5cba834642ep+1, 0x1.75f379fc7c09ep+1}, /* i=160 */
   {0x1.bed37c76759aep+0, 0x1.636fde2fffbcfp+1, 0x1.79c843607eb88p+1}, /* i=161 */
   {0x1.c199f7ab898f4p+0, 0x1.678db38d19a36p+1, 0x1.7da86a92cdae7p+1}, /* i=162 */
   {0x1.c46072e077394p+0, 0x1.6bb65a5c8d32ap+1, 0x1.81940d7983e49p+1}, /* i=163 */
   {0x1.c726ee164eee7p+0, 0x1.6fe9f2a86779ep+1, 0x1.858b4a489f88p+1}, /* i=164 */
   {0x1.c9ed694b2fb65p+0, 0x1.74289cca8ac8ep+1, 0x1.898e3f89453f8p+1}, /* i=165 */
   {0x1.ccb3e480392edp+0, 0x1.78727976a883bp+1, 0x1.8d9d0c232ac1ap+1}, /* i=166 */
   {0x1.cf7a5fb5974bbp+0, 0x1.7cc7a9b545461p+1, 0x1.91b7cf57e3659p+1}, /* i=167 */
   {0x1.d240dae9d91adp+0, 0x1.81284ee397c72p+1, 0x1.95dea8c2bdedcp+1}, /* i=168 */
   {0x1.d5075620df9cbp+0, 0x1.85948abf83204p+1, 0x1.9a11b86421699p+1}, /* i=169 */
   {0x1.d7cdd1570d6aap+0, 0x1.8a0c7f544dcdfp+1, 0x1.9e511e8f4070fp+1}, /* i=170 */
   {0x1.da944c88e925fp+0, 0x1.8e904f0784e47p+1, 0x1.a29cfbf64dbfap+1}, /* i=171 */
   {0x1.dd5ac7bf58651p+0, 0x1.93201cae3b542p+1, 0x1.a6f571bebc0a7p+1}, /* i=172 */
   {0x1.e02142f4932ddp+0, 0x1.97bc0b61fc71ap+1, 0x1.ab5aa158438ebp+1}, /* i=173 */
   {0x1.e2e7be299e7bbp+0, 0x1.9c643ea4d33ebp+1, 0x1.afccac9f2bf16p+1}, /* i=174 */
   {0x1.e5ae39608c02dp+0, 0x1.a118da5924bb9p+1, 0x1.b44bb5d48f894p+1}, /* i=175 */
   {0x1.e874b493acdccp+0, 0x1.a5da02b3b440fp+1, 0x1.b8d7df90fb665p+1}, /* i=176 */
   {0x1.eb3b2fc93a52ap+0, 0x1.aaa7dc5e3e90ap+1, 0x1.bd714ce584a5fp+1}, /* i=177 */
   {0x1.ee01aafe2f062p+0, 0x1.af828c54759d8p+1, 0x1.c218213a4d396p+1}, /* i=178 */
   {0x1.f0c82633a4455p+0, 0x1.b46a37fc058ebp+1, 0x1.c6cc80657dab2p+1}, /* i=179 */
   {0x1.f38ea16914852p+0, 0x1.b95f051bc725dp+1, 0x1.cb8e8ea2d80bp+1}, /* i=180 */
   {0x1.f6551c9e4b64cp+0, 0x1.be6119e040163p+1, 0x1.d05e709807eap+1}, /* i=181 */
   {0x1.f91b97d1e2291p+0, 0x1.c3709cda0a31ep+1, 0x1.d53c4b5319067p+1}, /* i=182 */
   {0x1.fbe2130947b75p+0, 0x1.c88db50dbf029p+1, 0x1.da284459c8bcep+1}, /* i=183 */
   {0x1.fea88e3e8c145p+0, 0x1.cdb889d3baa3fp+1, 0x1.df22818a80ea1p+1}, /* i=184 */
   {0x1.00b784ba3e917p+1, 0x1.d2f142fe39be8p+1, 0x1.e42b29410544fp+1}, /* i=185 */
   {0x1.021ac254a5d9cp+1, 0x1.d83808c1e1135p+1, 0x1.e942623fdb79p+1}, /* i=186 */
   {0x1.037dffed63178p+1, 0x1.dd8d03bea90bp+1, 0x1.ee6853b8dffbdp+1}, /* i=187 */
   {0x1.04e13d885a2cbp+1, 0x1.e2f05d1542eb3p+1, 0x1.f39d2561f3334p+1}, /* i=188 */
   {0x1.06447b2369b6bp+1, 0x1.e8623e3def6fp+1, 0x1.f8e0ff4d394e9p+1}, /* i=189 */
   {0x1.07a7b8c0374f8p+1, 0x1.ede2d12f904a5p+1, 0x1.fe340a0ed6654p+1}, /* i=190 */
   {0x1.090af658f0ed1p+1, 0x1.f372403468db1p+1, 0x1.01cb37498b1cfp+2}, /* i=191 */
   {0x1.0a6e33f586223p+1, 0x1.f910b640346ccp+1, 0x1.04842b38e28e4p+2}, /* i=192 */
   {0x1.0bd1718ddb58ep+1, 0x1.febe5e79c8d8ap+1, 0x1.0744f5bdc6b25p+2}, /* i=193 */
   {0x1.0d34af294667fp+1, 0x1.023db25dea6f1p+2, 0x1.0a0dac1b5bf2ap+2}, /* i=194 */
   {0x1.0e97ecc7026d8p+1, 0x1.0523fa9b4fe4p+2, 0x1.0cde63c1e9cb3p+2}, /* i=195 */
   {0x1.0ffb2a5f11a91p+1, 0x1.08121e39ac408p+2, 0x1.0fb7325068e72p+2}, /* i=196 */
   {0x1.115e67f8b4616p+1, 0x1.0b0833da17a6ep+2, 0x1.12982dc33148ap+2}, /* i=197 */
   {0x1.12c1a5932460ap+1, 0x1.0e06524a1d868p+2, 0x1.15816c4550a98p+2}, /* i=198 */
   {0x1.1424e32e67adp+1, 0x1.110c9096d99f4p+2, 0x1.18730443196a9p+2}, /* i=199 */
   {0x1.158820c8be7dep+1, 0x1.141b06080fb8cp+2, 0x1.1b6d0c655bbd3p+2}, /* i=200 */
   {0x1.16eb5e62beefap+1, 0x1.1731ca29db354p+2, 0x1.1e6f9b9ad3871p+2}, /* i=201 */
   {0x1.184e9bfdd0f4ap+1, 0x1.1a50f4ca3cb31p+2, 0x1.217ac915c8d9ep+2}, /* i=202 */
   {0x1.19b1d99911aebp+1, 0x1.1d789df2ce19bp+2, 0x1.248eac45ec0b8p+2}, /* i=203 */
   {0x1.1b151732b6cc8p+1, 0x1.20a8ddec6a4cap+2, 0x1.27ab5cdbe38a8p+2}, /* i=204 */
   {0x1.1c7854cd64758p+1, 0x1.23e1cd4c40082p+2, 0x1.2ad0f2d60f7b9p+2}, /* i=205 */
   {0x1.1ddb9268373b2p+1, 0x1.272384e25249p+2, 0x1.2dff866f76284p+2}, /* i=206 */
   {0x1.1f3ed002b7572p+1, 0x1.2a6e1dc333157p+2, 0x1.3137302942092p+2}, /* i=207 */
   {0x1.20a20d9d8188fp+1, 0x1.2dc1b14a68b4p+2, 0x1.347808cd1a04bp+2}, /* i=208 */
   {0x1.22054b3607618p+1, 0x1.311e5910fae94p+2, 0x1.37c22963e3803p+2}, /* i=209 */
   {0x1.236888d2681d4p+1, 0x1.34842f0af908ap+2, 0x1.3b15ab52a630bp+2}, /* i=210 */
   {0x1.24cbc66de92cfp+1, 0x1.37f34d58cd5acp+2, 0x1.3e72a82cdf2d9p+2}, /* i=211 */
   {0x1.262f0408b958cp+1, 0x1.3b6bce6e3c623p+2, 0x1.41d939daa9a4bp+2}, /* i=212 */
   {0x1.279241a4a0f3ap+1, 0x1.3eedcd0b6b68ep+2, 0x1.45497a91eecaep+2}, /* i=213 */
   {0x1.28f57f3d085c8p+1, 0x1.427964296aa6bp+2, 0x1.48c384c35142ap+2}, /* i=214 */
   {0x1.2a58bcd78d0ecp+1, 0x1.460eaf2532881p+2, 0x1.4c47734453954p+2}, /* i=215 */
   {0x1.2bbbfa723194p+1, 0x1.49adc99321096p+2, 0x1.4fd56123b3aa7p+2}, /* i=216 */
   {0x1.2d1f380cb7eecp+1, 0x1.4d56cf578ed59p+2, 0x1.536d69c184dc1p+2}, /* i=217 */
   {0x1.2e8275a6f7136p+1, 0x1.5109dca36ce98p+2, 0x1.570fa8cbdedd7p+2}, /* i=218 */
   {0x1.2fe5b3416fa8bp+1, 0x1.54c70df6b37ebp+2, 0x1.5abc3a4141d82p+2}, /* i=219 */
   {0x1.3148f0dcd9eap+1, 0x1.588e80201d939p+2, 0x1.5e733a7053bdep+2}, /* i=220 */
   {0x1.32ac2e78c52b8p+1, 0x1.5c60503a3fe98p+2, 0x1.6234c5f5039bcp+2}, /* i=221 */
   {0x1.340f6c12bbf09p+1, 0x1.603c9ba9d6b4ap+2, 0x1.6600f9b6dca5ap+2}, /* i=222 */
   {0x1.3572a9ac0cd71p+1, 0x1.6423802eb4c08p+2, 0x1.69d7f2f9ae55p+2}, /* i=223 */
   {0x1.36d5e746ca3adp+1, 0x1.68151bdca5352p+2, 0x1.6db9cf568ead1p+2}, /* i=224 */
   {0x1.383924e09c5a7p+1, 0x1.6c118d0d7d438p+2, 0x1.71a6acae227cap+2}, /* i=225 */
   {0x1.399c627d1f4bp+1, 0x1.7018f27f9b59bp+2, 0x1.759ea946a3bccp+2}, /* i=226 */
   {0x1.3affa017d0ab5p+1, 0x1.742b6b2e3cbbfp+2, 0x1.79a1e3a4d0dep+2}, /* i=227 */
   {0x1.3c62ddb17ac1cp+1, 0x1.7849167929b13p+2, 0x1.7db07ab2ff88bp+2}, /* i=228 */
   {0x1.3dc61b4bb6b75p+1, 0x1.7c721418f30fcp+2, 0x1.81ca8db589caep+2}, /* i=229 */
   {0x1.3f2958e713a58p+1, 0x1.80a6841a67dccp+2, 0x1.85f03c465437bp+2}, /* i=230 */
   {0x1.408c9680f57b1p+1, 0x1.84e686d8d3852p+2, 0x1.8a21a64f1e556p+2}, /* i=231 */
   {0x1.41efd41bdc3a9p+1, 0x1.89323d1835306p+2, 0x1.8e5eec235f0d2p+2}, /* i=232 */
   {0x1.435311b76c6d2p+1, 0x1.8d89c7ee048a4p+2, 0x1.92a82e695a84bp+2}, /* i=233 */
   {0x1.44b64f50b5692p+1, 0x1.91ed48c296996p+2, 0x1.96fd8e1b7f393p+2}, /* i=234 */
   {0x1.46198cec126ddp+1, 0x1.965ce1716e01ep+2, 0x1.9b5f2ca84e921p+2}, /* i=235 */
   {0x1.477cca862afeep+1, 0x1.9ad8b41ac6e31p+2, 0x1.9fcd2bc47e14bp+2}, /* i=236 */
   {0x1.48e008208f0ep+1, 0x1.9f60e34cb0c73p+2, 0x1.a447ad938e4e9p+2}, /* i=237 */
   {0x1.4a4345bbc5498p+1, 0x1.a3f591f13960dp+2, 0x1.a8ced4962ea4cp+2}, /* i=238 */
   {0x1.4ba683562fa36p+1, 0x1.a896e34b9498ap+2, 0x1.ad62c3a76c464p+2}, /* i=239 */
   {0x1.4d09c0f09bdd7p+1, 0x1.ad44fb0843e1dp+2, 0x1.b2039e0ca9efap+2}, /* i=240 */
   {0x1.4e6cfe8b2083bp+1, 0x1.b1fffd33d8628p+2, 0x1.b6b1876c7dad2p+2}, /* i=241 */
   {0x1.4fd03c25ec60ap+1, 0x1.b6c80e3eb63f3p+2, 0x1.bb6ca3d2692f4p+2}, /* i=242 */
   {0x1.513379c0818dap+1, 0x1.bb9d52fb80482p+2, 0x1.c03517ad498cp+2}, /* i=243 */
   {0x1.5296b75cdfc6cp+1, 0x1.c07ff0ab474bdp+2, 0x1.c50b07db66779p+2}, /* i=244 */
   {0x1.53f9f4f57a6e7p+1, 0x1.c5700cdb37b68p+2, 0x1.c9ee99887861fp+2}, /* i=245 */
   {0x1.555d329059debp+1, 0x1.ca6dcda9762ebp+2, 0x1.cedff271d45bcp+2}, /* i=246 */
   {0x1.56c0702b18b5fp+1, 0x1.cf79597f5c8bp+2, 0x1.d3df38a15bc2fp+2}, /* i=247 */
   {0x1.5823adc5863fbp+1, 0x1.d492d738493d4p+2, 0x1.d8ec9293e836cp+2}, /* i=248 */
   {0x1.5986eb604f18fp+1, 0x1.d9ba6e1e365a9p+2, 0x1.de082735ec626p+2}, /* i=249 */
   {0x1.5aea28fa942d9p+1, 0x1.def045e1fb639p+2, 0x1.e3321ddbc7e19p+2}, /* i=250 */
   {0x1.5c4d6694f75b1p+1, 0x1.e43486a7d1dccp+2, 0x1.e86a9e4e2eaa3p+2}, /* i=251 */
   {0x1.5db0a43009069p+1, 0x1.e9875902c8b02p+2, 0x1.edb1d0c5a776fp+2}, /* i=252 */
   {0x1.5f13e1ca674d5p+1, 0x1.eee8e5eea731cp+2, 0x1.f307dde47b841p+2}, /* i=253 */
   {0x1.60771f6594194p+1, 0x1.f45956e3cfac5p+2, 0x1.f86ceeca6e1d5p+2}, /* i=254 */
   {0x1.61da5d0052b67p+1, 0x1.f9d8d5c2a8d6fp+2, 0x1.fde12d0052237p+2}, /* i=255 */
};

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

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma(a, b, -*hi);
}

// Multiply a double with a double double : a * (bh + bl)
static inline void
s_mul (double *hi, double *lo, double a, double bh, double bl)
{
  double s;

  a_mul (hi, &s, a, bh); /* exact */
  *lo = __builtin_fma (a, bl, s);
  /* the error is bounded by ulp(lo), where |lo| < |a*bl| + ulp(hi) */
}

/* Put in hi+lo an approximation of (ah+al)*(bh+bl) */
static inline void
d_mul (double *hi, double *lo, double ah, double al, double bh, double bl)
{
  double s, t;

  a_mul (hi, &s, ah, bh); /* exact */
  t = __builtin_fma (al, bh, s);
  *lo = __builtin_fma (ah, bl, t);
}

/* the following is a degree-7 odd polynomial approximating sinh(x)
   for |x| < 0.00543 generated by Sollya (see Psinh.sollya), with
   maximal relative error 2^-74.818 and maximal absolute error 2^-83.263 */
static const double S[] = {
  1,                     /* degree 1 */
  0x1.5555555555555p-3,  /* degree 3 */
  0x1.11111111869d4p-7,  /* degree 5 */
  0x1.a01061b363a81p-13, /* degree 7 */
};

/* the following is a degree-6 even polynomial approximating cosh(x)
   for |x| < 0.00543 generated by Sollya (see Pcosh.sollya), with
   maximal relative/absolute error 2^-81.152 */
static const double C[] = {
  1,                     /* degree 0 */
  0x1p-1,                /* degree 2 */
  0x1.5555555554e2ep-5,  /* degree 4 */
  0x1.6c16d52a52a35p-10, /* degree 6 */
};

/* put in h+l a double-double approximation of sinh(w), for |w| < 0.00543 */
static void
eval_S (double *h, double *l, double w)
{
  double z = w * w;
  *h = S[3];
  *h = __builtin_fma (*h, z, S[2]);
  *h = __builtin_fma (*h, z, S[1]);
  *h = __builtin_fma (*h, z, S[0]);
  *l = 0;
  /* multiply by w */
  *h *= w;
}

/* put in h+l a double-double approximation of cosh(w), for |w| < 0.00543 */
static void
eval_C (double *h, double *l, double w)
{
  double z = w * w;
  *h = C[3];
  *h = __builtin_fma (*h, z, C[2]);
  *h = __builtin_fma (*h, z, C[1]);
  *h = __builtin_fma (*h, z, C[0]);
  *l = 0;
}

/* put in h+l a double-double approximation of sinh(x),
   for 0 <= x <= 0x1.633ce8fb9f87dp+9 */
static void
cr_sinh_fast (double *h, double *l, double x)
{
  /* magic is such that magic*x rounds to a number < 65535.5
     whatever the rounding mode */
  static const magic = 0x1.70f77fc88ae3cp6;
  int k = __builtin_round (magic * x); /* k <= 65535 */
  int i = k >> 8, j = k & 0xff;
  double v = x - T[i][0];
  /* since x = T[i][0] + v, we approximate sinh(x) as
     sinh(T[i][0])*cosh(v) + cosh(T[i][0])*sinh(v)
     = T[i][1]*cosh(v) + T[i][2]*sinh(v) */
  double w = v - U[j][0];
  /* since v = U[j][0] + w, we approximate sinh(v) as
     sinh(U[j][0])*cosh(w) + cosh(U[j][0])*sinh(w)
     = U[j][1]*cosh(w) + U[j][2]*sinh(w), and cosh(v) as
     sinh(U[j][0])*sinh(w) + cosh(U[j][0])*cosh(w)
     = U[j][1]*sinh(w) + U[j][2]*cosh(w) */

  /* we have |w| < 0.00543 */
  double swh, swl, cwh, cwl;
  eval_S (&swh, &swl, w);
  eval_C (&cwh, &cwl, w);

  double svh, svl, cvh, cvl, h1, l1, h2, l2;
  s_mul (&h1, &l1, U[j][1], cwh, cwl); /* U[j][1]*cosh(w) */
  s_mul (&h2, &l2, U[j][2], swh, swl); /* U[j][1]*sinh(w) */
  fast_two_sum (&svh, &svl, h1, h2);
  svl += l1 + l2; /* svh+svl approximates sinh(v) */
  s_mul (&h1, &l1, U[j][1], swh, swl); /* U[j][1]*sinh(w) */
  s_mul (&h2, &l2, U[j][2], cwh, cwl); /* U[j][2]*cosh(w) */
  fast_two_sum (&cvh, &cvl, h1, h2);
  cvl += l1 + l2; /* cvh+cvl approximates cosh(v) */

  s_mul (&h1, &l1, T[i][1], cvh, cvl); /* T[i][1]*cosh(v) */
  s_mul (&h2, &l2, T[i][2], svh, svl); /* T[i][2]*sinh(v) */
  fast_two_sum (h, l, h1, h2);
  l += l1 + l2;
}

#define MASK 0x7fffffffffffffff /* to mask the sign bit */

double
cr_sinh (double x)
{
  d64u64 v = {.f = x};
  int s = v >> 63; /* sign bit */
  v.u &= (uint64_t) MASK;
  int e = (v.u >> 52) - 0x3ff;
  /* 2^(e-1) <= |x| < 2^e */

  if (v.f >= 0x1.633ce8fb9f87ep+9) /* overflow */
  {
    /* this will return NaN for x=NaN, Inf with the correct sign for x=+Inf
       or -Inf, and Inf/DBL_MAX or -Inf/-DBL_MAX for other |x| >= 2^10 */
    return x * 0x1p1023;
  }
  
  double h, l;
  cr_sinh_fast (&h, &l, v.f);
  double sign[] = { 1.0, -1.0 };
  h *= sign[s];
  l *= sign[s];
  return h + l;
}
