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

#define TRACE -0x1.21a0d7b7f3335p-1010

#include <stdio.h>
#include <stdint.h>

/* For 0 <= i < 256, T[i] = {xi, si, ei} such that xi is near i*2^8/magic
   with magic = 0x1.70f77fc88ae3cp6, and si,si+ei approximate sinh(xi),cosh(xi)
   with accuracy >= 53+16 bits:
   |si - sinh(xi)| < 2^(-16-1) ulp(si), |ci - (si+ei)| < 2^(-16-1) ulp(si+ei).
   We have |xi - i*2^8/magic| < 2.36e-8.
   Generated with build_table_T(k=16) from accompanying file sinh.sage.
*/
static const double T[256][3] = {
   {0x0p+0, 0x0p+0, 0x1p+0}, /* i=0 */
   {0x1.633d9a9077741p+1, 0x1.ff678cbb5f806p+2, 0x1.fe9ad24cdfbp-5}, /* i=1 */
   {0x1.633d9a9a65199p+2, 0x1.0165a65ef5742p+7, 0x1.fd369d76ffbf2p-9}, /* i=2 */
   {0x1.0a6e33f3d334fp+3, 0x1.021ab2b881192p+11, 0x1.fbd36146e7db4p-13}, /* i=3 */
   {0x1.633d9a9a5efefp+3, 0x1.02cf40653010dp+15, 0x1.fa711ce962a12p-17}, /* i=4 */
   {0x1.bc0d0140f44d1p+3, 0x1.03844b629b49bp+19, 0x1.f90fcfafac308p-21}, /* i=5 */
   {0x1.0a6e33f3d3c0dp+4, 0x1.0439d504b756fp+23, 0x1.f7af78ec1be86p-25}, /* i=6 */
   {0x1.36d5e74711931p+4, 0x1.04efdda25e99p+27, 0x1.f65017f8026e2p-29}, /* i=7 */
   {0x1.633d9a9a6f0cep+4, 0x1.05a665982f42bp+31, 0x1.f4f1ac209af2ep-33}, /* i=8 */
   {0x1.8fa54dedb8507p+4, 0x1.065d6d3bee2fp+35, 0x1.f39434c13fdeap-37}, /* i=9 */
   {0x1.bc0d01410e054p+4, 0x1.0714f4e902ff7p+39, 0x1.f237b12b69573p-41}, /* i=10 */
   {0x1.e874b4944a5p+4, 0x1.07ccfcf68f1d5p+43, 0x1.f0dc20b99f389p-45}, /* i=11 */
   {0x1.0a6e33f3d8f7ep+5, 0x1.088585c2ce353p+47, 0x1.ef8182b9edc86p-49}, /* i=12 */
   {0x1.20a20d9d69758p+5, 0x1.093e8fa072122p+51, 0x1.ee27d690fa264p-53}, /* i=13 */
   {0x1.36d5e7471f37p+5, 0x1.09f81af32ab2ap+55, 0x1.eccf1b848b20ap-57}, /* i=14 */
   {0x1.4d09c0f0ca428p+5, 0x1.0ab2280f4e61bp+59, 0x1.eb7750f7fc962p-61}, /* i=15 */
   {0x1.633d9a9a660d4p+5, 0x1.0b6cb74f0d3efp+63, 0x1.ea2076449a7fap-65}, /* i=16 */
   {0x1.797174440f35p+5, 0x1.0c27c9112cb1ap+67, 0x1.e8ca8abc38f11p-69}, /* i=17 */
   {0x1.8fa54dedab4bap+5, 0x1.0ce35dad821c2p+71, 0x1.e7758dbe4debbp-73}, /* i=18 */
   {0x1.a5d9279764a4ep+5, 0x1.0d9f7585265dbp+75, 0x1.e6217e9a601e9p-77}, /* i=19 */
   {0x1.bc0d014104f08p+5, 0x1.0e5c10ecbff98p+79, 0x1.e4ce5cb766df7p-81}, /* i=20 */
   {0x1.d240daeac98d2p+5, 0x1.0f19304871d39p+83, 0x1.e37c276159f25p-85}, /* i=21 */
   {0x1.e874b49457c96p+5, 0x1.0fd6d3e886d21p+87, 0x1.e22ade089a00dp-89}, /* i=22 */
   {0x1.fea88e3dee0c7p+5, 0x1.1094fc31c92fep+91, 0x1.e0da7ff9cd127p-93}, /* i=23 */
   {0x1.0a6e33f3c8752p+6, 0x1.1153a981a366fp+95, 0x1.df8b0c8fbf3f3p-97}, /* i=24 */
   {0x1.158820c8916e5p+6, 0x1.1212dc3132448p+99, 0x1.de3c832da8b54p-101}, /* i=25 */
   {0x1.20a20d9d75dabp+6, 0x1.12d294a762155p+103, 0x1.dceee31f86c9bp-105}, /* i=26 */
   {0x1.2bbbfa722e1f4p+6, 0x1.1392d32e8c449p+107, 0x1.dba22be3e7e32p-109}, /* i=27 */
   {0x1.36d5e74724076p+6, 0x1.14539840f4ed5p+111, 0x1.da565ca7422d8p-113}, /* i=28 */
   {0x1.41efd41bf286ep+6, 0x1.1514e42179819p+115, 0x1.d90b74f693782p-117}, /* i=29 */
   {0x1.4d09c0f0c42cdp+6, 0x1.15d6b739e3cc5p+119, 0x1.d7c1741c9f285p-121}, /* i=30 */
   {0x1.5823adc5b1a95p+6, 0x1.169911ef8312dp+123, 0x1.d678596d06228p-125}, /* i=31 */
   {0x1.633d9a9a6e96ep+6, 0x1.175bf48c67da3p+127, 0x1.d530246a4a2ep-129}, /* i=32 */
   {0x1.6e57876f3b02cp+6, 0x1.181f5f8116732p+131, 0x1.d3e8d456ce04fp-133}, /* i=33 */
   {0x1.79717444087aep+6, 0x1.18e35328fe587p+135, 0x1.d2a268997b6cep-137}, /* i=34 */
   {0x1.848b6118f520fp+6, 0x1.19a7cfec098ap+139, 0x1.d15ce0855703p-141}, /* i=35 */
   {0x1.8fa54dedb6ee6p+6, 0x1.1a6cd615c1f5ap+143, 0x1.d0183b9d2f6d1p-145}, /* i=36 */
   {0x1.9abf3ac293105p+6, 0x1.1b3266195629bp+147, 0x1.ced479232ad2dp-149}, /* i=37 */
   {0x1.a5d9279760ec8p+6, 0x1.1bf8804bfbddep+151, 0x1.cd91988bb0f73p-153}, /* i=38 */
   {0x1.b0f3146c33ea5p+6, 0x1.1cbf2513bc758p+155, 0x1.cc4f99306b3b6p-157}, /* i=39 */
   {0x1.bc0d0140fcb9cp+6, 0x1.1d8654cd45e9fp+159, 0x1.cb0e7a7b1985ap-161}, /* i=40 */
   {0x1.c726ee15db1fcp+6, 0x1.1e4e0fe2a8059p+163, 0x1.c9ce3bc0c9986p-165}, /* i=41 */
   {0x1.d240daeab5408p+6, 0x1.1f1656ae27543p+167, 0x1.c88edc70cfd77p-169}, /* i=42 */
   {0x1.dd5ac7bf84a89p+6, 0x1.1fdf298fad145p+171, 0x1.c7505bf22c3d7p-173}, /* i=43 */
   {0x1.e874b4945270ap+6, 0x1.20a888ebc3c4p+175, 0x1.c612b9a55e9a7p-177}, /* i=44 */
   {0x1.f38ea169333c2p+6, 0x1.2172752a84c61p+179, 0x1.c4d5f4e643d37p-181}, /* i=45 */
   {0x1.fea88e3df7e1p+6, 0x1.223ceea126d0fp+183, 0x1.c39a0d2f255f1p-185}, /* i=46 */
   {0x1.04e13d896bc0ep+7, 0x1.2307f5c2133e9p+187, 0x1.c25f01cd80441p-189}, /* i=47 */
   {0x1.0a6e33f3ce1bp+7, 0x1.23d38ae110ea2p+191, 0x1.c124d23f54d41p-193}, /* i=48 */
   {0x1.0ffb2a5e2bd66p+7, 0x1.249fae66766b7p+195, 0x1.bfeb7de37bb8ep-197}, /* i=49 */
   {0x1.158820c8a300ep+7, 0x1.256c60c71601fp+199, 0x1.beb30406cc819p-201}, /* i=50 */
   {0x1.1b15173311defp+7, 0x1.2639a2538fedp+203, 0x1.bd7b642e5292bp-205}, /* i=51 */
   {0x1.20a20d9d796d3p+7, 0x1.2707737095da9p+207, 0x1.bc449dc111294p-209}, /* i=52 */
   {0x1.262f0407db5bep+7, 0x1.27d5d48388dcdp+211, 0x1.bb0eb025dca16p-213}, /* i=53 */
   {0x1.2bbbfa7251d8ap+7, 0x1.28a4c6004353cp+215, 0x1.b9d99aaed4e44p-217}, /* i=54 */
   {0x1.3148f0dcb5cbap+7, 0x1.2974483524df3p+219, 0x1.b8a55ce6ef373p-221}, /* i=55 */
   {0x1.36d5e7471f2bcp+7, 0x1.2a445b9550519p+223, 0x1.b771f6230ce14p-225}, /* i=56 */
   {0x1.3c62ddb162742p+7, 0x1.2b15006ceaba9p+227, 0x1.b63f65f274e58p-229}, /* i=57 */
   {0x1.41efd41bf5431p+7, 0x1.2be637667b99ep+231, 0x1.b50dab5ae217dp-233}, /* i=58 */
   {0x1.477cca866558ep+7, 0x1.2cb800a562c48p+235, 0x1.b3dcc628f34b6p-237}, /* i=59 */
   {0x1.4d09c0f0d62d6p+7, 0x1.2d8a5ca4b916ep+239, 0x1.b2acb5a9842f7p-241}, /* i=60 */
   {0x1.5296b75b2e80cp+7, 0x1.2e5d4bbc3de1ep+243, 0x1.b17d795d94c81p-245}, /* i=61 */
   {0x1.5823adc59f02ep+7, 0x1.2f30ce6f96f2dp+247, 0x1.b04f1087e34c7p-249}, /* i=62 */
   {0x1.5db0a430012bfp+7, 0x1.3004e50f38b3ep+251, 0x1.af217ab57937p-253}, /* i=63 */
   {0x1.633d9a9a69ca5p+7, 0x1.30d9900ef44bdp+255, 0x1.adf4b741981d2p-257}, /* i=64 */
   {0x1.68ca9104ce404p+7, 0x1.31aecfd04334ep+259, 0x1.acc8c5a269d97p-261}, /* i=65 */
   {0x1.6e57876f46315p+7, 0x1.3284a4c957bap+263, 0x1.ab9da531c9c73p-265}, /* i=66 */
   {0x1.73e47dd9a70cep+7, 0x1.335b0f491672ap+267, 0x1.aa73558154d58p-269}, /* i=67 */
   {0x1.797174441a449p+7, 0x1.34320fd0f7294p+271, 0x1.a949d5dce2af9p-273}, /* i=68 */
   {0x1.7efe6aae85e4bp+7, 0x1.3509a6ba52952p+275, 0x1.a82125c8e51a2p-277}, /* i=69 */
   {0x1.848b6118e4679p+7, 0x1.35e1d46afc842p+279, 0x1.a6f944b9165b1p-281}, /* i=70 */
   {0x1.8a1857833affdp+7, 0x1.36ba994f8a757p+283, 0x1.a5d23218b41c3p-285}, /* i=71 */
   {0x1.8fa54dedc7a43p+7, 0x1.3793f5f76ff94p+287, 0x1.a4abed24c8e1fp-289}, /* i=72 */
   {0x1.953244581c997p+7, 0x1.386dea8a0ebfap+291, 0x1.a38675a8010adp-293}, /* i=73 */
   {0x1.9abf3ac293e81p+7, 0x1.394877a893217p+295, 0x1.a261cac934c41p-297}, /* i=74 */
   {0x1.a04c312cf16eap+7, 0x1.3a239d98eb75ep+299, 0x1.a13dec2acd45bp-301}, /* i=75 */
   {0x1.a5d927975fb78p+7, 0x1.3aff5ce014458p+303, 0x1.a01ad91b9d5fap-305}, /* i=76 */
   {0x1.ab661e01bd3edp+7, 0x1.3bdbb5d4b28dap+307, 0x1.9ef89128e067ap-309}, /* i=77 */
   {0x1.b0f3146c3d309p+7, 0x1.3cb8a901dfaa7p+311, 0x1.9dd7139b75655p-313}, /* i=78 */
   {0x1.b6800ad69a0e2p+7, 0x1.3d9636a87f998p+315, 0x1.9cb6601e560d1p-317}, /* i=79 */
   {0x1.bc0d014113f9cp+7, 0x1.3e745f5c66b04p+319, 0x1.9b9675f0d615bp-321}, /* i=80 */
   {0x1.c199f7ab7ae41p+7, 0x1.3f53236c2b206p+323, 0x1.9a7754ad2d91fp-325}, /* i=81 */
   {0x1.c726ee15daf59p+7, 0x1.4032834c0425p+327, 0x1.9958fbbd65b72p-329}, /* i=82 */
   {0x1.ccb3e4803620cp+7, 0x1.41127f6a22eefp+331, 0x1.983b6a9428fe1p-333}, /* i=83 */
   {0x1.d240daeaab289p+7, 0x1.41f3184726aa3p+335, 0x1.971ea08d9e55ep-337}, /* i=84 */
   {0x1.d7cdd1551632bp+7, 0x1.42d44e3a308abp+339, 0x1.96029d3b653adp-341}, /* i=85 */
   {0x1.dd5ac7bf8a4a2p+7, 0x1.43b621bd248f9p+343, 0x1.94e76003c4cbp-345}, /* i=86 */
   {0x1.e2e7be29f16a5p+7, 0x1.449893304d214p+347, 0x1.93cce86df1207p-349}, /* i=87 */
   {0x1.e874b49461158p+7, 0x1.457ba30fc656cp+351, 0x1.92b335df1222bp-353}, /* i=88 */
   {0x1.ee01aafeb7465p+7, 0x1.465f51b4ba0cap+355, 0x1.919a47e86da9ep-357}, /* i=89 */
   {0x1.f38ea16920f8fp+7, 0x1.47439faae2f2bp+359, 0x1.90821ddd919d4p-361}, /* i=90 */
   {0x1.f91b97d38f1bep+7, 0x1.48288d581c307p+363, 0x1.8f6ab741a722cp-365}, /* i=91 */
   {0x1.fea88e3dfb5c1p+7, 0x1.490e1b281297fp+367, 0x1.8e541391484e1p-369}, /* i=92 */
   {0x1.021ac2542b57ap+8, 0x1.49f4498117eb1p+371, 0x1.8d3e325023a62p-373}, /* i=93 */
   {0x1.04e13d89669afp+8, 0x1.4adb18efecfbp+375, 0x1.8c2912d45e3e8p-377}, /* i=94 */
   {0x1.07a7b8be9cbp+8, 0x1.4bc289cd0245ep+379, 0x1.8b14b4b3e4113p-381}, /* i=95 */
   {0x1.0a6e33f3c7e49p+8, 0x1.4caa9c81d4106p+383, 0x1.8a0117708ad51p-385}, /* i=96 */
   {0x1.0d34af2904c3p+8, 0x1.4d9351a4b5669p+387, 0x1.88ee3a57e1467p-389}, /* i=97 */
   {0x1.0ffb2a5e3be0dp+8, 0x1.4e7ca988cda97p+391, 0x1.87dc1d07d6127p-393}, /* i=98 */
   {0x1.12c1a593675abp+8, 0x1.4f66a4983541p+395, 0x1.86cabf038f504p-397}, /* i=99 */
   {0x1.158820c8a4297p+8, 0x1.5051436b0a8c1p+399, 0x1.85ba1f9965986p-401}, /* i=100 */
   {0x1.184e9bfdd73bdp+8, 0x1.513c8650539d3p+403, 0x1.84aa3e6d9e5fdp-405}, /* i=101 */
   {0x1.1b1517330fbap+8, 0x1.52286dcec99ccp+407, 0x1.839b1ae47a883p-409}, /* i=102 */
   {0x1.1ddb92683cca3p+8, 0x1.5314fa43489cep+411, 0x1.828cb4932c526p-413}, /* i=103 */
   {0x1.20a20d9d62b2ap+8, 0x1.54022c26cd9bap+415, 0x1.817f0aef2b816p-417}, /* i=104 */
   {0x1.236888d2ad2e5p+8, 0x1.54f00427421a2p+419, 0x1.80721d331d57dp-421}, /* i=105 */
   {0x1.262f0407e1ef2p+8, 0x1.55de826b40941p+423, 0x1.7f65eb338a8c8p-425}, /* i=106 */
   {0x1.28f57f3d176d2p+8, 0x1.56cda785050f2p+427, 0x1.7e5a744bbed31p-429}, /* i=107 */
   {0x1.2bbbfa72450e5p+8, 0x1.57bd73ddb8a4p+431, 0x1.7d4fb80606bfbp-433}, /* i=108 */
   {0x1.2e8275a78cb2cp+8, 0x1.58ade817e6a06p+435, 0x1.7c45b5adcd7p-437}, /* i=109 */
   {0x1.3148f0dcb4446p+8, 0x1.599f045ac3dcfp+439, 0x1.7b3c6d17a0728p-441}, /* i=110 */
   {0x1.340f6c11efe59p+8, 0x1.5a90c96232de3p+443, 0x1.7a33dd74a1acap-445}, /* i=111 */
   {0x1.36d5e7473131fp+8, 0x1.5b833790deb43p+447, 0x1.792c065928ac8p-449}, /* i=112 */
   {0x1.399c627c55ep+8, 0x1.5c764f2e79f3p+451, 0x1.7824e776e6f41p-453}, /* i=113 */
   {0x1.3c62ddb191ff5p+8, 0x1.5d6a10f862a1bp+455, 0x1.771e8000d46aap-457}, /* i=114 */
   {0x1.3f2958e6b4791p+8, 0x1.5e5e7d22a0575p+459, 0x1.7618cfbee7c29p-461}, /* i=115 */
   {0x1.41efd41bffd3ep+8, 0x1.5f53947f7cddp+463, 0x1.7513d5d02b18cp-465}, /* i=116 */
   {0x1.44b64f512f715p+8, 0x1.604957289acdp+467, 0x1.740f92197aab1p-469}, /* i=117 */
   {0x1.477cca865bc1ep+8, 0x1.613fc5b751a06p+471, 0x1.730c03f813941p-473}, /* i=118 */
   {0x1.4a4345bb9179fp+8, 0x1.6236e0b56df96p+475, 0x1.72092adae528dp-477}, /* i=119 */
   {0x1.4d09c0f0be881p+8, 0x1.632ea8828fd0fp+479, 0x1.7107065dc6604p-481}, /* i=120 */
   {0x1.4fd03c25fdf46p+8, 0x1.64271dbd23ec4p+483, 0x1.700595dbc7bfcp-485}, /* i=121 */
   {0x1.5296b75b3a49ap+8, 0x1.652040c0afca3p+487, 0x1.6f04d8f642eaap-489}, /* i=122 */
   {0x1.555d32906a0c3p+8, 0x1.661a11f97a90cp+491, 0x1.6e04cf3d73e68p-493}, /* i=123 */
   {0x1.5823adc594d31p+8, 0x1.671491ebedc4fp+495, 0x1.6d057829706c1p-497}, /* i=124 */
   {0x1.5aea28fac3111p+8, 0x1.680fc11e1af36p+499, 0x1.6c06d3318b33ep-501}, /* i=125 */
   {0x1.5db0a43005669p+8, 0x1.690ba02212f6p+503, 0x1.6b08dfc1f69a9p-505}, /* i=126 */
   {0x1.60771f6539733p+8, 0x1.6a082f425d5f9p+507, 0x1.6a0b9d8f9ecap-509}, /* i=127 */
   {0x1.633d9a9a69e5ep+8, 0x1.6b056f0935b4dp+511, 0x1.690f0c0fb1258p-513}, /* i=128 */
   {0x1.660415cfa1401p+8, 0x1.6c03600117277p+515, 0x1.68132ab81b8e4p-517}, /* i=129 */
   {0x1.68ca9104c92edp+8, 0x1.6d0202862e0afp+519, 0x1.6717f92d58117p-521}, /* i=130 */
   {0x1.6b910c39fd65cp+8, 0x1.6e01573c2a145p+523, 0x1.661d76cde40ecp-525}, /* i=131 */
   {0x1.6e57876f44364p+8, 0x1.6f015ea8d87f1p+527, 0x1.6523a316d3c88p-529}, /* i=132 */
   {0x1.711e02a47361fp+8, 0x1.7002190ca8316p+531, 0x1.642a7dc93925ep-533}, /* i=133 */
   {0x1.73e47dd9a8ccp+8, 0x1.7103870faf39dp+535, 0x1.63320641c22b3p-537}, /* i=134 */
   {0x1.76aaf90ed015fp+8, 0x1.7205a9122f658p+539, 0x1.623a3c237009cp-541}, /* i=135 */
   {0x1.797174440e50dp+8, 0x1.73087fc7aeb12p+543, 0x1.61431ec207e4cp-545}, /* i=136 */
   {0x1.7c37ef79511d8p+8, 0x1.740c0b940facp+547, 0x1.604cadbe7d0bfp-549}, /* i=137 */
   {0x1.7efe6aae77193p+8, 0x1.75104cc561995p+551, 0x1.5f56e8ce6b2dp-553}, /* i=138 */
   {0x1.81c4e5e3b13acp+8, 0x1.76154421da255p+555, 0x1.5e61cf368bd84p-557}, /* i=139 */
   {0x1.848b6118e828cp+8, 0x1.771af206da4abp+559, 0x1.5d6d609f5c213p-561}, /* i=140 */
   {0x1.8751dc4e17d2ep+8, 0x1.782156ee129f2p+563, 0x1.5c799c9712bb5p-565}, /* i=141 */
   {0x1.8a1857835936p+8, 0x1.7928737c494c9p+567, 0x1.5b868284ba0a1p-569}, /* i=142 */
   {0x1.8cded2b8886ddp+8, 0x1.7a3047fd14f6p+571, 0x1.5a94122279c3cp-573}, /* i=143 */
   {0x1.8fa54dedc53b2p+8, 0x1.7b38d520067a1p+575, 0x1.59a24acef0194p-577}, /* i=144 */
   {0x1.926bc92300d5fp+8, 0x1.7c421b506cfcdp+579, 0x1.58b12c282db95p-581}, /* i=145 */
   {0x1.9532445835dbdp+8, 0x1.7d4c1b07b3af9p+583, 0x1.57c0b5bfbfbc6p-585}, /* i=146 */
   {0x1.97f8bf8d4e408p+8, 0x1.7e56d4a6a6d09p+587, 0x1.56d0e73dc19d9p-589}, /* i=147 */
   {0x1.9abf3ac28734p+8, 0x1.7f62490ab8443p+591, 0x1.55e1bfdb25c35p-593}, /* i=148 */
   {0x1.9d85b5f7cd9fp+8, 0x1.806e789a27483p+595, 0x1.54f33f3cffaaep-597}, /* i=149 */
   {0x1.a04c312cf523ap+8, 0x1.817b63952f70ep+599, 0x1.54056529ff08fp-601}, /* i=150 */
   {0x1.a312ac6226e77p+8, 0x1.82890abccb21ep+603, 0x1.531830f73afdep-605}, /* i=151 */
   {0x1.a5d927975c685p+8, 0x1.83976e8ad72ap+607, 0x1.522ba239a7b71p-609}, /* i=152 */
   {0x1.a89fa2cc98b87p+8, 0x1.84a68f87fc04ap+611, 0x1.513fb879d35d9p-613}, /* i=153 */
   {0x1.ab661e01c1dccp+8, 0x1.85b66e1111c07p+615, 0x1.50547366d2d7cp-617}, /* i=154 */
   {0x1.ae2c9936ff4bdp+8, 0x1.86c70ae6c4697p+619, 0x1.4f69d259ed44ep-621}, /* i=155 */
   {0x1.b0f3146c56cebp+8, 0x1.87d86697389c5p+623, 0x1.4e7fd4d94b9bcp-625}, /* i=156 */
   {0x1.b3b98fa16ad76p+8, 0x1.88ea811890d7p+627, 0x1.4d967aecee226p-629}, /* i=157 */
   {0x1.b6800ad6a1ad3p+8, 0x1.89fd5b8d1fddbp+631, 0x1.4cadc39d52a48p-633}, /* i=158 */
   {0x1.b946860bc2896p+8, 0x1.8b10f623d80bap+635, 0x1.4bc5aec2dc44ep-637}, /* i=159 */
   {0x1.bc0d014107222p+8, 0x1.8c2551bc3ef5ep+639, 0x1.4ade3ba177178p-641}, /* i=160 */
   {0x1.bed37c7645439p+8, 0x1.8d3a6e9c28fc6p+643, 0x1.49f769fef7085p-645}, /* i=161 */
   {0x1.c199f7ab6ecd7p+8, 0x1.8e504d34c18dfp+647, 0x1.4911397ccf60dp-649}, /* i=162 */
   {0x1.c46072e0b3b59p+8, 0x1.8f66ee58326eap+651, 0x1.482ba96cf8907p-653}, /* i=163 */
   {0x1.c726ee15d4135p+8, 0x1.907e522ad27b7p+655, 0x1.4746b9b16949ap-657}, /* i=164 */
   {0x1.c9ed694b152f2p+8, 0x1.919679a1308eep+659, 0x1.46626981b2547p-661}, /* i=165 */
   {0x1.ccb3e4802fc96p+8, 0x1.92af64d45c65ap+663, 0x1.457eb8c957ba5p-665}, /* i=166 */
   {0x1.cf7a5fb571fb4p+8, 0x1.93c914c81c7cap+667, 0x1.449ba6b5cec5bp-669}, /* i=167 */
   {0x1.d240daeab5f89p+8, 0x1.94e389caa233ep+671, 0x1.43b933087d695p-673}, /* i=168 */
   {0x1.d507561fd3c41p+8, 0x1.95fec4265eb8dp+675, 0x1.42d75d855481p-677}, /* i=169 */
   {0x1.d7cdd155087f8p+8, 0x1.971ac4c64114cp+679, 0x1.41f62570d8267p-681}, /* i=170 */
   {0x1.da944c8a53fb3p+8, 0x1.98378c34e13ap+683, 0x1.41158a5d98b25p-685}, /* i=171 */
   {0x1.dd5ac7bf76ab3p+8, 0x1.99551a97e8315p+687, 0x1.40358c2db2bfp-689}, /* i=172 */
   {0x1.e02142f4b1afdp+8, 0x1.9a7370e2923a3p+691, 0x1.3f562a222f4a4p-693}, /* i=173 */
   {0x1.e2e7be29ed5fbp+8, 0x1.9b928f7ae1d91p+695, 0x1.3e7763ebbacb3p-697}, /* i=174 */
   {0x1.e5ae395f1d65p+8, 0x1.9cb276d913b39p+699, 0x1.3d99392cee083p-701}, /* i=175 */
   {0x1.e874b49444378p+8, 0x1.9dd3278d7ca88p+703, 0x1.3cbba97632507p-705}, /* i=176 */
   {0x1.eb3b2fc9955fdp+8, 0x1.9ef4a27866fe3p+707, 0x1.3bdeb41bbad2p-709}, /* i=177 */
   {0x1.ee01aafebaa3ep+8, 0x1.a016e79b80653p+711, 0x1.3b02591c3b8a5p-713}, /* i=178 */
   {0x1.f0c82633fa2ep+8, 0x1.a139f7f634ed1p+715, 0x1.3a2697b5a0529p-717}, /* i=179 */
   {0x1.f38ea1692528bp+8, 0x1.a25dd3ca5e50dp+719, 0x1.394b6fb6402a1p-721}, /* i=180 */
   {0x1.f6551c9e61309p+8, 0x1.a3827be3c5d1p+723, 0x1.3870e0850c3edp-725}, /* i=181 */
   {0x1.f91b97d385c43p+8, 0x1.a4a7f08ee7d16p+727, 0x1.3796e9e8d54e8p-729}, /* i=182 */
   {0x1.fbe21308bf692p+8, 0x1.a5ce32a403b1cp+731, 0x1.36bd8b40f7409p-733}, /* i=183 */
   {0x1.fea88e3df7c15p+8, 0x1.a6f5428e06818p+735, 0x1.35e4c43ea7d6ep-737}, /* i=184 */
   {0x1.00b784b9a1438p+9, 0x1.a81d20fd97d03p+739, 0x1.350c946049907p-741}, /* i=185 */
   {0x1.021ac2543dfd3p+9, 0x1.a945ce47f68e1p+743, 0x1.3434fb6795889p-745}, /* i=186 */
   {0x1.037dffeed1c99p+9, 0x1.aa6f4afcca835p+747, 0x1.335df8ebcd9ccp-749}, /* i=187 */
   {0x1.04e13d89795e1p+9, 0x1.ab99980ce954bp+751, 0x1.32878c3f18264p-753}, /* i=188 */
   {0x1.06447b23e6da7p+9, 0x1.acc4b50591e5fp+755, 0x1.31b1b5b3316f1p-757}, /* i=189 */
   {0x1.07a7b8be88971p+9, 0x1.adf0a3e9c7685p+759, 0x1.30dc73d7cc921p-761}, /* i=190 */
   {0x1.090af65935ee7p+9, 0x1.af1d64c4707d9p+763, 0x1.3007c6a6688a2p-765}, /* i=191 */
   {0x1.0a6e33f3cf0e9p+9, 0x1.b04af7bd2c18cp+767, 0x1.2f33ae02d5debp-769}, /* i=192 */
   {0x1.0bd1718e63208p+9, 0x1.b1795d99ed6efp+771, 0x1.2e6029615f33fp-773}, /* i=193 */
   {0x1.0d34af2902bafp+9, 0x1.b2a89726771b9p+775, 0x1.2d8d38339e80ap-777}, /* i=194 */
   {0x1.0e97ecc39d1dap+9, 0x1.b3d8a4be0a7f8p+779, 0x1.2cbada3a375bcp-781}, /* i=195 */
   {0x1.0ffb2a5e3423cp+9, 0x1.b50986fb39673p+783, 0x1.2be90f0a0e12p-785}, /* i=196 */
   {0x1.115e67f8cf85ap+9, 0x1.b63b3e8d22521p+787, 0x1.2b17d62aa7342p-789}, /* i=197 */
   {0x1.12c1a59373a55p+9, 0x1.b76dcc1820f9dp+791, 0x1.2a472f2bc7d41p-793}, /* i=198 */
   {0x1.1424e32df90a1p+9, 0x1.b8a12faa2c562p+795, 0x1.29771a0383254p-797}, /* i=199 */
   {0x1.158820c8a00bap+9, 0x1.b9d56ab68199ep+799, 0x1.28a795b6991c5p-801}, /* i=200 */
   {0x1.16eb5e632e232p+9, 0x1.bb0a7d09ea956p+803, 0x1.27d8a267ae7bdp-805}, /* i=201 */
   {0x1.184e9bfdc6899p+9, 0x1.bc4067b4c58e9p+807, 0x1.270a3f601c96cp-809}, /* i=202 */
   {0x1.19b1d99876159p+9, 0x1.bd772b7b4886ap+811, 0x1.263c6c1ddaefap-813}, /* i=203 */
   {0x1.1b15173310cffp+9, 0x1.beaec85c56067p+815, 0x1.256f28a1b1d2fp-817}, /* i=204 */
   {0x1.1c7854cda66f8p+9, 0x1.bfe73f2663c0dp+819, 0x1.24a2746325159p-821}, /* i=205 */
   {0x1.1ddb926848ea3p+9, 0x1.c12090b0ccdp+823, 0x1.23d64ed543f32p-825}, /* i=206 */
   {0x1.1f3ed002ab6b5p+9, 0x1.c25abc86838d3p+827, 0x1.230ab8435e567p-829}, /* i=207 */
   {0x1.20a20d9d78e83p+9, 0x1.c395c59a2513fp+831, 0x1.223faec540e9p-833}, /* i=208 */
   {0x1.22054b3817321p+9, 0x1.c4d1aa6712e23p+835, 0x1.217533567e2e6p-837}, /* i=209 */
   {0x1.236888d2a282bp+9, 0x1.c60e6bea37535p+839, 0x1.20ab4553bee02p-841}, /* i=210 */
   {0x1.24cbc66d3f314p+9, 0x1.c74c0b3edc4e3p+843, 0x1.1fe1e408632bdp-845}, /* i=211 */
   {0x1.262f0407df163p+9, 0x1.c88a88ce0693ap+847, 0x1.1f190f3236a28p-849}, /* i=212 */
   {0x1.279241a278889p+9, 0x1.c9c9e510adf14p+851, 0x1.1e50c684e2966p-853}, /* i=213 */
   {0x1.28f57f3d0e0b9p+9, 0x1.cb0a20ab7d94fp+855, 0x1.1d890998f1103p-857}, /* i=214 */
   {0x1.2a58bcd7ba70dp+9, 0x1.cc4b3c9b0ecfcp+859, 0x1.1cc1d7d12d3ep-861}, /* i=215 */
   {0x1.2bbbfa725bcb9p+9, 0x1.cd8d39026d4a2p+863, 0x1.1bfb311811006p-865}, /* i=216 */
   {0x1.2d1f380ce9114p+9, 0x1.ced0165db5cb6p+867, 0x1.1b3515207904dp-869}, /* i=217 */
   {0x1.2e8275a77cdf1p+9, 0x1.d013d5aa5cb63p+871, 0x1.1a6f834ea16bcp-873}, /* i=218 */
   {0x1.2fe5b3422585fp+9, 0x1.d15877ba9c4cep+875, 0x1.19aa7b22bded1p-877}, /* i=219 */
   {0x1.3148f0dcb10b7p+9, 0x1.d29dfc774d77dp+879, 0x1.18e5fcaab5a8ep-881}, /* i=220 */
   {0x1.32ac2e775b7fp+9, 0x1.d3e46559f7644p+883, 0x1.18220702ad64ap-885}, /* i=221 */
   {0x1.340f6c11deebp+9, 0x1.d52bb2026285cp+887, 0x1.175e9a644370cp-889}, /* i=222 */
   {0x1.3572a9ac80d98p+9, 0x1.d673e40eb1acfp+891, 0x1.169bb5d80af48p-893}, /* i=223 */
   {0x1.36d5e7471c12ap+9, 0x1.d7bcfb96cda81p+895, 0x1.15d959503dee2p-897}, /* i=224 */
   {0x1.383924e1baecap+9, 0x1.d906f9614e92ep+899, 0x1.1517845784344p-901}, /* i=225 */
   {0x1.399c627c43e55p+9, 0x1.da51ddb0d30f4p+903, 0x1.145636c6693a7p-905}, /* i=226 */
   {0x1.3affa016ee36p+9, 0x1.db9da9f35d4dp+907, 0x1.13956fc74f98bp-909}, /* i=227 */
   {0x1.3c62ddb186072p+9, 0x1.dcea5e0ab2e18p+911, 0x1.12d52f6c0a77bp-913}, /* i=228 */
   {0x1.3dc61b4c1be4cp+9, 0x1.de37fad6718a1p+915, 0x1.1215753306795p-917}, /* i=229 */
   {0x1.3f2958e6bcc4ep+9, 0x1.df868129e5737p+919, 0x1.115640a2fd7e6p-921}, /* i=230 */
   {0x1.408c96815eadp+9, 0x1.e0d5f183283c7p+923, 0x1.10979174122e8p-925}, /* i=231 */
   {0x1.41efd41beba8cp+9, 0x1.e2264c3340c75p+927, 0x1.0fd96777dc177p-929}, /* i=232 */
   {0x1.435311b68cf34p+9, 0x1.e377927969615p+931, 0x1.0f1bc1f9dd497p-933}, /* i=233 */
   {0x1.44b64f5128431p+9, 0x1.e4c9c4974a55p+935, 0x1.0e5ea0d56a4d5p-937}, /* i=234 */
   {0x1.46198cebc2b9ap+9, 0x1.e61ce3453a3e8p+939, 0x1.0da203a34e6ebp-941}, /* i=235 */
   {0x1.477cca8650cafp+9, 0x1.e770eefcb6d87p+943, 0x1.0ce5ea1fbcc0cp-945}, /* i=236 */
   {0x1.48e0082107d4ap+9, 0x1.e8c5e92f159d8p+947, 0x1.0c2a537ee6437p-949}, /* i=237 */
   {0x1.4a4345bb8fffbp+9, 0x1.ea1bd133ad9eep+951, 0x1.0b6f401d7ac7ap-953}, /* i=238 */
   {0x1.4ba6835624f44p+9, 0x1.eb72a894ecc9ap+955, 0x1.0ab4af23147d2p-957}, /* i=239 */
   {0x1.4d09c0f0e2a21p+9, 0x1.ecca706615bd7p+959, 0x1.09fa9ffad3118p-961}, /* i=240 */
   {0x1.4e6cfe8b5de8ap+9, 0x1.ee2327b2a5d4bp+963, 0x1.094113289bb3cp-965}, /* i=241 */
   {0x1.4fd03c2614f72p+9, 0x1.ef7cd10932112p+967, 0x1.0888074bfa2b1p-969}, /* i=242 */
   {0x1.513379c0aebbcp+9, 0x1.f0d76bba92cd8p+971, 0x1.07cf7cc337287p-973}, /* i=243 */
   {0x1.5296b75b552fbp+9, 0x1.f232f9125ca33p+975, 0x1.071772dd80154p-977}, /* i=244 */
   {0x1.53f9f4f5c77afp+9, 0x1.f38f78bd7f123p+979, 0x1.065fe9c626d43p-981}, /* i=245 */
   {0x1.555d329061013p+9, 0x1.f4ececca0707dp+983, 0x1.05a8e067e8e7ep-985}, /* i=246 */
   {0x1.56c0702b02a56p+9, 0x1.f64b5569a5c5ap+987, 0x1.04f256a98a4afp-989}, /* i=247 */
   {0x1.5823adc57929fp+9, 0x1.f7aab27dfd69ap+991, 0x1.043c4c9a46dd9p-993}, /* i=248 */
   {0x1.5986eb603df06p+9, 0x1.f90b06900e4dbp+995, 0x1.0386c0ea762bp-997}, /* i=249 */
   {0x1.5aea28fac317fp+9, 0x1.fa6c501cfa356p+999, 0x1.02d1b461ae4dap-1001}, /* i=250 */
   {0x1.5c4d669567a1p+9, 0x1.fbce9147fe08cp+1003, 0x1.021d25e6f520bp-1005}, /* i=251 */
   {0x1.5db0a42ff2c16p+9, 0x1.fd31c9dd14708p+1007, 0x1.01691594e77fap-1009}, /* i=252 */
   {0x1.5f13e1caaf837p+9, 0x1.fe95fbb417b5ap+1011, 0x1.00b5827cbc444p-1013}, /* i=253 */
   {0x1.60771f653ce2fp+9, 0x1.fffb25f86bd03p+1015, 0x1.00026d09aca19p-1017}, /* i=254 */
   {0x1.61da5cffe0366p+9, 0x1.00b0a5367476ep+1020, 0x1.fe9fa8b27f923p-1022}, /* i=255 */
};

/* For 0 <= i < 256, U[i] = {xi, si, ci} such that xi is near i/magic
   with magic = 0x1.70f77fc88ae3cp6, and si,ci approximate sinh(xi),cosh(xi)
   with accuracy >= 53+11 bits:
   |si - sinh(xi)| < 2^(-11-1) ulp(si), |ci - cosh(xi)| < 2^(-11-1) ulp(ci).
   We have |xi - i/magic| < 3.75e-9.
   Generated with build_table_U(k=11) from accompanying file sinh.sage. */
static const double U[256][3] = {                                                                   
   {0x0p+0, 0x0p+0, 0x1p+0}, /* i=0 */                                                              
   {0x1.633d9a98c86a7p-7, 0x1.633f62a20b32fp-7, 0x1.0003d9ea4bffdp+0}, /* i=1 */                    
   {0x1.633d9a9a4ead4p-6, 0x1.6344bac795757p-6, 0x1.000f67c6da1acp+0}, /* i=2 */                    
   {0x1.0a6e33f133c45p-5, 0x1.0a7a3a54c281dp-5, 0x1.0022a9eea8607p+0}, /* i=3 */                    
   {0x1.633d9a9abb863p-5, 0x1.635a1bd39116cp-5, 0x1.003da0f60f4f2p+0}, /* i=4 */                    
   {0x1.bc0d01432248dp-5, 0x1.bc44ae671ce73p-5, 0x1.00604dacbfc8ap+0}, /* i=5 */                    
   {0x1.0a6e33f18500ep-4, 0x1.0a9e4f73f095p-4, 0x1.008ab11dca39cp+0}, /* i=6 */                     
   {0x1.36d5e748550f8p-4, 0x1.37224d5979f5ep-4, 0x1.00bccc8fb8be6p+0}, /* i=7 */                    
   {0x1.633d9a99ec38bp-4, 0x1.63afa7b941a63p-4, 0x1.00f6a184710c1p+0}, /* i=8 */                    
   {0x1.8fa54de9d9c0dp-4, 0x1.9047b5c025f26p-4, 0x1.013831b968986p+0}, /* i=9 */                    
   {0x1.bc0d0141c4103p-4, 0x1.bcebcef396b7ap-4, 0x1.01817f27af044p+0}, /* i=10 */                   
   {0x1.e874b499f0ee4p-4, 0x1.e99d4b24533cap-4, 0x1.01d28c03d5be2p+0}, /* i=11 */                   
   {0x1.0a6e33f4bf69cp-3, 0x1.0b2ec14481f99p-3, 0x1.022b5abe13b18p+0}, /* i=12 */                   
   {0x1.20a20d9e92c5p-3, 0x1.2196e6efc60a5p-3, 0x1.028bee028a825p+0}, /* i=13 */                    
   {0x1.36d5e7467f673p-3, 0x1.3807c3266eec3p-3, 0x1.02f448b907212p+0}, /* i=14 */                   
   {0x1.4d09c0f1d1d0dp-3, 0x1.4e8202c76aefp-3, 0x1.03646e056c539p+0}, /* i=15 */                    
   {0x1.633d9a9a1f225p-3, 0x1.650652ee1d78cp-3, 0x1.03dc6147662e8p+0}, /* i=16 */                   
   {0x1.797174385396fp-3, 0x1.7b956102b7a21p-3, 0x1.045c261ab1099p+0}, /* i=17 */                   
   {0x1.8fa54debe01acp-3, 0x1.922fdae75988ap-3, 0x1.04e3c0582a49dp+0}, /* i=18 */                   
   {0x1.a5d9279934588p-3, 0x1.a8d66e9a1cfa1p-3, 0x1.05733413aabb8p+0}, /* i=19 */                   
   {0x1.bc0d01411268p-3, 0x1.bf89ca9310f7dp-3, 0x1.060a859e20405p+0}, /* i=20 */                    
   {0x1.d240daebd31bbp-3, 0x1.d64a9db483486p-3, 0x1.06a9b98548d42p+0}, /* i=21 */                   
   {0x1.e874b49405c15p-3, 0x1.ed19973a5f9ddp-3, 0x1.0750d4933e734p+0}, /* i=22 */                   
   {0x1.fea88e3b3238p-3, 0x1.01fbb36a58897p-2, 0x1.07ffdbcf2acb2p+0}, /* i=23 */                    
   {0x1.0a6e33f461f4ap-2, 0x1.0d725e55f20f1p-2, 0x1.08b6d47d71475p+0}, /* i=24 */                   
   {0x1.158820c3591a4p-2, 0x1.18f124a0c1da8p-2, 0x1.0975c41eabf04p+0}, /* i=25 */                   
   {0x1.20a20d9bd0299p-2, 0x1.24785ee755796p-2, 0x1.0a3cb072b3ff4p+0}, /* i=26 */                   
   {0x1.2bbbfa746d378p-2, 0x1.300865ebba992p-2, 0x1.0b0b9f751f0cap+0}, /* i=27 */                   
   {0x1.36d5e7457f38ep-2, 0x1.3ba192b565882p-2, 0x1.0be2975f3bbc4p+0}, /* i=28 */                   
   {0x1.41efd4120cd8ap-2, 0x1.47443e9d65df8p-2, 0x1.0cc19ea900cbdp+0}, /* i=29 */                   
   {0x1.4d09c0f323ef2p-2, 0x1.52f0c35d2e861p-2, 0x1.0da8bc0a55676p+0}, /* i=30 */                   
   {0x1.5823adc60e321p-2, 0x1.5ea77abb3ae2bp-2, 0x1.0e97f67494f5ap+0}, /* i=31 */                   
   {0x1.633d9a9a861e2p-2, 0x1.6a68bf01925f5p-2, 0x1.0f8f551ba9fe7p+0}, /* i=32 */                   
   {0x1.6e57876f71e8ep-2, 0x1.7634eab9d1244p-2, 0x1.108edf70d95d1p+0}, /* i=33 */                   
   {0x1.7971744450862p-2, 0x1.820c58c22f55fp-2, 0x1.11969d24622c8p+0}, /* i=34 */                   
   {0x1.848b6119e6182p-2, 0x1.8def6450f8d2fp-2, 0x1.12a69625cbf0ap+0}, /* i=35 */                   
   {0x1.8fa54deb08b2bp-2, 0x1.99de68ef911e6p-2, 0x1.13bed2a36ffc9p+0}, /* i=36 */                   
   {0x1.9abf3ac319177p-2, 0x1.a5d9c29561345p-2, 0x1.14df5b0cf2c6cp+0}, /* i=37 */                   
   {0x1.a5d927953725ap-2, 0x1.b1e1cd7eab67p-2, 0x1.1608380f709efp+0}, /* i=38 */                    
   {0x1.b0f3146e0b7e8p-2, 0x1.bdf6e6651348cp-2, 0x1.1739729ad10d5p+0}, /* i=39 */                   
   {0x1.bc0d0141b532dp-2, 0x1.ca196a4c251e6p-2, 0x1.187313dcd26cap+0}, /* i=40 */                   
   {0x1.c726ee1426721p-2, 0x1.d649b6aff6ce8p-2, 0x1.19b5254595f8dp+0}, /* i=41 */                   
   {0x1.d240dae5e95ep-2, 0x1.e288297317a9fp-2, 0x1.1affb085e0a8cp+0}, /* i=42 */                    
   {0x1.dd5ac7c1b4d51p-2, 0x1.eed520f063cp-2, 0x1.1c52bf90fce33p+0}, /* i=43 */                     
   {0x1.e874b49449e79p-2, 0x1.fb30fbd155be6p-2, 0x1.1dae5c9836286p+0}, /* i=44 */                   
   {0x1.f38ea1697dde3p-2, 0x1.03ce0caa4344fp-1, 0x1.1f129212720fdp+0}, /* i=45 */                   
   {0x1.fea88e3ee44b9p-2, 0x1.0a0b6c8f251bep-1, 0x1.207f6ab715a7cp+0}, /* i=46 */                   
   {0x1.04e13d8931209p-1, 0x1.1050cda5d776ep-1, 0x1.21f4f18010b63p+0}, /* i=47 */                   
   {0x1.0a6e33f32e8e2p-1, 0x1.169e603d12f39p-1, 0x1.237331aac51c7p+0}, /* i=48 */                   
   {0x1.0ffb2a5d59b2ap-1, 0x1.1cf454e14119fp-1, 0x1.24fa36b774d0ep+0}, /* i=49 */                   
   {0x1.158820c5db152p-1, 0x1.2352dc5d51c23p-1, 0x1.268a0c696b231p+0}, /* i=50 */                   
   {0x1.1b15173180aa6p-1, 0x1.29ba27c5e2fa8p-1, 0x1.2822bec9bc085p+0}, /* i=51 */                   
   {0x1.20a20d9d14bdep-1, 0x1.302a6869e30e4p-1, 0x1.29c45a237b4afp+0}, /* i=52 */                   
   {0x1.262f04046a142p-1, 0x1.36a3cfdc09902p-1, 0x1.2b6eeb05fe1fdp+0}, /* i=53 */                   
   {0x1.2bbbfa6c92e0ep-1, 0x1.3d2690004ffcp-1, 0x1.2d227e485dac6p+0}, /* i=54 */                    
   {0x1.3148f0df7976p-1, 0x1.43b2db088aee1p-1, 0x1.2edf2108c0d4p+0}, /* i=55 */                     
   {0x1.36d5e74898cf9p-1, 0x1.4a48e34c1d8f7p-1, 0x1.30a4e0a1887adp+0}, /* i=56 */                   
   {0x1.3c62ddb2d88fdp-1, 0x1.50e8db927b60bp-1, 0x1.3273cabd15bebp+0}, /* i=57 */                   
   {0x1.41efd41bcc36p-1, 0x1.5792f6dfe4274p-1, 0x1.344bed4832e5ap+0}, /* i=58 */                    
   {0x1.477cca862041ap-1, 0x1.5e47688cc6ff4p-1, 0x1.362d56785d9b9p+0}, /* i=59 */                   
   {0x1.4d09c0f2c4e23p-1, 0x1.6506643f291b9p-1, 0x1.381814c9fd313p+0}, /* i=60 */                   
   {0x1.5296b759425f3p-1, 0x1.6bd01de2babefp-1, 0x1.3a0c36fe0a17p+0}, /* i=61 */                    
   {0x1.5823adc4f6df3p-1, 0x1.72a4c9ce2392ap-1, 0x1.3c09cc24c4e9dp+0}, /* i=62 */           
   {0x1.5db0a431bdf8ep-1, 0x1.79849c9971258p-1, 0x1.3e10e391cbfebp+0}, /* i=63 */                   
   {0x1.633d9aa2d495fp-1, 0x1.806fcb3bb01f2p-1, 0x1.40218ce4a48d1p+0}, /* i=64 */                   
   {0x1.68ca910445a8p-1, 0x1.87668ae653b8p-1, 0x1.423bd7fdc2f6dp+0}, /* i=65 */                     
   {0x1.6e57876c233abp-1, 0x1.8e6911587dd09p-1, 0x1.445fd5179acddp+0}, /* i=66 */                   
   {0x1.73e47de15534cp-1, 0x1.9577949969fc2p-1, 0x1.468d94b1e690ap+0}, /* i=67 */                   
   {0x1.7971743d66c5dp-1, 0x1.9c924ad2a5a87p-1, 0x1.48c5278577eep+0}, /* i=68 */                    
   {0x1.7efe6aaba5d7fp-1, 0x1.a3b96af416c7ap-1, 0x1.4b069eb73f6d5p+0}, /* i=69 */                   
   {0x1.848b610f2c6c3p-1, 0x1.aaed2bf0c5425p-1, 0x1.4d520b9b7a8ddp+0}, /* i=70 */                   
   {0x1.8a1857842123fp-1, 0x1.b22dc566c1b64p-1, 0x1.4fa77fea7e103p+0}, /* i=71 */                   
   {0x1.8fa54deda7092p-1, 0x1.b97b6f0ccdcf1p-1, 0x1.52070d9216085p+0}, /* i=72 */                   
   {0x1.9532445346abcp-1, 0x1.c0d6612d9e95p-1, 0x1.5470c6dcf94e4p+0}, /* i=73 */                    
   {0x1.9abf3ac13f51dp-1, 0x1.c83ed480d0c65p-1, 0x1.56e4be6685c46p+0}, /* i=74 */                   
   {0x1.a04c312dfad74p-1, 0x1.cfb5020909debp-1, 0x1.5963070f76e61p+0}, /* i=75 */                   
   {0x1.a5d92799beeffp-1, 0x1.d739233fb2613p-1, 0x1.5bebb40c4164p+0}, /* i=76 */                    
   {0x1.ab661dff533f1p-1, 0x1.decb720227066p-1, 0x1.5e7ed8ded58afp+0}, /* i=77 */                   
   {0x1.b0f3146be7bcdp-1, 0x1.e66c28b42647dp-1, 0x1.611c89625ba27p+0}, /* i=78 */                   
   {0x1.b6800ad49347ep-1, 0x1.ee1b820771d39p-1, 0x1.63c4d9b7fa2d7p+0}, /* i=79 */                   
   {0x1.bc0d0141b5be9p-1, 0x1.f5d9b93912d7dp-1, 0x1.6677de5bccbeep+0}, /* i=80 */                   
   {0x1.c199f7ab79a4ap-1, 0x1.fda709e211184p-1, 0x1.6935ac1491c72p+0}, /* i=81 */                   
   {0x1.c726ee14596edp-1, 0x1.02c1d80f0acc7p+0, 0x1.6bfe580113e3ap+0}, /* i=82 */                   
   {0x1.ccb3e482e518dp-1, 0x1.06b7f4426da4p+0, 0x1.6ed1f795ff547p+0}, /* i=83 */                    
   {0x1.d240daeac2fbap-1, 0x1.0ab5f80583b64p+0, 0x1.71b0a092b58b9p+0}, /* i=84 */                   
   {0x1.d7cdd154dfe0cp-1, 0x1.0ebc021eb6936p+0, 0x1.749a691677c4dp+0}, /* i=85 */                   
   {0x1.dd5ac7be36a6dp-1, 0x1.12ca3189b0129p+0, 0x1.778f679012bfap+0}, /* i=86 */                   
   {0x1.e2e7be2bbc903p-1, 0x1.16e0a586ac6bcp+0, 0x1.7a8fb2c8ecfc9p+0}, /* i=87 */                   
   {0x1.e874b491a8b4ap-1, 0x1.1aff7d893a0e7p+0, 0x1.7d9b61d85de19p+0}, /* i=88 */                   
   {0x1.ee01aafc407e3p-1, 0x1.1f26d957537dfp+0, 0x1.80b28c3a94754p+0}, /* i=89 */                   
   {0x1.f38ea16b4ebcp-1, 0x1.2356d8ef60fc2p+0, 0x1.83d549bd7be6p+0}, /* i=90 */                     
   {0x1.f91b97d0d8451p-1, 0x1.278f9c87d1714p+0, 0x1.8703b28031fb2p+0}, /* i=91 */                   
   {0x1.fea88e3cbaefp-1, 0x1.2bd144b13a4adp+0, 0x1.8a3ddf0cb9b06p+0}, /* i=92 */                    
   {0x1.021ac252eac31p+0, 0x1.301bf22dae5a3p+0, 0x1.8d83e839627abp+0}, /* i=93 */                   
   {0x1.04e13d85e05bfp+0, 0x1.346fc60b551cfp+0, 0x1.90d5e73cccbcp+0}, /* i=94 */                    
   {0x1.07a7b8c29a97ap+0, 0x1.38cce1b10eb63p+0, 0x1.9433f5b7d585p+0}, /* i=95 */                    
   {0x1.0a6e33f8842bcp+0, 0x1.3d3366a152823p+0, 0x1.979e2d8666738p+0}, /* i=96 */                   
   {0x1.0d34af2904b1ep+0, 0x1.41a376c31539ap+0, 0x1.9b14a8f78b123p+0}, /* i=97 */                   
   {0x1.0ffb2a5fffa63p+0, 0x1.461d3457dac3bp+0, 0x1.9e9782c6389eap+0}, /* i=98 */                   
   {0x1.12c1a593e7b51p+0, 0x1.4aa0c1c9536e9p+0, 0x1.a226d5f1cd6eap+0}, /* i=99 */                   
   {0x1.158820c8f4c5p+0, 0x1.4f2e41e2eb096p+0, 0x1.a5c2bdeb17d8cp+0}, /* i=100 */                   
   {0x1.184e9c0184efdp+0, 0x1.53c5d7b9d6532p+0, 0x1.a96b5681a631dp+0}, /* i=101 */                  
   {0x1.1b151733df05ep+0, 0x1.5867a69cd6832p+0, 0x1.ad20bbd6acf2fp+0}, /* i=102 */                  
   {0x1.1ddb926b842abp+0, 0x1.5d13d24cbef1dp+0, 0x1.b0e30a8a351b6p+0}, /* i=103 */                  
   {0x1.20a20d9ca2012p+0, 0x1.61ca7eb316c2ep+0, 0x1.b4b25f8055ba4p+0}, /* i=104 */                  
   {0x1.236888d4c1915p+0, 0x1.668bd035396bp+0, 0x1.b88ed8241c453p+0}, /* i=105 */                   
   {0x1.262f04082ee1ep+0, 0x1.6b57eb5f8031ep+0, 0x1.bc789222fa6c8p+0}, /* i=106 */                  
   {0x1.28f57f3c81e23p+0, 0x1.702ef52ef1207p+0, 0x1.c06faba87ac76p+0}, /* i=107 */                  
   {0x1.2bbbfa6f0cf1dp+0, 0x1.751112e68081p+0, 0x1.c474433b879c1p+0}, /* i=108 */                   
   {0x1.2e8275a259035p+0, 0x1.79fe6a2796cbcp+0, 0x1.c88677d275fb5p+0}, /* i=109 */                  
   {0x1.3148f0dc8f6efp+0, 0x1.7ef720f0c1916p+0, 0x1.cca668d214e81p+0}, /* i=110 */                  
   {0x1.340f6c0d332bep+0, 0x1.83fb5d6f6edfp+0, 0x1.d0d435e7282c1p+0}, /* i=111 */                   
   {0x1.36d5e7468fa39p+0, 0x1.890b466873bddp+0, 0x1.d50fff5d503fcp+0}, /* i=112 */                  
   {0x1.399c627a3052ap+0, 0x1.8e2702bf37f28p+0, 0x1.d959e5ba6ac08p+0}, /* i=113 */                  
   {0x1.3c62ddaeed33fp+0, 0x1.934eb9d8ffcccp+0, 0x1.ddb20a1170228p+0}, /* i=114 */                  
   {0x1.3f2958e83e783p+0, 0x1.9882937171072p+0, 0x1.e2188dde36353p+0}, /* i=115 */                  
   {0x1.41efd41b411f1p+0, 0x1.9dc2b7869dee3p+0, 0x1.e68d92f45c9c7p+0}, /* i=116 */                  
   {0x1.44b64f4ebf6dp+0, 0x1.a30f4e966a5d3p+0, 0x1.eb113bb346f59p+0}, /* i=117 */                   
   {0x1.477cca8b40f6fp+0, 0x1.a868818275ebp+0, 0x1.efa3aaee8cb9bp+0}, /* i=118 */                   
   {0x1.4a4345be511bdp+0, 0x1.adce795960664p+0, 0x1.f44503befd4cp+0}, /* i=119 */                   
   {0x1.4d09c0ee1d3edp+0, 0x1.b3415fbb647bp+0, 0x1.f8f569d897643p+0}, /* i=120 */                   
   {0x1.4fd03c2534368p+0, 0x1.b8c15eb58334ap+0, 0x1.fdb5016b59824p+0}, /* i=121 */                  
   {0x1.5296b75c5057dp+0, 0x1.be4ea09690717p+0, 0x1.0141f77f0e053p+1}, /* i=122 */                  
   {0x1.555d3292dd7dep+0, 0x1.c3e950208ea57p+0, 0x1.03b12bcc83ad2p+1}, /* i=123 */                  
   {0x1.5823adc56efe6p+0, 0x1.c99198770f432p+0, 0x1.0628305b2edbbp+1}, /* i=124 */                  
   {0x1.5aea28fef35b9p+0, 0x1.cf47a543c8eb2p+0, 0x1.08a7183135f32p+1}, /* i=125 */                  
   {0x1.5db0a42e51e6ap+0, 0x1.d50ba260d7478p+0, 0x1.0b2df67850bf8p+1}, /* i=126 */                  
   {0x1.60771f64205d7p+0, 0x1.daddbc59526d5p+0, 0x1.0dbcdeb603887p+1}, /* i=127 */             
   {0x1.633d9a99376b7p+0, 0x1.e0be1ff3d8625p+0, 0x1.1053e49834f43p+1}, /* i=128 */                  
   {0x1.660415d058621p+0, 0x1.e6acfa79c93d4p+0, 0x1.12f31c1466062p+1}, /* i=129 */                  
   {0x1.68ca910a2cfcdp+0, 0x1.ecaa799f9f357p+0, 0x1.159a995d51ebbp+1}, /* i=130 */                  
   {0x1.6b910c384e772p+0, 0x1.f2b6cb69b02a4p+0, 0x1.184a70d6c7f64p+1}, /* i=131 */                  
   {0x1.6e57877224038p+0, 0x1.f8d21ea0dcb9fp+0, 0x1.1b02b74983fa8p+1}, /* i=132 */                  
   {0x1.711e02a61d2a4p+0, 0x1.fefca22821f5p+0, 0x1.1dc38197727f2p+1}, /* i=133 */                   
   {0x1.73e47ddba7616p+0, 0x1.029b42c7276fcp+1, 0x1.208ce4fc7cedcp+1}, /* i=134 */                  
   {0x1.76aaf90ead5b7p+0, 0x1.05bffc5f87d82p+1, 0x1.235ef6eb5dd67p+1}, /* i=135 */                  
   {0x1.79717445f3768p+0, 0x1.08ec961b4105p+1, 0x1.2639cd24b5d35p+1}, /* i=136 */                   
   {0x1.7c37ef793a5fep+0, 0x1.0c212863eb8aap+1, 0x1.291d7d9d486a7p+1}, /* i=137 */                  
   {0x1.7efe6aaecdfaep+0, 0x1.0f5dcbf131853p+1, 0x1.2c0a1e9d09a77p+1}, /* i=138 */                  
   {0x1.81c4e5e4e30ap+0, 0x1.12a299af96f72p+1, 0x1.2effc6a86c7ccp+1}, /* i=139 */                   
   {0x1.848b61182609p+0, 0x1.15efaac8999aep+1, 0x1.31fe8c87ab73ep+1}, /* i=140 */                   
   {0x1.8751dc4d51c43p+0, 0x1.194518aeedf79p+1, 0x1.35068751e707cp+1}, /* i=141 */                  
   {0x1.8a18578252a1bp+0, 0x1.1ca2fd0d9f326p+1, 0x1.3817ce5dd78a6p+1}, /* i=142 */                  
   {0x1.8cded2b8bfc87p+0, 0x1.200971d552106p+1, 0x1.3b32794dd7846p+1}, /* i=143 */                  
   {0x1.8fa54deb7a686p+0, 0x1.237891305de27p+1, 0x1.3e56a00507227p+1}, /* i=144 */                  
   {0x1.926bc923efed9p+0, 0x1.26f0759ddf2fbp+1, 0x1.41845ac00e88dp+1}, /* i=145 */                  
   {0x1.9532445a670a2p+0, 0x1.2a7139cb35ba6p+1, 0x1.44bbc1f1e0e15p+1}, /* i=146 */                  
   {0x1.97f8bf8e0fca8p+0, 0x1.2dfaf8b293106p+1, 0x1.47fcee5faf504p+1}, /* i=147 */                  
   {0x1.9abf3ac404a2bp+0, 0x1.318dcd9aef4bep+1, 0x1.4b47f920f0419p+1}, /* i=148 */                  
   {0x1.9d85b5f6a068cp+0, 0x1.3529d403752bdp+1, 0x1.4e9cfb8c6b03dp+1}, /* i=149 */                  
   {0x1.a04c312d61163p+0, 0x1.38cf27c32beb3p+1, 0x1.51fc0f5568c03p+1}, /* i=150 */                  
   {0x1.a312ac62cba11p+0, 0x1.3c7de4e7d9654p+1, 0x1.55654e6d2266fp+1}, /* i=151 */                  
   {0x1.a5d9279637a5cp+0, 0x1.403627cdf6428p+1, 0x1.58d8d318dc8a9p+1}, /* i=152 */                  
   {0x1.a89fa2cb93b1cp+0, 0x1.43f80d2189342p+1, 0x1.5c56b7f2bec04p+1}, /* i=153 */                  
   {0x1.ab661e00b7c32p+0, 0x1.47c3b1d0a73dap+1, 0x1.5fdf17dd4b702p+1}, /* i=154 */                  
   {0x1.ae2c993667599p+0, 0x1.4b99331869535p+1, 0x1.63720e0f66da1p+1}, /* i=155 */                  
   {0x1.b0f3146b13af2p+0, 0x1.4f78ae7ea778fp+1, 0x1.670fb60e82129p+1}, /* i=156 */                  
   {0x1.b3b98fa1dcdfcp+0, 0x1.536241dc9ec12p+1, 0x1.6ab82bb88b28cp+1}, /* i=157 */                  
   {0x1.b6800ad6924d5p+0, 0x1.57560b4f0096ap+1, 0x1.6e6b8b3508169p+1}, /* i=158 */                  
   {0x1.b946860a8a5bcp+0, 0x1.5b542948e52abp+1, 0x1.7229f106cba3dp+1}, /* i=159 */                  
   {0x1.bc0d0130f6eb5p+0, 0x1.5f5cba78612bep+1, 0x1.75f379f23f6d2p+1}, /* i=160 */                  
   {0x1.bed37c7ae67ddp+0, 0x1.636fde368d748p+1, 0x1.79c84366a934p+1}, /* i=161 */                   
   {0x1.c199f7acb0c3bp+0, 0x1.678db38ed1beap+1, 0x1.7da86a946c4c6p+1}, /* i=162 */                  
   {0x1.c46072e077394p+0, 0x1.6bb65a5c8d32ap+1, 0x1.81940d7983e49p+1}, /* i=163 */                  
   {0x1.c726ee164eee7p+0, 0x1.6fe9f2a86779ep+1, 0x1.858b4a489f88p+1}, /* i=164 */                   
   {0x1.c9ed694a8b281p+0, 0x1.74289cc98dcf1p+1, 0x1.898e3f885606bp+1}, /* i=165 */                  
   {0x1.ccb3e480392edp+0, 0x1.78727976a883bp+1, 0x1.8d9d0c232ac1ap+1}, /* i=166 */                  
   {0x1.cf7a5fb5a188cp+0, 0x1.7cc7a9b555571p+1, 0x1.91b7cf57f2a03p+1}, /* i=167 */                  
   {0x1.d240daedeccd4p+0, 0x1.81284eea0e7c3p+1, 0x1.95dea8c8e0317p+1}, /* i=168 */                  
   {0x1.d507561c8dd03p+0, 0x1.85948ab897d23p+1, 0x1.9a11b85d8e9c3p+1}, /* i=169 */                  
   {0x1.d7cdd15deffe4p+0, 0x1.8a0c7f5f72673p+1, 0x1.9e511e99d97e3p+1}, /* i=170 */                  
   {0x1.da944c8c4ac62p+0, 0x1.8e904f0d0c5edp+1, 0x1.a29cfbfb916eap+1}, /* i=171 */                  
   {0x1.dd5ac7bf3a27cp+0, 0x1.93201cae095e1p+1, 0x1.a6f571be8c6c2p+1}, /* i=172 */                  
   {0x1.e02142f36e6bep+0, 0x1.97bc0b6013ba8p+1, 0x1.ab5aa15671475p+1}, /* i=173 */                  
   {0x1.e2e7be2c42c0dp+0, 0x1.9c643ea947ebcp+1, 0x1.afccaca36d597p+1}, /* i=174 */                  
   {0x1.e5ae3960d58d2p+0, 0x1.a118da59a210ep+1, 0x1.b44bb5d5075abp+1}, /* i=175 */                  
   {0x1.e874b493acdccp+0, 0x1.a5da02b3b440fp+1, 0x1.b8d7df90fb665p+1}, /* i=176 */                  
   {0x1.eb3b2fc8d77bp+0, 0x1.aaa7dc5d92941p+1, 0x1.bd714ce4dfea5p+1}, /* i=177 */                   
   {0x1.ee01aafe2f062p+0, 0x1.af828c54759d8p+1, 0x1.c218213a4d396p+1}, /* i=178 */                  
   {0x1.f0c826336ecf7p+0, 0x1.b46a37fba694fp+1, 0x1.c6cc806522883p+1}, /* i=179 */                  
   {0x1.f38ea16b17512p+0, 0x1.b95f051f63484p+1, 0x1.cb8e8ea64f9b8p+1}, /* i=180 */                  
   {0x1.f6551c9dc6b4bp+0, 0x1.be6119df4f662p+1, 0x1.d05e7097208cfp+1}, /* i=181 */                  
   {0x1.f91b97d1b17d6p+0, 0x1.c3709cd9b0fbep+1, 0x1.d53c4b52c3329p+1}, /* i=182 */                  
   {0x1.fbe2130c93ca6p+0, 0x1.c88db513da62ap+1, 0x1.da28445faa121p+1}, /* i=183 */                  
   {0x1.fea88e3eddd3p+0, 0x1.cdb889d453a2cp+1, 0x1.df22818b14596p+1}, /* i=184 */                   
   {0x1.00b784b8af01dp+1, 0x1.d2f142f85260bp+1, 0x1.e42b293b53acfp+1}, /* i=185 */                  
   {0x1.021ac254a5d9cp+1, 0x1.d83808c1e1135p+1, 0x1.e942623fdb79p+1}, /* i=186 */                   
   {0x1.037dffeb2a036p+1, 0x1.dd8d03b612f18p+1, 0x1.ee6853b094d3ep+1}, /* i=187 */                  
   {0x1.04e13d8c06674p+1, 0x1.e2f05d2398d9dp+1, 0x1.f39d256fcea56p+1}, /* i=188 */                  
   {0x1.06447b2369b6bp+1, 0x1.e8623e3def6fp+1, 0x1.f8e0ff4d394e9p+1}, /* i=189 */                   
   {0x1.07a7b8c043a8dp+1, 0x1.ede2d12fc1833p+1, 0x1.fe340a0f060b2p+1}, /* i=190 */
   {0x1.090af65703af7p+1, 0x1.f372402ca6117p+1, 0x1.01cb3745c8d1bp+2}, /* i=191 */
   {0x1.0a6e33f139dadp+1, 0x1.f910b62eb5aa9p+1, 0x1.04842b3067cdap+2}, /* i=192 */
   {0x1.0bd1718d207fep+1, 0x1.febe5e76c83b9p+1, 0x1.0744f5bc51eb1p+2}, /* i=193 */
   {0x1.0d34af2f4f791p+1, 0x1.023db26a75ebbp+2, 0x1.0a0dac278921ep+2}, /* i=194 */
   {0x1.0e97ecc7026d8p+1, 0x1.0523fa9b4fe4p+2, 0x1.0cde63c1e9cb3p+2}, /* i=195 */
   {0x1.0ffb2a5fb7c3ep+1, 0x1.08121e3b0cdbp+2, 0x1.0fb73251bf95fp+2}, /* i=196 */
   {0x1.115e67f8b4616p+1, 0x1.0b0833da17a6ep+2, 0x1.12982dc33148ap+2}, /* i=197 */
   {0x1.12c1a59294a3dp+1, 0x1.0e065248e5e69p+2, 0x1.15816c4421702p+2}, /* i=198 */
   {0x1.1424e32ec969fp+1, 0x1.110c9097afc46p+2, 0x1.18730443e9e92p+2}, /* i=199 */
   {0x1.158820c75ceep+1, 0x1.141b060500d77p+2, 0x1.1b6d0c6261145p+2}, /* i=200 */
   {0x1.16eb5e668213bp+1, 0x1.1731ca324681ap+2, 0x1.1e6f9ba308565p+2}, /* i=201 */
   {0x1.184e9bfec0769p+1, 0x1.1a50f4cc5a5c3p+2, 0x1.217ac917d91bap+2}, /* i=202 */
   {0x1.19b1d9951ff31p+1, 0x1.1d789de9ca403p+2, 0x1.248eac3d20184p+2}, /* i=203 */
   {0x1.1b15173098d25p+1, 0x1.20a8dde786607p+2, 0x1.27ab5cd71d4c9p+2}, /* i=204 */
   {0x1.1c7854ccc0d4bp+1, 0x1.23e1cd4ac20aap+2, 0x1.2ad0f2d49a5b5p+2}, /* i=205 */
   {0x1.1ddb926ab6f24p+1, 0x1.272384e8379abp+2, 0x1.2dff86753931dp+2}, /* i=206 */
   {0x1.1f3ed002b7572p+1, 0x1.2a6e1dc333157p+2, 0x1.3137302942092p+2}, /* i=207 */
   {0x1.20a20d99d4db2p+1, 0x1.2dc1b1418dbfap+2, 0x1.347808c470648p+2}, /* i=208 */
   {0x1.22054b3e2be91p+1, 0x1.311e5924d004ap+2, 0x1.37c229774c797p+2}, /* i=209 */
   {0x1.236888d5bf2fdp+1, 0x1.34842f1331e11p+2, 0x1.3b15ab5ab3286p+2}, /* i=210 */
   {0x1.24cbc66de92cfp+1, 0x1.37f34d58cd5acp+2, 0x1.3e72a82cdf2d9p+2}, /* i=211 */
   {0x1.262f04030bef6p+1, 0x1.3b6bce5ff5dd6p+2, 0x1.41d939ccac1b4p+2}, /* i=212 */
   {0x1.279241a7c8703p+1, 0x1.3eedcd136f7abp+2, 0x1.45497a99cac08p+2}, /* i=213 */
   {0x1.28f57f3789badp+1, 0x1.4279641b4dc3cp+2, 0x1.48c384b5797dep+2}, /* i=214 */
   {0x1.2a58bcd78d0ecp+1, 0x1.460eaf2532881p+2, 0x1.4c47734453954p+2}, /* i=215 */
   {0x1.2bbbfa6cde099p+1, 0x1.49adc9852784p+2, 0x1.4fd56115fbb55p+2}, /* i=216 */
   {0x1.2d1f380c3b6e7p+1, 0x1.4d56cf5644af5p+2, 0x1.536d69c040a1dp+2}, /* i=217 */
   {0x1.2e8275a6f7136p+1, 0x1.5109dca36ce98p+2, 0x1.570fa8cbdedd7p+2}, /* i=218 */
   {0x1.2fe5b344c1cd7p+1, 0x1.54c70dffb26c3p+2, 0x1.5abc3a4a1933dp+2}, /* i=219 */
   {0x1.3148f0dc27dcbp+1, 0x1.588e801e3616ep+2, 0x1.5e733a6e7473ep+2}, /* i=220 */
   {0x1.32ac2e78c52b8p+1, 0x1.5c60503a3fe98p+2, 0x1.6234c5f5039bcp+2}, /* i=221 */
   {0x1.340f6c12dd013p+1, 0x1.603c9baa332f6p+2, 0x1.6600f9b737a3p+2}, /* i=222 */
   {0x1.3572a9add8d49p+1, 0x1.64238033c9198p+2, 0x1.69d7f2feae2dap+2}, /* i=223 */
   {0x1.36d5e7479cc8p+1, 0x1.68151bdefecd6p+2, 0x1.6db9cf58defcep+2}, /* i=224 */
   {0x1.383924e2f5214p+1, 0x1.6c118d14443efp+2, 0x1.71a6acb4cf446p+2}, /* i=225 */
   {0x1.399c627a29929p+1, 0x1.7018f276f7a3p+2, 0x1.759ea93e20b6dp+2}, /* i=226 */
   {0x1.3affa017d0ab5p+1, 0x1.742b6b2e3cbbfp+2, 0x1.79a1e3a4d0dep+2}, /* i=227 */
   {0x1.3c62ddadf36aep+1, 0x1.7849166ea3fb2p+2, 0x1.7db07aa89ff5bp+2}, /* i=228 */
   {0x1.3dc61b4bb6b75p+1, 0x1.7c721418f30fcp+2, 0x1.81ca8db589caep+2}, /* i=229 */
   {0x1.3f2958e576856p+1, 0x1.80a684157d51cp+2, 0x1.85f03c417abdep+2}, /* i=230 */
   {0x1.408c9680f57b1p+1, 0x1.84e686d8d3852p+2, 0x1.8a21a64f1e556p+2}, /* i=231 */
   {0x1.41efd41bdc3a9p+1, 0x1.89323d1835306p+2, 0x1.8e5eec235f0d2p+2}, /* i=232 */
   {0x1.435311b86ca8dp+1, 0x1.8d89c7f12a966p+2, 0x1.92a82e6c7651ap+2}, /* i=233 */
   {0x1.44b64f4db3c4ep+1, 0x1.91ed48b9076fdp+2, 0x1.96fd8e120e81ep+2}, /* i=234 */
   {0x1.46198cec980d5p+1, 0x1.965ce1731b736p+2, 0x1.9b5f2ca9f6c8fp+2}, /* i=235 */
   {0x1.477cca89d86d8p+1, 0x1.9ad8b426b914cp+2, 0x1.9fcd2bd04bd4ep+2}, /* i=236 */
   {0x1.48e0082203f63p+1, 0x1.9f60e3517931ep+2, 0x1.a447ad984871ap+2}, /* i=237 */
   {0x1.4a4345ba8fefep+1, 0x1.a3f591ed36b41p+2, 0x1.a8ced49237afep+2}, /* i=238 */
   {0x1.4ba6835a9b848p+1, 0x1.a896e35a69928p+2, 0x1.ad62c3b616d64p+2}, /* i=239 */
   {0x1.4d09c0f47afd4p+1, 0x1.ad44fb156485bp+2, 0x1.b2039e19a5d66p+2}, /* i=240 */
   {0x1.4e6cfe8b2083bp+1, 0x1.b1fffd33d8628p+2, 0x1.b6b1876c7dad2p+2}, /* i=241 */
   {0x1.4fd03c2059898p+1, 0x1.b6c80e2b674fap+2, 0x1.bb6ca3bf4e011p+2}, /* i=242 */
   {0x1.513379bf03e12p+1, 0x1.bb9d52f647ce1p+2, 0x1.c03517a81ec3ep+2}, /* i=243 */
   {0x1.5296b75cf200fp+1, 0x1.c07ff0ab87cf5p+2, 0x1.c50b07dba6557p+2}, /* i=244 */
   {0x1.53f9f4f3aef89p+1, 0x1.c5700cd4cbf33p+2, 0x1.c9ee99821cc09p+2}, /* i=245 */
   {0x1.555d329059debp+1, 0x1.ca6dcda9762ebp+2, 0x1.cedff271d45bcp+2}, /* i=246 */
   {0x1.56c0702b18b5fp+1, 0x1.cf79597f5c8bp+2, 0x1.d3df38a15bc2fp+2}, /* i=247 */
   {0x1.5823adc5ef038p+1, 0x1.d492d739cc50fp+2, 0x1.d8ec929567baep+2}, /* i=248 */
   {0x1.5986eb5c46241p+1, 0x1.d9ba6e0f24a65p+2, 0x1.de082726fd69p+2}, /* i=249 */
   {0x1.5aea28f9ba286p+1, 0x1.def045dec45eap+2, 0x1.e3321dd8981cfp+2}, /* i=250 */
   {0x1.5c4d669417e8cp+1, 0x1.e43486a47d3edp+2, 0x1.e86a9e4ae1664p+2}, /* i=251 */
   {0x1.5db0a42fc7bdbp+1, 0x1.e9875901cce2cp+2, 0x1.edb1d0c4adc98p+2}, /* i=252 */
   {0x1.5f13e1ca674d5p+1, 0x1.eee8e5eea731cp+2, 0x1.f307dde47b841p+2}, /* i=253 */
   {0x1.60771f67fbd14p+1, 0x1.f45956ed4a1cdp+2, 0x1.f86ceed3d4f1dp+2}, /* i=254 */
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
  0x1p0,                 /* degree 1 */
  0x1.5555555555555p-3,  /* degree 3 */
  0x1.11111111869d4p-7,  /* degree 5 */
  0x1.a01061b363a81p-13, /* degree 7 */
};

/* the following is a degree-6 even polynomial approximating cosh(x)
   for |x| < 0.00543 generated by Sollya (see Pcosh.sollya), with
   maximal relative/absolute error 2^-81.152 */
static const double C[] = {
  0x1p0,                 /* degree 0 */
  0x1p-1,                /* degree 2 */
  0x1.5555555554e2ep-5,  /* degree 4 */
  0x1.6c16d52a52a35p-10, /* degree 6 */
};

/* put in h+l a double-double approximation of sinh(w), for |w| < 0.00543,
   with maximal relative error 2^-67.99 (see analyze_eval_S_all(rel=true)
   from accompanying file sinh.sage) */
static void
eval_S (double *h, double *l, double w)
{
  double z = w * w;
  *h = __builtin_fma (S[3], z, S[2]);
  *h = __builtin_fma (*h, z, S[1]);
  *h = *h * z; /* h approximates w^2*(S[1]+w^2*S[2]+w^4*S[2]) */
  /* we use the fact that S[0]=1 here, thus we add w + w*h */
  fast_two_sum (h, l, w, *h * w);
}

/* put in h+l a double-double approximation of cosh(w), for |w| < 0.00543,
   with maximal absolute error 2^-68.04 (see analyze_eval_C() from
   accompanying file sinh.sage). Since |cosh(w)| > 1, this is also a bound
   on the relative error. */
static void
eval_C (double *h, double *l, double w)
{
  double z = w * w;
  *h = __builtin_fma (C[3], z, C[2]);
  *h = __builtin_fma (*h, z, C[1]);
  *h = *h * z; /* h approximates w^2*(C[1]+w^2*C[2]+w^4*C[2]) */
  /* we use the fact that C[0]=1 here, thus we add 1 + h */
  fast_two_sum (h, l, 1.0, *h);
}

/* Put in h+l a double-double approximation of sinh(x),
   for 0 <= x <= 0x1.633ce8fb9f87dp+9.
   Return the absolute error bound:
   |h + l - sin(x)| < err. */
static double
cr_sinh_fast (double *h, double *l, double x)
{
  /* magic is such that magic*x rounds to a number < 65535.5
     whatever the rounding mode */
  static const double magic = 0x1.70f77fc88ae3cp6;
  int k = __builtin_round (magic * x); /* k <= 65535 */
  /* |k - magic*x| <= 1/2 + |magic*x - round(magic*x)|
                   <= 1/2 + ulp(magic*x) <= 1/2 + 2^-37
     thus |x - k/magic| <= 0.00542055 */
  int i = k >> 8, j = k & 0xff;
  // if (x == TRACE) printf ("k=%d i=%d j=%d\n", k, i, j);
  double v = x - T[i][0];
  /* since x = T[i][0] + v, we approximate sinh(x) as
     sinh(T[i][0])*cosh(v) + cosh(T[i][0])*sinh(v)
     = T[i][1]*cosh(v) + T[i][2]*sinh(v) */
  double w = v - U[j][0];
  // if (x == TRACE) printf ("x=%la v=%la w=%la\n", x, v, w);
  /* since v = U[j][0] + w, we approximate sinh(v) as
     sinh(U[j][0])*cosh(w) + cosh(U[j][0])*sinh(w)
     = U[j][1]*cosh(w) + U[j][2]*sinh(w), and cosh(v) as
     sinh(U[j][0])*sinh(w) + cosh(U[j][0])*cosh(w)
     = U[j][1]*sinh(w) + U[j][2]*cosh(w) */

  /* since |T[i][0] - i*2^8/magic| < 2.36e-8 and
           |U[j][0] - j/magic| < 3.75e-9, we have:
     |x - T[i][0] - U[j][0]| < 0.00542055 + 2.36e-8 + 3.75e-9 < 0.00543 */

  // if (x == TRACE) printf ("T[i][0]=%la U[j][0]=%la w=%la\n", T[i][0], U[j][0], w);

  /* we have |w| < 0.00543 */
  double swh, swl, cwh, cwl;
  eval_S (h, l, w);
  // if (x == TRACE) printf ("swh=%la swl=%la\n", *h, *l);
  /* |h + l - sinh(w)| < 2^-67.99*|h| */

  if (k == 0)
    return __builtin_fma (0x1.02p-68, *h, 0x1p-1074);
  /* 2^-67.99 < 0x1.02p-68, and we add 2^-1074 to workaround cases
     when 0x1.02p-68 * h is rounded to zero */

  eval_C (&cwh, &cwl, w);
  // if (x == TRACE) printf ("cwh=%la cwl=%la\n", cwh, cwl);
  /* |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| */
  
  swh = *h;
  swl = *l;
  double svh, svl, cvh, cvl, h1, l1, h2, l2;
  s_mul (&h1, &l1, U[j][1], cwh, cwl); /* U[j][1]*cosh(w) */
  /* |U[j][1] - sinh(U[j][0])| < 2^-12 ulp(U[j][1]) <= 2^-64 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| thus
     |h1+l1-sinh(U[j][0])*cosh(w)| < 2^-63.91*|h1+l1| */
  // if (x == TRACE) printf ("h1=%la l1=%la\n", h1, l1);
  s_mul (&h2, &l2, U[j][2], swh, swl); /* U[j][1]*sinh(w) */
  // if (x == TRACE) printf ("h2=%la l2=%la\n", h2, l2);
  /* |U[j][2] - cosh(U[j][0])| < 2^-12 ulp(U[j][1]) <= 2^-64 |U[j][1]|
     and |swh + swl - sinh(w)| < 2^-67.99*|swh+swl| thus
     |h2+l2-cosh(U[j][0])*sinh(w)| < 2^-63.91*|h2+l2| */
  /* since h1+l1 and h2+l2 have the same relative error bound, that bounds
     holds for the sum of their absolute values, but we might have
     cancellation, the worst case being for j=1 and w=-0.00543,
     where h1+l1 >= 0.0108414. and h2+l2 >= -0.0054303,
     thus (|h1+l1| + |h2+l2|)/((|h1+l1| - |h2+l2|) < 3.008 */
     
  fast_two_sum (h, l, h1, h2);
  *l += l1 + l2; /* h+l approximates sinh(v) */
  /* we neglect the rounding errors in fast_two_sum and *l += l1 + l2,
     which would not change the error bound below, since the above
     bounds are upper approximations */
  // if (x == TRACE) printf ("h=%la l=%la\n", *h, *l);

  if (i == 0)
    return 0x1.9ap-63 * *h; /* 3.008*2^-63.91 < 0x1.9ap-63 */

  svh = *h;
  svl = *l;
  // if (x == TRACE) printf ("v=%la svh=%la svl=%la\n", v, svh, svl);
  s_mul (&h1, &l1, U[j][1], swh, swl); /* U[j][1]*sinh(w) */
  /* |U[j][1] - sinh(U[j][0])| < 2^-12 ulp(U[j][1]) <= 2^-64 |U[j][1]|
     and |swh + swl - sinh(w)| < 2^-67.99*|swh+swl| thus
     |h1+l1-sinh(U[j][0])*sinh(w)| < 2^-63.91*|h1+l1| */
  // if (x == TRACE) printf ("u=%la h1=%la l1=%la\n", U[j][0], h1, l1);
  s_mul (&h2, &l2, U[j][2], cwh, cwl); /* U[j][2]*cosh(w) */
  // if (x == TRACE) printf ("h2=%la l2=%la\n", h2, l2);
  /* |U[j][2] - cosh(U[j][0])| < 2^-12 ulp(U[j][1]) <= 2^-64 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| thus
     |h2+l2-cosh(U[j][0])*cosh(w)| < 2^-63.91*|h2+l2| */
  /* since h1+l1 and h2+l2 have the same relative error bound, that bounds
     holds for the sum of their absolute values, but we might have
     cancellation, the worst case being for j=1 and w=-0.00543,
     where h1+l1 >= 1.0000735. and h1+l1 >= -0.0000589
     thus (|h1+l1| + |h2+l2|)/((|h1+l1| - |h2+l2|) < 1.000118 */
  fast_two_sum (&cvh, &cvl, h2, h1); /* here |h2| > |h1| */
  cvl += l1 + l2; /* cvh+cvl approximates cosh(v) */
  /* we neglect the rounding errors in fast_two_sum and *l += l1 + l2,
     which would not change the error bound below, since the above
     bounds are upper approximations */
  // if (x == TRACE) printf ("cvh=%la cvl=%la\n", cvh, cvl);

  /* at this point cvh+cvl approximates cosh(v) with relative error bounded
     by 1.000118*2^-63.91 < 2^-63.90, svh+svl approximates sinh(v) with
     relative error bounded by 3.008*2^-63.91 < 2^-62.32, T[i][1] approximates
     sinh(T[i][0]) with relative error bounded by 2^-68, T[i][1]+T[i][2]
     approximates cosh(T[i][0]) with relative error bounded by 2^-68,
     and we have to compute:
     T[i][1]*(cvh+cvl) + (T[i][1]+T[i][2])*(svh+svl) =
     T[i][1]*(cvh+cvl+svh+svl) + T[i][2]*(svh+svl) */

  /* since |x - k/magic| <= 0.00542055, |T[i][0] - i*2^8/magic| < 2.36e-8,
     k = i*2^8+j:
     |v| = |x - T[i][0]|
       <= |x - k/magic| + |k/magic - i*2^8/magic| + |i*2^8/magic - T[i][0]|
       <= 0.00542055 + j/magic + 2.36e-8
       <= 0.00542055 + 255/magic + 2.36e-8 <= 2.77.
     We also have |v| >= 1/magic - (0.00542055 + 2.36e-8) > 0.00542.
     Thus |sinh(v)| < sinh(2.77) < 7.95 and |cosh(v)| < cosh(2.77) < 8.02,
     the absolute error on svh+svl is bounded by 3.008*2^-63.91*7.95
     < 2^-59.33, and the absolute error on cvh+cvl is bounded by
     1.000118*2^-63.91*8.02 < 2^-60.90. */

  fast_two_sum (&cvh, l, cvh, svh);
  cvl = *l + (cvl + svl);
  /* absolute error on cvh+cvl bounded by (neglecting the error in fast_two_sum
     and cvl = *l + (cvl + svl)):
     2^-59.33+2^-60.90 < 2^-58.91.
     Since |v| > 0.00542, the cancellation factor
     (cosh(v)+sinh(v))/(cosh(v)-sinh(v)) is bounded by 1.0109,
     thus the relative error on cvh+cvl is < 1.0109*2^-58.91 < 2^-58.89. */

  s_mul (&h1, &l1, T[i][1], cvh, cvl); /* T[i][1]*(cvh+cvl+svh+svl) */
  /* |T[i][1] - sinh(T[i][0])| < 2^-17 ulp(T[i][1]) <= 2^-69 |T[i][1]|
     and |cvh + cvl - (cosh(v)+sinh(v))| < 2^-58.89*|cvh + cvl| thus
     |h1+l1-sinh(T[i][0])*(cosh(v)+sinh(v))| < 2^-58.88*|h1+l1| */
  s_mul (&h2, &l2, T[i][2], svh, svl); /* T[i][2]*(svh+svl) */
  /* |T[i][2] - exp(T[i][0])| < 2^-17 ulp(T[i][2]) <= 2^-69 |T[i][2]|
     and |svh + svl - sinh(v)| < 2^-62.32*|svh + svl| thus
     |h2+l2-exp(T[i][0])*sinh(v)| < 2^-62.30*|h2+l2| */
  fast_two_sum (h, l, h1, h2);
  *l += l1 + l2;

  /* 2^-58.88 < 0x1.17p-59 and 2^-62.32 < 0x1.9bp-63 */
  return 0x1.17p-59 * h1 + 0x1.9bp-63 * (h2 > 0 ? h2 : -h2);
}

#define MASK 0x7fffffffffffffff /* to mask the sign bit */

double
cr_sinh (double x)
{
  d64u64 v = {.f = x};
  int e = (v.u >> 52) - 0x3ff;
  int s = v.u >> 63; /* sign bit */
  v.u &= (uint64_t) MASK; /* get absolute value */

  if (e == 0x400 || e == 0xc00 || v.f >= 0x1.633ce8fb9f87ep+9)
    /* NaN or overflow */
  {
    /* this will return NaN for x=NaN, Inf with the correct sign for x=+Inf
       or -Inf, and Inf/DBL_MAX or -Inf/-DBL_MAX for other |x| >= 2^10 */
    return x * 0x1p1023;
  }
  
  double h, l;
  double err = cr_sinh_fast (&h, &l, v.f);
  // if (x == TRACE) printf ("h=%la l=%la err=%la\n", h, l, err);
  double sign[] = { 1.0, -1.0 };
  h *= sign[s];
  l *= sign[s];
  
  double left  = h + (l - err);
  double right = h + (l + err);
  if (left == right)
    return left;

  return 0;
}
