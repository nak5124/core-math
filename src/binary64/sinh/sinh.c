/* Correctly rounded sinh for binary64 values.

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
   with accuracy >= 53+14 bits:
   |si - sinh(xi)| < 2^(-14-1) ulp(si), |ci - cosh(xi)| < 2^(-14-1) ulp(ci).
   We have |xi - i/magic| < 2.30e-07.
   Generated with ./buildu 14 with accompanying file buildu.c. */
static const double U[256][3] = {
   {0x0p+0, 0x0p+0, 0x1p+0}, /* i=0 */
   {0x1.633d9a6f0b004p-7, 0x1.633f62784d28p-7, 0x1.0003d9ea4b182p+0}, /* i=1 */
   {0x1.633d9a8bd1d79p-6, 0x1.6344bab917c09p-6, 0x1.000f67c6d8d91p+0}, /* i=2 */
   {0x1.0a6e34acf683p-5, 0x1.0a7a3b109eadp-5, 0x1.0022a9eed93cfp+0}, /* i=3 */
   {0x1.633d9b02cb4bdp-5, 0x1.635a1c3bb9e99p-5, 0x1.003da0f6336bcp+0}, /* i=4 */
   {0x1.bc0d00bf620d4p-5, 0x1.bc44ade32b1b9p-5, 0x1.00604dac869f7p+0}, /* i=5 */
   {0x1.0a6e350a8c8bp-4, 0x1.0a9e508d905fap-4, 0x1.008ab11eeee96p+0}, /* i=6 */
   {0x1.36d5e671c3771p-4, 0x1.37224c824a1f1p-4, 0x1.00bccc8eb3f6cp+0}, /* i=7 */
   {0x1.633d99380441bp-4, 0x1.63afa65604bbp-4, 0x1.00f6a18285549p+0}, /* i=8 */
   {0x1.8fa54dc51ba9p-4, 0x1.9047b59b3b0bcp-4, 0x1.013831b92f251p+0}, /* i=9 */
   {0x1.bc0d00f20d1a3p-4, 0x1.bcebcea367b7ep-4, 0x1.01817f2724798p+0}, /* i=10 */
   {0x1.e874b4b1b76bp-4, 0x1.e99d4b3c450dbp-4, 0x1.01d28c040336ep+0}, /* i=11 */
   {0x1.0a6e34c652ff4p-3, 0x1.0b2ec217dc345p-3, 0x1.022b5ac17e9e5p+0}, /* i=12 */
   {0x1.20a20d40df20ep-3, 0x1.2196e69123c74p-3, 0x1.028bee00e2863p+0}, /* i=13 */
   {0x1.36d5e77726522p-3, 0x1.3807c357a591ep-3, 0x1.02f448b9f454cp+0}, /* i=14 */
   {0x1.4d09c0e447a93p-3, 0x1.4e8202b9b2d92p-3, 0x1.03646e05258f1p+0}, /* i=15 */
   {0x1.633d9a5b2f6b7p-3, 0x1.650652ae3ac4dp-3, 0x1.03dc614607172p+0}, /* i=16 */
   {0x1.797174c5471b2p-3, 0x1.7b95619211b0dp-3, 0x1.045c261df504ap+0}, /* i=17 */
   {0x1.8fa54e324fd6cp-3, 0x1.922fdb2f21ad9p-3, 0x1.04e3c059e4ec5p+0}, /* i=18 */
   {0x1.a5d928074d2f3p-3, 0x1.a8d66f0a8dd88p-3, 0x1.057334168590fp+0}, /* i=19 */
   {0x1.bc0d015c8eaf7p-3, 0x1.bf89caaf334a2p-3, 0x1.060a859ee0738p+0}, /* i=20 */
   {0x1.d240db04ad09p-3, 0x1.d64a9dce02cb2p-3, 0x1.06a9b985ff718p+0}, /* i=21 */
   {0x1.e874b515da9c2p-3, 0x1.ed1997bfea48fp-3, 0x1.0750d49726c2dp+0}, /* i=22 */
   {0x1.fea88e881397cp-3, 0x1.01fbb391fcb96p-2, 0x1.07ffdbd19699ep+0}, /* i=23 */
   {0x1.0a6e3484a77c7p-2, 0x1.0d725eeb20cc5p-2, 0x1.08b6d486eedf9p+0}, /* i=24 */
   {0x1.1588209aa7ff5p-2, 0x1.18f124768fcd8p-2, 0x1.0975c41be16fap+0}, /* i=25 */
   {0x1.20a20cc654634p-2, 0x1.24785e0950413p-2, 0x1.0a3cb06375a65p+0}, /* i=26 */
   {0x1.2bbbfaad0523ep-2, 0x1.30086626c39e7p-2, 0x1.0b0b9f795270ep+0}, /* i=27 */
   {0x1.36d5e645fff22p-2, 0x1.3ba191aa09a3fp-2, 0x1.0be2974b8b8e7p+0}, /* i=28 */
   {0x1.41efd4926f6ffp-2, 0x1.47443f242e2fdp-2, 0x1.0cc19eb342ce6p+0}, /* i=29 */
   {0x1.4d09c11390d7p-2, 0x1.52f0c37f5655p-2, 0x1.0da8bc0d044bfp+0}, /* i=30 */
   {0x1.5823ad1e5ba9cp-2, 0x1.5ea77a09f90b3p-2, 0x1.0e97f66639b8ep+0}, /* i=31 */
   {0x1.633d9ad484133p-2, 0x1.6a68bf3f16adep-2, 0x1.0f8f5520cb8bep+0}, /* i=32 */
   {0x1.6e5787aac8c07p-2, 0x1.7634eaf8fe874p-2, 0x1.108edf76452f4p+0}, /* i=33 */
   {0x1.797174cddd6cap-2, 0x1.820c59552f84cp-2, 0x1.11969d3158fdep+0}, /* i=34 */
   {0x1.848b5fc227a5bp-2, 0x1.8def62e02f4d7p-2, 0x1.12a6960466b42p+0}, /* i=35 */
   {0x1.8fa54e19d6296p-2, 0x1.99de6921fab7ep-2, 0x1.13bed2a81eeb5p+0}, /* i=36 */
   {0x1.9abf3a3f273cap-2, 0x1.a5d9c206ad55dp-2, 0x1.14df5aff5bf73p+0}, /* i=37 */
   {0x1.a5d9273eee604p-2, 0x1.b1e1cd20f59b8p-2, 0x1.160838064ccbep+0}, /* i=38 */
   {0x1.b0f314604d197p-2, 0x1.bdf6e65615b31p-2, 0x1.1739729951f9ep+0}, /* i=39 */
   {0x1.bc0d01a423a2cp-2, 0x1.ca196ab7fa281p-2, 0x1.187313e7d4a24p+0}, /* i=40 */
   {0x1.c726ee1426721p-2, 0x1.d649b6aff6ce8p-2, 0x1.19b5254595f8dp+0}, /* i=41 */
   {0x1.d240dafb57c69p-2, 0x1.e288298ac8b0ep-2, 0x1.1affb08866fd3p+0}, /* i=42 */
   {0x1.dd5ac800192cep-2, 0x1.eed52135af3c4p-2, 0x1.1c52bf98867c6p+0}, /* i=43 */
   {0x1.e874b45b14378p-2, 0x1.fb30fb917e023p-2, 0x1.1dae5c9120a4p+0}, /* i=44 */
   {0x1.f38ea183274ecp-2, 0x1.03ce0cb8a6addp-1, 0x1.1f129215b3726p+0}, /* i=45 */
   {0x1.fea88e28e7ebp-2, 0x1.0a0b6c82c1acep-1, 0x1.207f6ab43a80dp+0}, /* i=46 */
   {0x1.04e13d938253fp-1, 0x1.1050cdb18701p-1, 0x1.21f4f182cf184p+0}, /* i=47 */
   {0x1.0a6e33f32e8e2p-1, 0x1.169e603d12f39p-1, 0x1.237331aac51c7p+0}, /* i=48 */
   {0x1.0ffb2a1b141bep-1, 0x1.1cf4549568f3ep-1, 0x1.24fa36a503b4p+0}, /* i=49 */
   {0x1.158822b8eadccp-1, 0x1.2352de9b8300bp-1, 0x1.268a0cf76641cp+0}, /* i=50 */
   {0x1.1b1517f85c89ap-1, 0x1.29ba28abec322p-1, 0x1.2822bf038d749p+0}, /* i=51 */
   {0x1.20a20e5b34bf4p-1, 0x1.302a694708036p-1, 0x1.29c45a5bf4ab2p+0}, /* i=52 */
   {0x1.262f03e111bc2p-1, 0x1.36a3cfb2b210fp-1, 0x1.2b6eeafb4539bp+0}, /* i=53 */
   {0x1.2bbbf9f463266p-1, 0x1.3d268f72efac2p-1, 0x1.2d227e23245f4p+0}, /* i=54 */
   {0x1.3148f074438bbp-1, 0x1.43b2da89b3e1ep-1, 0x1.2edf20e6dcd89p+0}, /* i=55 */
   {0x1.36d5e750aae26p-1, 0x1.4a48e355b8387p-1, 0x1.30a4e0a422e0fp+0}, /* i=56 */
   {0x1.3c62dce1be859p-1, 0x1.50e8da982bac4p-1, 0x1.3273ca7849a6bp+0}, /* i=57 */
   {0x1.41efd58eff412p-1, 0x1.5792f89eeba9p-1, 0x1.344bedc4be853p+0}, /* i=58 */
   {0x1.477cca35a487cp-1, 0x1.5e47682b42eb4p-1, 0x1.362d565cd5b4cp+0}, /* i=59 */
   {0x1.4d09c29ecd609p-1, 0x1.65066648fbb9ap-1, 0x1.3818155f39d3ap+0}, /* i=60 */
   {0x1.5296b73d51261p-1, 0x1.6bd01dc073899p-1, 0x1.3a0c36f41ca63p+0}, /* i=61 */
   {0x1.5823ad23ec1bap-1, 0x1.72a4c907541f8p-1, 0x1.3c09cbea7a9bcp+0}, /* i=62 */
   {0x1.5db0a4e2985e8p-1, 0x1.79849d752c1abp-1, 0x1.3e10e3d2ff4dep+0}, /* i=63 */
   {0x1.633d9aa6eddfap-1, 0x1.806fcb40d044cp-1, 0x1.40218ce62e7b4p+0}, /* i=64 */
   {0x1.68ca90c6740c3p-1, 0x1.87668a98839e9p-1, 0x1.423bd7e621fcdp+0}, /* i=65 */
   {0x1.6e57875c4aa11p-1, 0x1.8e69114469bf9p-1, 0x1.445fd511707dep+0}, /* i=66 */
   {0x1.73e47de3bb7c6p-1, 0x1.9577949c798f1p-1, 0x1.468d94b2d9cc3p+0}, /* i=67 */
   {0x1.7971741ea503fp-1, 0x1.9c924aab25b83p-1, 0x1.48c5277913942p+0}, /* i=68 */
   {0x1.7efe6b1e54387p-1, 0x1.a3b96b886135ap-1, 0x1.4b069ee6410b5p+0}, /* i=69 */
   {0x1.848b617db07ecp-1, 0x1.aaed2c80aa79dp-1, 0x1.4d520bc98e186p+0}, /* i=70 */
   {0x1.8a18576fcbe3p-1, 0x1.b22dc54c18d88p-1, 0x1.4fa77fe1df058p+0}, /* i=71 */
   {0x1.8fa54d96985dp-1, 0x1.b97b6e99da09dp-1, 0x1.52070d6c8d784p+0}, /* i=72 */
   {0x1.9532447cf6b5cp-1, 0x1.c0d661650ebfbp-1, 0x1.5470c6ef3f0cep+0}, /* i=73 */
   {0x1.9abf3a9ec54edp-1, 0x1.c83ed452a2f8p-1, 0x1.56e4be5729518p+0}, /* i=74 */
   {0x1.a04c320307fd9p-1, 0x1.cfb503287b017p-1, 0x1.5963076ff1414p+0}, /* i=75 */
   {0x1.a5d927b80c96p-1, 0x1.d7392368e187ep-1, 0x1.5bebb41a334cp+0}, /* i=76 */
   {0x1.ab661e4dcd2fbp-1, 0x1.decb726d989fep-1, 0x1.5e7ed903870e3p+0}, /* i=77 */
   {0x1.b0f3143cf5e5cp-1, 0x1.e66c28736592ap-1, 0x1.611c894c0ee58p+0}, /* i=78 */
   {0x1.b6800bbb04412p-1, 0x1.ee1b8347b1afep-1, 0x1.63c4da272bde9p+0}, /* i=79 */
   {0x1.bc0d0083b9e94p-1, 0x1.f5d9b82f0bb69p-1, 0x1.6677ddfeb0e3ep+0}, /* i=80 */
   {0x1.c199f826ee45bp-1, 0x1.fda70a9042719p-1, 0x1.6935ac5203a3bp+0}, /* i=81 */
   {0x1.c726edc456b4fp-1, 0x1.02c1d7d6291e9p+0, 0x1.6bfe57d8a4393p+0}, /* i=82 */
   {0x1.ccb3e6dde3caap-1, 0x1.06b7f5f270fe6p+0, 0x1.6ed1f8cb68502p+0}, /* i=83 */
   {0x1.d240dabda1318p-1, 0x1.0ab5f7e4ed4a5p+0, 0x1.71b0a07b32f34p+0}, /* i=84 */
   {0x1.d7cdd152c608bp-1, 0x1.0ebc021d2f2a3p+0, 0x1.749a69155b5e6p+0}, /* i=85 */
   {0x1.dd5ac752265aap-1, 0x1.12ca313a6bde6p+0, 0x1.778f6756135ap+0}, /* i=86 */
   {0x1.e2e7bf50cc625p-1, 0x1.16e0a65f5b5a5p+0, 0x1.7a8fb3688d299p+0}, /* i=87 */
   {0x1.e874b5df6fb7ep-1, 0x1.1aff7e81fff22p+0, 0x1.7d9b6290db0d2p+0}, /* i=88 */
   {0x1.ee01ab8c1d302p-1, 0x1.1f26d9c36b2e9p+0, 0x1.80b28c8b43958p+0}, /* i=89 */
   {0x1.f38ea22c82baep-1, 0x1.2356d981ba455p+0, 0x1.83d54a2b6bb91p+0}, /* i=90 */
   {0x1.f91b9771f089ep-1, 0x1.278f9c3f56ca6p+0, 0x1.8703b24968d81p+0}, /* i=91 */
   {0x1.fea88e813aa91p-1, 0x1.2bd144e5f8dc1p+0, 0x1.8a3ddf34d646ep+0}, /* i=92 */
   {0x1.021ac24954c52p+0, 0x1.301bf21ecbcefp+0, 0x1.8d83e82dff514p+0}, /* i=93 */
   {0x1.04e13d9c92e68p+0, 0x1.346fc62edf0c8p+0, 0x1.90d5e7582573bp+0}, /* i=94 */
   {0x1.07a7b811ae1bp+0, 0x1.38cce099b598ep+0, 0x1.9433f4dfa7b51p+0}, /* i=95 */
   {0x1.0a6e349fabf1fp+0, 0x1.3d3367ab7a06dp+0, 0x1.979e2e5584443p+0}, /* i=96 */
   {0x1.0d34aedcc3078p+0, 0x1.41a37648a1a5bp+0, 0x1.9b14a897bc0a7p+0}, /* i=97 */
   {0x1.0ffb2a42eb99p+0, 0x1.461d3428c320bp+0, 0x1.9e9782a12dc4ap+0}, /* i=98 */
   {0x1.12c1a4fa438e8p+0, 0x1.4aa0c0ce5e191p+0, 0x1.a226d52b5f5a5p+0}, /* i=99 */
   {0x1.158821891b2dbp+0, 0x1.4f2e431f7c5f9p+0, 0x1.a5c2bee6acd49p+0}, /* i=100 */
   {0x1.184e9c4cc620bp+0, 0x1.53c5d836e51bp+0, 0x1.a96b56e587ae3p+0}, /* i=101 */
   {0x1.1b15175448e0cp+0, 0x1.5867a6d32c0bfp+0, 0x1.ad20bc024850cp+0}, /* i=102 */
   {0x1.1ddb927f0cadap+0, 0x1.5d13d26dc6a1ap+0, 0x1.b0e30aa4d7b91p+0}, /* i=103 */
   {0x1.20a20c6aa517dp+0, 0x1.61ca7ca91ed28p+0, 0x1.b4b25dd975f4cp+0}, /* i=104 */
   {0x1.2368888c49d68p+0, 0x1.668bcfb883326p+0, 0x1.b88ed7be9d421p+0}, /* i=105 */
   {0x1.262f045c66d9bp+0, 0x1.6b57ebf1b8ee1p+0, 0x1.bc78929a82b57p+0}, /* i=106 */
   {0x1.28f57fb629cdep+0, 0x1.702ef6040c0e7p+0, 0x1.c06fac57727b1p+0}, /* i=107 */
   {0x1.2bbbfae83fc8ep+0, 0x1.751113bcb54fcp+0, 0x1.c47443ec26c4fp+0}, /* i=108 */
   {0x1.2e8275d99a639p+0, 0x1.79fe6a8a20455p+0, 0x1.c88678240c2bep+0}, /* i=109 */
   {0x1.3148f0ac8a94fp+0, 0x1.7ef7209a59a2ep+0, 0x1.cca6688a3f4ffp+0}, /* i=110 */
   {0x1.340f6c2adca41p+0, 0x1.83fb5da54a9f7p+0, 0x1.d0d436141c7ccp+0}, /* i=111 */
   {0x1.36d5e74c4479p+0, 0x1.890b4672e8642p+0, 0x1.d50fff66131bcp+0}, /* i=112 */
   {0x1.399c625e184aep+0, 0x1.8e27028b45aeap+0, 0x1.d959e58eb91c8p+0}, /* i=113 */
   {0x1.3c62dd7d54b1p+0, 0x1.934eb97c7422dp+0, 0x1.ddb209c34dcbfp+0}, /* i=114 */
   {0x1.3f295899d0e98p+0, 0x1.988292ddbf3adp+0, 0x1.e2188d610f991p+0}, /* i=115 */
   {0x1.41efd3bc6ccf9p+0, 0x1.9dc2b6d2626ep+0, 0x1.e68d925b17f77p+0}, /* i=116 */
   {0x1.44b64eaa1de81p+0, 0x1.a30f4d5a9d7e7p+0, 0x1.eb113aa5c8c16p+0}, /* i=117 */
   {0x1.477cca90e99d2p+0, 0x1.a868818d6aa28p+0, 0x1.efa3aaf7ee5c5p+0}, /* i=118 */
   {0x1.4a4345a0b34b6p+0, 0x1.adce791f802f3p+0, 0x1.f445038d43f2bp+0}, /* i=119 */
   {0x1.4d09bf422e1b5p+0, 0x1.b3415c6f4b72bp+0, 0x1.f8f5670102c5ep+0}, /* i=120 */
   {0x1.4fd03c27e95fp+0, 0x1.b8c15ebae7504p+0, 0x1.fdb5017002eb7p+0}, /* i=121 */
   {0x1.5296b7d08e393p+0, 0x1.be4ea18030985p+0, 0x1.0141f7e461c68p+1}, /* i=122 */
   {0x1.555d323f8e424p+0, 0x1.c3e94f778902dp+0, 0x1.03b12b82fb6b9p+1}, /* i=123 */
   {0x1.5823ad0b8ef69p+0, 0x1.c99196fa5e576p+0, 0x1.06282fb51189bp+1}, /* i=124 */
   {0x1.5aea29253a95ap+0, 0x1.cf47a592edcb1p+0, 0x1.08a71853d8b76p+1}, /* i=125 */
   {0x1.5db0a41c2a7a4p+0, 0x1.d50ba23af286ap+0, 0x1.0b2df667af392p+1}, /* i=126 */
   {0x1.60771f6b5dacap+0, 0x1.daddbc6893f46p+0, 0x1.0dbcdebcba6d6p+1}, /* i=127 */
   {0x1.633d9a84779aep+0, 0x1.e0be1fc7b32dbp+0, 0x1.1053e484b96c3p+1}, /* i=128 */
   {0x1.660415d058621p+0, 0x1.e6acfa79c93d4p+0, 0x1.12f31c1466062p+1}, /* i=129 */
   {0x1.68ca910cae844p+0, 0x1.ecaa79a50e8b6p+0, 0x1.159a995fbb397p+1}, /* i=130 */
   {0x1.6b910b91e7a06p+0, 0x1.f2b6c9fd4e6ecp+0, 0x1.184a7034b286fp+1}, /* i=131 */
   {0x1.6e578721042cbp+0, 0x1.f8d21ded7e9c7p+0, 0x1.1b02b6f987583p+1}, /* i=132 */
   {0x1.711e02af0b799p+0, 0x1.fefca23c1235bp+0, 0x1.1dc381a05c484p+1}, /* i=133 */
   {0x1.73e47e0880078p+0, 0x1.029b42f9b3d94p+1, 0x1.208ce529ca781p+1}, /* i=134 */
   {0x1.76aaf8a2ab691p+0, 0x1.05bffbe49990dp+1, 0x1.235ef67ceeda8p+1}, /* i=135 */
   {0x1.797174109567cp+0, 0x1.08ec95ddeaf36p+1, 0x1.2639cced7b822p+1}, /* i=136 */
   {0x1.7c37f15f629d9p+0, 0x1.0c212a98283d3p+1, 0x1.291d7f9a7985p+1}, /* i=137 */
   {0x1.7efe69ca938e5p+0, 0x1.0f5dcae5b408dp+1, 0x1.2c0a1dab1c2ffp+1}, /* i=138 */
   {0x1.81c4e69bd9ff9p+0, 0x1.12a29a88251b1p+1, 0x1.2effc76cb506p+1}, /* i=139 */
   {0x1.848b6325176e2p+0, 0x1.15efad3c0f2e5p+1, 0x1.31fe8ec1981c6p+1}, /* i=140 */
   {0x1.8751dc303f355p+0, 0x1.1945188bd5d35p+1, 0x1.35068731f5d02p+1}, /* i=141 */
   {0x1.8a1857876e063p+0, 0x1.1ca2fd13d90e7p+1, 0x1.3817ce63852e4p+1}, /* i=142 */
   {0x1.8cded2771f5abp+0, 0x1.2009718484b8fp+1, 0x1.3b327904009d2p+1}, /* i=143 */
   {0x1.8fa54db4ff135p+0, 0x1.237890ec9e3f5p+1, 0x1.3e569fc6ff481p+1}, /* i=144 */
   {0x1.926bc894050bep+0, 0x1.26f074e91f42cp+1, 0x1.41845a1a3f9b8p+1}, /* i=145 */
   {0x1.953244211d893p+0, 0x1.2a7139828aaf2p+1, 0x1.44bbc1af17faep+1}, /* i=146 */
   {0x1.97f8c030d80fap+0, 0x1.2dfaf98321b5ap+1, 0x1.47fcef1fb45f8p+1}, /* i=147 */
   {0x1.9abf3c99be63bp+0, 0x1.318dcffaca8a1p+1, 0x1.4b47fb5196c34p+1}, /* i=148 */
   {0x1.9d85b616d38ap+0, 0x1.3529d42d8bap+1, 0x1.4e9cfbb34dfddp+1}, /* i=149 */
   {0x1.a04c314964588p+0, 0x1.38cf27e827ca6p+1, 0x1.51fc0f77a364p+1}, /* i=150 */
   {0x1.a312acd2acb44p+0, 0x1.3c7de57d0c7ap+1, 0x1.55654ef7733fbp+1}, /* i=151 */
   {0x1.a5d927b0d2725p+0, 0x1.403627f1cccdfp+1, 0x1.58d8d33a23ab1p+1}, /* i=152 */
   {0x1.a89fa2ed12739p+0, 0x1.43f80d4f1cdc4p+1, 0x1.5c56b81d22235p+1}, /* i=153 */
   {0x1.ab661e49a7e19p+0, 0x1.47c3b234e8076p+1, 0x1.5fdf183aade8ap+1}, /* i=154 */
   {0x1.ae2c98bdad4c6p+0, 0x1.4b993270c989bp+1, 0x1.63720d73060bep+1}, /* i=155 */
   {0x1.b0f31409ed2aap+0, 0x1.4f78adf6647eep+1, 0x1.670fb58f32dfp+1}, /* i=156 */
   {0x1.b3b98f8d3e99ap+0, 0x1.536241bf681cbp+1, 0x1.6ab82b9d35a8p+1}, /* i=157 */
   {0x1.b680080bd6e86p+0, 0x1.5756074ffc709p+1, 0x1.6e6b877676cd6p+1}, /* i=158 */
   {0x1.b9468664f02b3p+0, 0x1.5b5429cb9b2p+1, 0x1.7229f181715c1p+1}, /* i=159 */
   {0x1.bc0d015d3faeep+0, 0x1.5f5cbab9114fp+1, 0x1.75f37a2f073bbp+1}, /* i=160 */
   {0x1.bed37b61a6a9ep+0, 0x1.636fdc9782725p+1, 0x1.79c841e02acb3p+1}, /* i=161 */
   {0x1.c199f6de3162ap+0, 0x1.678db25af6521p+1, 0x1.7da86972655ecp+1}, /* i=162 */
   {0x1.c460733bd971fp+0, 0x1.6bb65ae630c39p+1, 0x1.81940dfb59432p+1}, /* i=163 */
   {0x1.c726ed59f7ac1p+0, 0x1.6fe9f189d0686p+1, 0x1.858b4939f2525p+1}, /* i=164 */
   {0x1.c9ed692eb45b3p+0, 0x1.74289c9ec1968p+1, 0x1.898e3f5fdd7a6p+1}, /* i=165 */
   {0x1.ccb3e411b470ep+0, 0x1.787278cb00dcap+1, 0x1.8d9d0b80a65f4p+1}, /* i=166 */
   {0x1.cf7a5fb981325p+0, 0x1.7cc7a9bb6977bp+1, 0x1.91b7cf5db5a55p+1}, /* i=167 */
   {0x1.d240db51193b6p+0, 0x1.81284f874a087p+1, 0x1.95dea95e15a0bp+1}, /* i=168 */
   {0x1.d5075667d41b8p+0, 0x1.85948b312b9ccp+1, 0x1.9a11b8d01c1a2p+1}, /* i=169 */
   {0x1.d7cdd1f3a32d6p+0, 0x1.8a0c8051b99d4p+1, 0x1.9e511f804693ep+1}, /* i=170 */
   {0x1.da944cc155361p+0, 0x1.8e904f63c7f03p+1, 0x1.a29cfc4e258efp+1}, /* i=171 */
   {0x1.dd5ac7cab4322p+0, 0x1.93201cc0ff8c1p+1, 0x1.a6f571d09efbp+1}, /* i=172 */
   {0x1.e02142e5dd656p+0, 0x1.97bc0b496e072p+1, 0x1.ab5aa140d5beep+1}, /* i=173 */
   {0x1.e2e7bc4344e63p+0, 0x1.9c643b707d979p+1, 0x1.afcca98fb555dp+1}, /* i=174 */
   {0x1.e5ae37622082ep+0, 0x1.a118d6f33eb4ep+1, 0x1.b44bb294f0e3dp+1}, /* i=175 */
   {0x1.e874b50c1d9fep+0, 0x1.a5da03831ba06p+1, 0x1.b8d7e0577368bp+1}, /* i=176 */
   {0x1.eb3b3016372bbp+0, 0x1.aaa7dce43428cp+1, 0x1.bd714d65d3e27p+1}, /* i=177 */
   {0x1.ee01a9dda5dadp+0, 0x1.af828a59294e2p+1, 0x1.c2181f53f326cp+1}, /* i=178 */
   {0x1.f0c826f66d997p+0, 0x1.b46a395612348p+1, 0x1.c6cc81b18d60ep+1}, /* i=179 */
   {0x1.f38ea0999ea5ep+0, 0x1.b95f03a75b482p+1, 0x1.cb8e8d3d28fd2p+1}, /* i=180 */
   {0x1.f6551c6975ca6p+0, 0x1.be6119806970cp+1, 0x1.d05e703be7bcdp+1}, /* i=181 */
   {0x1.f91b963605abfp+0, 0x1.c37099e71e421p+1, 0x1.d53c487cce6b9p+1}, /* i=182 */
   {0x1.fbe213bca0067p+0, 0x1.c88db659ecbb3p+1, 0x1.da284599a150cp+1}, /* i=183 */
   {0x1.fea88df97bf9p+0, 0x1.cdb889527831dp+1, 0x1.df22810df1203p+1}, /* i=184 */
   {0x1.00b785dd10f72p+1, 0x1.d2f1474a47642p+1, 0x1.e42b2d65ef6c6p+1}, /* i=185 */
   {0x1.021ac25d1ef66p+1, 0x1.d83808e24427ap+1, 0x1.e942625f1dc82p+1}, /* i=186 */
   {0x1.037dfe6f26636p+1, 0x1.dd8cfdfa3f40dp+1, 0x1.ee684e26ccb32p+1}, /* i=187 */
   {0x1.04e13df1efa17p+1, 0x1.e2f05eb161312p+1, 0x1.f39d26f050449p+1}, /* i=188 */
   {0x1.06447c012d64ap+1, 0x1.e86241a8a7bf8p+1, 0x1.f8e1029b5d8e5p+1}, /* i=189 */
   {0x1.07a7b9219ffefp+1, 0x1.ede2d2b3d5001p+1, 0x1.fe340b86b0359p+1}, /* i=190 */
   {0x1.090af607bd6c6p+1, 0x1.f3723eed54353p+1, 0x1.01cb36ab1f7cp+2}, /* i=191 */
   {0x1.0a6e32e5fcfedp+1, 0x1.f910b1eee68f4p+1, 0x1.048429212b4fep+2}, /* i=192 */
   {0x1.0bd171a7c237ep+1, 0x1.febe5ee4557ddp+1, 0x1.0744f5f173e58p+2}, /* i=193 */
   {0x1.0d34aed4e9f14p+1, 0x1.023db1ae9145dp+2, 0x1.0a0dab7128eap+2}, /* i=194 */
   {0x1.0e97ecbb7e9f9p+1, 0x1.0523fa831feb4p+2, 0x1.0cde63aa6bccbp+2}, /* i=195 */
   {0x1.0ffb2b3d436fp+1, 0x1.0812201157a5ap+2, 0x1.0fb7341acf047p+2}, /* i=196 */
   {0x1.115e67d11ce25p+1, 0x1.0b0833852848ep+2, 0x1.12982d7098bc2p+2}, /* i=197 */
   {0x1.12c1a66c9c096p+1, 0x1.0e065421965dcp+2, 0x1.15816e1013d07p+2}, /* i=198 */
   {0x1.1424e2277b96fp+1, 0x1.110c8e56c8f06p+2, 0x1.187302123c16bp+2}, /* i=199 */
   {0x1.158820c75ceep+1, 0x1.141b060500d77p+2, 0x1.1b6d0c6261145p+2}, /* i=200 */
   {0x1.16eb5e7cd41f3p+1, 0x1.1731ca643949ap+2, 0x1.1e6f9bd3b7da1p+2}, /* i=201 */
   {0x1.184e9c729b846p+1, 0x1.1a50f5d25e14ap+2, 0x1.217aca1760fdp+2}, /* i=202 */
   {0x1.19b1dae0ff6abp+1, 0x1.1d78a0e051f75p+2, 0x1.248eaf2148676p+2}, /* i=203 */
   {0x1.1b1517caef957p+1, 0x1.20a8df4c08f6dp+2, 0x1.27ab5e332c21ep+2}, /* i=204 */
   {0x1.1c78553e32347p+1, 0x1.23e1ce53972d7p+2, 0x1.2ad0f3d74a3abp+2}, /* i=205 */
   {0x1.1ddb909a1ae42p+1, 0x1.272380a0092a2p+2, 0x1.2dff8245f0a3fp+2}, /* i=206 */
   {0x1.1f3ed0ff5a611p+1, 0x1.2a6e201d9c828p+2, 0x1.3137327646f57p+2}, /* i=207 */
   {0x1.20a20e33237e3p+1, 0x1.2dc1b2b302bd1p+2, 0x1.34780a2ddb4dbp+2}, /* i=208 */
   {0x1.22054ab2c5f21p+1, 0x1.311e57d14ad38p+2, 0x1.37c2282b0276ap+2}, /* i=209 */
   {0x1.236887f3161fep+1, 0x1.34842ce53f769p+2, 0x1.3b15a938624a4p+2}, /* i=210 */
   {0x1.24cbc5d6d5733p+1, 0x1.37f34be0f1061p+2, 0x1.3e72a6bcae166p+2}, /* i=211 */
   {0x1.262f03bb3e1ebp+1, 0x1.3b6bcdab69dcep+2, 0x1.41d9391bbb221p+2}, /* i=212 */
   {0x1.27924159e429dp+1, 0x1.3eedcc4d7d297p+2, 0x1.454979d7b6ec7p+2}, /* i=213 */
   {0x1.28f57f256b7e3p+1, 0x1.427963ecc49b9p+2, 0x1.48c38487d43eap+2}, /* i=214 */
   {0x1.2a58bc6b08b4p+1, 0x1.460eae0b7eaap+2, 0x1.4c47722fe60bap+2}, /* i=215 */
   {0x1.2bbbfb166f167p+1, 0x1.49adcb420bd0dp+2, 0x1.4fd562cab8c29p+2}, /* i=216 */
   {0x1.2d1f382ea6d2ep+1, 0x1.4d56cfb18a874p+2, 0x1.536d6a19e35d1p+2}, /* i=217 */
   {0x1.2e8275d236b97p+1, 0x1.5109dd1756c31p+2, 0x1.570fa93dbfc5cp+2}, /* i=218 */
   {0x1.2fe5b33be487ep+1, 0x1.54c70de7af432p+2, 0x1.5abc3a327faa1p+2}, /* i=219 */
   {0x1.3148f023a6bbep+1, 0x1.588e7e250eea5p+2, 0x1.5e73387dcc053p+2}, /* i=220 */
   {0x1.32ac2dde293e1p+1, 0x1.5c604e8e68ef6p+2, 0x1.6234c45037596p+2}, /* i=221 */
   {0x1.340f6c3a6b2c3p+1, 0x1.603c9c18d51cep+2, 0x1.6600fa2411541p+2}, /* i=222 */
   {0x1.3572a9d7666c7p+1, 0x1.642380a9408aap+2, 0x1.69d7f3724b827p+2}, /* i=223 */
   {0x1.36d5e767ad2adp+1, 0x1.68151c3a9c12dp+2, 0x1.6db9cfb3125c8p+2}, /* i=224 */
   {0x1.383924e2f5214p+1, 0x1.6c118d14443efp+2, 0x1.71a6acb4cf446p+2}, /* i=225 */
   {0x1.399c621b78e6dp+1, 0x1.7018f162936fap+2, 0x1.759ea82dd254fp+2}, /* i=226 */
   {0x1.3affa08f3fc66p+1, 0x1.742b6c8e99087p+2, 0x1.79a1e5001448ap+2}, /* i=227 */
   {0x1.3c62dd7ccac4ep+1, 0x1.784915dc0d362p+2, 0x1.7db07a181c7c6p+2}, /* i=228 */
   {0x1.3dc61b5627925p+1, 0x1.7c7214386b082p+2, 0x1.81ca8dd492233p+2}, /* i=229 */
   {0x1.3f295953c0049p+1, 0x1.80a6856577ac9p+2, 0x1.85f03d8ce6b52p+2}, /* i=230 */
   {0x1.408c9701ca0c6p+1, 0x1.84e6886583b36p+2, 0x1.8a21a7d68ab45p+2}, /* i=231 */
   {0x1.41efd4028b69ap+1, 0x1.89323cc96b20bp+2, 0x1.8e5eebd59afcp+2}, /* i=232 */
   {0x1.435312322918ep+1, 0x1.8d89c9701e5a8p+2, 0x1.92a82fe68bcfcp+2}, /* i=233 */
   {0x1.44b65018c33d4p+1, 0x1.91ed4b3eaec19p+2, 0x1.96fd908fad69p+2}, /* i=234 */
   {0x1.46198afe18c5dp+1, 0x1.965cdb3ddf20ap+2, 0x1.9b5f2688144bbp+2}, /* i=235 */
   {0x1.477cca9bd1a24p+1, 0x1.9ad8b4611bdd1p+2, 0x1.9fcd2c09fc7fbp+2}, /* i=236 */
   {0x1.48e008ed9cc19p+1, 0x1.9f60e5edf8922p+2, 0x1.a447b02cfbf34p+2}, /* i=237 */
   {0x1.4a43451bdb06fp+1, 0x1.a3f58fde7effp+2, 0x1.a8ced28983037p+2}, /* i=238 */
   {0x1.4ba6833ba59e7p+1, 0x1.a896e2f28d89ep+2, 0x1.ad62c34f63cd1p+2}, /* i=239 */
   {0x1.4d09c0f47afd4p+1, 0x1.ad44fb156485bp+2, 0x1.b2039e19a5d66p+2}, /* i=240 */
   {0x1.4e6cfe8b2083bp+1, 0x1.b1fffd33d8628p+2, 0x1.b6b1876c7dad2p+2}, /* i=241 */
   {0x1.4fd03b80f7ac6p+1, 0x1.b6c80c0343573p+2, 0x1.bb6ca19cf2073p+2}, /* i=242 */
   {0x1.513379a1812abp+1, 0x1.bb9d528ef212p+2, 0x1.c0351741d81b1p+2}, /* i=243 */
   {0x1.5296b7b80f8fcp+1, 0x1.c07ff1ee06455p+2, 0x1.c50b091ae8dc6p+2}, /* i=244 */
   {0x1.53f9f4b49a8e3p+1, 0x1.c5700bf31f7ap+2, 0x1.c9ee98a2a7482p+2}, /* i=245 */
   {0x1.555d326c961cbp+1, 0x1.ca6dcd28210bcp+2, 0x1.cedff1f1bd3b8p+2}, /* i=246 */
   {0x1.56c06f428bc7cp+1, 0x1.cf79562d54d49p+2, 0x1.d3df355751889p+2}, /* i=247 */
   {0x1.5823acf21a17dp+1, 0x1.d492d42b23ab3p+2, 0x1.d8ec8f8df23cap+2}, /* i=248 */
   {0x1.5986ec789f781p+1, 0x1.d9ba7235145d2p+2, 0x1.de082b435dac5p+2}, /* i=249 */
   {0x1.5aea27b72ad74p+1, 0x1.def0411d1d4d4p+2, 0x1.e3321921ab733p+2}, /* i=250 */
   {0x1.5c4d67081d394p+1, 0x1.e434885f32266p+2, 0x1.e86aa001c51bap+2}, /* i=251 */
   {0x1.5db0a47657d24p+1, 0x1.e9875a11f5d89p+2, 0x1.edb1d1d28ad56p+2}, /* i=252 */
   {0x1.5f13e18f7e376p+1, 0x1.eee8e508fae79p+2, 0x1.f307dd00b4c35p+2}, /* i=253 */
   {0x1.60771f92133bap+1, 0x1.f45957932a221p+2, 0x1.f86cef785dca5p+2}, /* i=254 */
   {0x1.61da5cdc19ac9p+1, 0x1.f9d8d5325e4ccp+2, 0x1.fde12c712bbddp+2}, /* i=255 */
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

// Add a + (bh + bl), assuming |a| >= |bh|
static inline void fast_sum(double *hi, double *lo, double a, double bh,
                            double bl) {
  fast_two_sum(hi, lo, a, bh);
  /* |(a+bh)-(hi+lo)| <= 2^-105 |hi| and |lo| < ulp(hi) */
  *lo += bl;
  /* |(a+bh+bl)-(hi+lo)| <= 2^-105 |hi| + ulp(lo),
     where |lo| <= ulp(hi) + |bl|. */
}

// Add (ah + al) + (bh + bl), assuming |ah| >= |bh|
static inline void fast_sum2(double *hi, double *lo, double ah, double al,
                             double bh, double bl) {
  fast_two_sum (hi, lo, ah, bh);
  *lo += al + bl;
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

/* The following is a degree-9 odd polynomial approximating sinh(x)
   for |x| < 0.00543 generated by Sollya (see Psinh2.sollya), with
   double-double coefficients and maximal relative error 2^-108.33.
   The code below assumes that the degree-1 coefficient is 1. */
static const double S2[][2] = {
  {0x1p0, 0},                                     /* degree 1 */
  {0x1.5555555555555p-3, 0x1.55555554062b9p-57},  /* degree 3 */
  {0x1.1111111111111p-7, 0x1.126bf9abf837p-63},   /* degree 5 */
  {0x1.a01a01a01989fp-13, 0},                     /* degree 7 */
  {0x1.71de4b3a00401p-19, 0},                     /* degree 9 */
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

/* The following is a degree-8 even polynomial approximating cosh(x)
   for |x| < 0.00543 generated by Sollya (see Pcosh2.sollya), with
   double-double coefficients and maximal absolute/relative error 2^-105.803.
   The code below assumes that the constant coefficient is 1. */
static const double C2[][2] = {
  {0x1p0, 0},                                     /* degree 0 */
  {0x1p-1, -0x1.27726p-86},                       /* degree 2 */
  {0x1.5555555555555p-5, 0x1.560cce697b2a2p-59},  /* degree 4 */
  {0x1.6c16c16c1633p-10, 0},                      /* degree 6 */
  {0x1.a01a1776b8d0bp-16, 0},                     /* degree 8 */
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

/* put in h+l a double-double approximation of sinh(w), for |w| < 0.00543 */
static void
eval_S2 (double *h, double *l, double w)
{
  double zh, zl;
  a_mul (&zh, &zl, w, w); /* zh+zl = w^2 */
  *h = __builtin_fma (S2[4][0], zh, S2[3][0]);
  if (w == TRACE) printf ("h3=%la\n", *h);
  /* We neglect at input S2[4][0]*zl*w^7 which relatively to sinh(w) ~ w
     is less than S2[4][0]*zl*w^7/w < 2^-18.46*2^-68*2^-45.14 = 2^-131.6.
     We have |h| < 2^-12.29, we neglect in output ulp(h)*w^7, which relatively
     to w gives ulp(2^-12.29)*w^6 < 2^-110.14 */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  fast_sum2 (h, l, S2[2][0], S2[2][1], *h, *l); /* add S2[2] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  fast_sum2 (h, l, S2[1][0], S2[1][1], *h, *l); /* add S2[1] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  s_mul (h, l, w, *h, *l);                      /* multiply by w */
  fast_sum (h, l, w, *h, *l);                   /* add w */
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

/* put in h+l a double-double approximation of cosh(w), for |w| < 0.00543 */
static void
eval_C2 (double *h, double *l, double w)
{
  double zh, zl;
  a_mul (&zh, &zl, w, w); /* zh+zl = w^2 */
  *h = __builtin_fma (C2[4][0], zh, C2[3][0]);
  s_mul (h, l, *h, zh, zl);                     /* multiply by w^2 */
  fast_sum (h, l, C2[3][0], *h, *l);            /* add C2[3] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  fast_sum2 (h, l, C2[2][0], C2[2][1], *h, *l); /* add C2[2] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  fast_sum2 (h, l, C2[1][0], C2[1][1], *h, *l); /* add C2[1] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  fast_sum (h, l, 1.0, *h, *l);                 /* add 1 */
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
           |U[j][0] - j/magic| < 1.92e-8, we have:
     |x - T[i][0] - U[j][0]| < 0.00542055 + 2.36e-8 + 1.92e-8 < 0.00543 */

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
  /* |U[j][1] - sinh(U[j][0])| < 2^-13 ulp(U[j][1]) <= 2^-65 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| thus
     |h1+l1-sinh(U[j][0])*cosh(w)| < 2^-64.82*|h1+l1| */
  // if (x == TRACE) printf ("h1=%la l1=%la\n", h1, l1);
  s_mul (&h2, &l2, U[j][2], swh, swl); /* U[j][1]*sinh(w) */
  // if (x == TRACE) printf ("h2=%la l2=%la\n", h2, l2);
  /* |U[j][2] - cosh(U[j][0])| < 2^-13 ulp(U[j][1]) <= 2^-65 |U[j][1]|
     and |swh + swl - sinh(w)| < 2^-67.99*|swh+swl| thus
     |h2+l2-cosh(U[j][0])*sinh(w)| < 2^-64.82*|h2+l2| */
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
    return 0x1.b5p-64 * *h; /* 3.008*2^-64.82 < 0x1.b5p-64 */

  svh = *h;
  svl = *l;
  // if (x == TRACE) printf ("v=%la svh=%la svl=%la\n", v, svh, svl);
  s_mul (&h1, &l1, U[j][1], swh, swl); /* U[j][1]*sinh(w) */
  /* |U[j][1] - sinh(U[j][0])| < 2^-13 ulp(U[j][1]) <= 2^-65 |U[j][1]|
     and |swh + swl - sinh(w)| < 2^-67.99*|swh+swl| thus
     |h1+l1-sinh(U[j][0])*sinh(w)| < 2^-64.82*|h1+l1| */
  // if (x == TRACE) printf ("u=%la h1=%la l1=%la\n", U[j][0], h1, l1);
  s_mul (&h2, &l2, U[j][2], cwh, cwl); /* U[j][2]*cosh(w) */
  // if (x == TRACE) printf ("h2=%la l2=%la\n", h2, l2);
  /* |U[j][2] - cosh(U[j][0])| < 2^-13 ulp(U[j][1]) <= 2^-65 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| thus
     |h2+l2-cosh(U[j][0])*cosh(w)| < 2^-64.82*|h2+l2| */
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
     by 1.000118*2^-64.82 < 2^-64.81, svh+svl approximates sinh(v) with
     relative error bounded by 3.008*2^-64.82 < 2^-63.23, T[i][1] approximates
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
     the absolute error on svh+svl is bounded by 3.008*2^-64.82*7.95
     < 2^-60.24, and the absolute error on cvh+cvl is bounded by
     1.000118*2^-64.82*8.02 < 2^-61.81. */

  fast_two_sum (&cvh, l, cvh, svh);
  cvl = *l + (cvl + svl);
  /* absolute error on cvh+cvl bounded by (neglecting the error in fast_two_sum
     and cvl = *l + (cvl + svl)):
     2^-60.24+2^-61.81 < 2^-59.82.
     Since |v| > 0.00542, the cancellation factor
     (cosh(v)+sinh(v))/(cosh(v)-sinh(v)) is bounded by 1.0109,
     thus the relative error on cvh+cvl is < 1.0109*2^-59.82 < 2^-59.80. */

  s_mul (&h1, &l1, T[i][1], cvh, cvl); /* T[i][1]*(cvh+cvl+svh+svl) */
  /* |T[i][1] - sinh(T[i][0])| < 2^-17 ulp(T[i][1]) <= 2^-69 |T[i][1]|
     and |cvh + cvl - (cosh(v)+sinh(v))| < 2^-59.80*|cvh + cvl| thus
     |h1+l1-sinh(T[i][0])*(cosh(v)+sinh(v))| < 2^-59.79*|h1+l1| */
  s_mul (&h2, &l2, T[i][2], svh, svl); /* T[i][2]*(svh+svl) */
  /* |T[i][2] - exp(T[i][0])| < 2^-17 ulp(T[i][2]) <= 2^-69 |T[i][2]|
     and |svh + svl - sinh(v)| < 2^-63.23*|svh + svl| thus
     |h2+l2-exp(T[i][0])*sinh(v)| < 2^-63.20*|h2+l2| */
  fast_two_sum (h, l, h1, h2);
  *l += l1 + l2;

  /* 2^-59.79 < 0x1.29p-60 and 2^-63.20 < 0x1.bep-64 */
  return 0x1.29p-60 * h1 + 0x1.bep-64 * (h2 > 0 ? h2 : -h2);
}

static void
cr_sinh_accurate (double *h, double *l, double x)
{
  static const double magic = 0x1.70f77fc88ae3cp6;
  int k = __builtin_round (magic * x);
  int i = k >> 8, j = k & 0xff;
  double v = x - T[i][0];
  double w = v - U[j][0];
  eval_S2 (h, l, w);
  if (k == 0)
  {
    static double exceptions[][3] = {
      {0x1.1bd15d167005p-11, 0x1.1bd15dff0122ap-11, 0x1.0000000000001p-64},
      {0x1.92a2ee78ed49cp-23, 0x1.92a2ee78ed4c6p-23, -0x1.0000000000001p-76},
      {0x1.bcee70ebe7ec9p-25, 0x1.bcee70ebe7ecdp-25, -0x1.fffffffffffffp-79},
      {0x1.e72460254649ap-19, 0x1.e72460254ae19p-19, 0x1.fffffffffffffp-73},
    };
    for (int i = 0; i < 4; i++)
      if (x == exceptions[i][0])
      {
        *h = exceptions[i][1];
        *l = exceptions[i][2];
        break;
      }
    return;
  }
  *h = *l = 0;
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

  //  if (x == TRACE) printf ("fast path failed\n");

  /* FIXME: update error analysis of fast path with new U[] (k=13) */

  /* Special case for small numbers, to avoid underflow issues in the accurate
     path: for |x| <= 0x1.7137449123ef6p-26, |sinh(x) - x| < ulp(x)/2,
     thus sinh(x) rounds to the same value than x + 2^-54*x. */
  if (v.f <= 0x1.7137449123ef6p-26)
    return __builtin_fma (x, 0x1p-54, x);

  cr_sinh_accurate (&h, &l, v.f);
  h *= sign[s];
  l *= sign[s];
  return h + l;
}
