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
       Hewlett-Packard Professional Books, 2000, Chapter 16.

   The argument reduction here differs from that proposed in [1]:
   - we write x = t + u + w where sinh(t), cosh(t), sinh(u), cosh(u)
     are approximated by a binary64 numbers from a table (produced
     using Gal's technique, thus with some extra bits of accuracy)
   - if t=u=0, approximate directly sinh(w) using a minimax polynomial
   - if t=0, write v=u+w, and approximate sinh(v) as sinh(u)*cosh(w)
     +cosh(u)*sinh(w) where sinh(u), cosh(u) are from the table, and cosh(w),
     sinh(w) are approximated using minimax polynomials
   - otherwise, approximate sinh(t+v) as sinh(t)*cosh(v)+cosh(t)*sinh(v)
     where sinh(v) is obtained as above, cosh(v) is obtained similarly,
     and sinh(t), cosh(t) are read from a table (a little variant is used
     here since it is not always possible to find binary64 values t such that
     both sinh(t) and cosh(t) have some extra bits of accuracy)
*/       

#include <stdio.h>
#include <stdint.h>

/* For 0 <= i < 256, T[i] = {xi, si, ei} such that xi is near i*2^8/magic
   with magic = 0x1.70f77fc88ae3cp6, and si,si+ei approximate sinh(xi),cosh(xi)
   with accuracy >= 53+16 bits:
   |si - sinh(xi)| < 2^(-16-1) ulp(si), |ci - (si+ei)| < 2^(-16-1) ulp(si+ei).
   Thus si approximates sinh(xi) with relative error < 2^-17 ulp(si)/si
   <= 2^-69 (same for ci).
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

/* For 0 <= j < 256, U[j] = {xj, sj, cj} such that xj is near j/magic
   with magic = 0x1.70f77fc88ae3cp6, and sj,cj approximate sinh(xj),cosh(xj)
   with accuracy >= 53+15 bits:
   |sj - sinh(xj)| < 2^(-15-1) ulp(sj), |cj - cosh(xj)| < 2^(-15-1) ulp(cj).
   Thus sj approximates sinh(xj) with relative error < 2^-16 ulp(sj)/sj
   <= 2^-68 (same for cj).
   We have |xj - j/magic| < 9.14e-07.
   Generated with ./buildu 15 with accompanying file buildu.c. */
static const double U[256][3] = {
   {0x0p+0, 0x0p+0, 0x1p+0}, /* i=0 */
   {0x1.633d9a6f0b004p-7, 0x1.633f62784d28p-7, 0x1.0003d9ea4b182p+0}, /* i=1 */
   {0x1.633d9acbdb14p-6, 0x1.6344baf924d78p-6, 0x1.000f67c6de66fp+0}, /* i=2 */
   {0x1.0a6e368e80b08p-5, 0x1.0a7a3cf26a0e8p-5, 0x1.0022a9ef568cdp+0}, /* i=3 */
   {0x1.633d9b02cb4bdp-5, 0x1.635a1c3bb9e99p-5, 0x1.003da0f6336bcp+0}, /* i=4 */
   {0x1.bc0d006159f75p-5, 0x1.bc44ad84ffa61p-5, 0x1.00604dac5dd3ap+0}, /* i=5 */
   {0x1.0a6e361bd3d1bp-4, 0x1.0a9e519f6bb3dp-4, 0x1.008ab1200b867p+0}, /* i=6 */
   {0x1.36d5e3bb9d063p-4, 0x1.372249ca23bf5p-4, 0x1.00bccc8b68511p+0}, /* i=7 */
   {0x1.633da0e4c393fp-4, 0x1.63afae0a28e4dp-4, 0x1.00f6a18d2f261p+0}, /* i=8 */
   {0x1.8fa5496393047p-4, 0x1.9047b1345aaeep-4, 0x1.013831b255857p+0}, /* i=9 */
   {0x1.bc0d00f20d1a3p-4, 0x1.bcebcea367b7ep-4, 0x1.01817f2724798p+0}, /* i=10 */
   {0x1.e874b52990dddp-4, 0x1.e99d4bb4f8ebep-4, 0x1.01d28c04e86eep+0}, /* i=11 */
   {0x1.0a6e32cbb27bep-3, 0x1.0b2ec018f0a36p-3, 0x1.022b5ab93b969p+0}, /* i=12 */
   {0x1.20a20c4a280bap-3, 0x1.2196e597f8691p-3, 0x1.028bedfc862dcp+0}, /* i=13 */
   {0x1.36d5e58589abap-3, 0x1.3807c1604adb7p-3, 0x1.02f448b07a3ccp+0}, /* i=14 */
   {0x1.4d09c141406dep-3, 0x1.4e820317e7014p-3, 0x1.03646e070b7e3p+0}, /* i=15 */
   {0x1.633d9624bd62dp-3, 0x1.65064e6785028p-3, 0x1.03dc612e86f8ap+0}, /* i=16 */
   {0x1.797174c5471b2p-3, 0x1.7b95619211b0dp-3, 0x1.045c261df504ap+0}, /* i=17 */
   {0x1.8fa54f6191cbap-3, 0x1.922fdc642e75ap-3, 0x1.04e3c06156a56p+0}, /* i=18 */
   {0x1.a5d928074d2f3p-3, 0x1.a8d66f0a8dd88p-3, 0x1.057334168590fp+0}, /* i=19 */
   {0x1.bc0d015c8eaf7p-3, 0x1.bf89caaf334a2p-3, 0x1.060a859ee0738p+0}, /* i=20 */
   {0x1.d240db3c3d1b9p-3, 0x1.d64a9e0705148p-3, 0x1.06a9b98797bc5p+0}, /* i=21 */
   {0x1.e874b515da9c2p-3, 0x1.ed1997bfea48fp-3, 0x1.0750d49726c2dp+0}, /* i=22 */
   {0x1.fea88e881397cp-3, 0x1.01fbb391fcb96p-2, 0x1.07ffdbd19699ep+0}, /* i=23 */
   {0x1.0a6e3565253b2p-2, 0x1.0d725fd342cccp-2, 0x1.08b6d495b3666p+0}, /* i=24 */
   {0x1.15881f937a43ep-2, 0x1.18f12365a865p-2, 0x1.0975c409d450ep+0}, /* i=25 */
   {0x1.20a20fd756ca6p-2, 0x1.24786139b6da2p-2, 0x1.0a3cb09b83241p+0}, /* i=26 */
   {0x1.2bbbfaad0523ep-2, 0x1.30086626c39e7p-2, 0x1.0b0b9f795270ep+0}, /* i=27 */
   {0x1.36d5ea564c65p-2, 0x1.3ba195eaa22abp-2, 0x1.0be2979bb5791p+0}, /* i=28 */
   {0x1.41efd2eab50dcp-2, 0x1.47443d675696fp-2, 0x1.0cc19e9167cdep+0}, /* i=29 */
   {0x1.4d09c163f00b8p-2, 0x1.52f0c3d3ff58bp-2, 0x1.0da8bc13aae09p+0}, /* i=30 */
   {0x1.5823ac33ced8p-2, 0x1.5ea779120d4b3p-2, 0x1.0e97f65225596p+0}, /* i=31 */
   {0x1.633d96515c0d2p-2, 0x1.6a68ba75b88cep-2, 0x1.0f8f54ba96989p+0}, /* i=32 */
   {0x1.6e5786b185726p-2, 0x1.7634e9ef9be77p-2, 0x1.108edf5f7f71ap+0}, /* i=33 */
   {0x1.797173ad3c024p-2, 0x1.820c5820b9929p-2, 0x1.11969d1624e9p+0}, /* i=34 */
   {0x1.848b5e2de96afp-2, 0x1.8def612e7da57p-2, 0x1.12a695dd20ccp+0}, /* i=35 */
   {0x1.8fa54e19d6296p-2, 0x1.99de6921fab7ep-2, 0x1.13bed2a81eeb5p+0}, /* i=36 */
   {0x1.9abf3a3f273cap-2, 0x1.a5d9c206ad55dp-2, 0x1.14df5aff5bf73p+0}, /* i=37 */
   {0x1.a5d9269d7c228p-2, 0x1.b1e1cc719e5d9p-2, 0x1.160837f532c1dp+0}, /* i=38 */
   {0x1.b0f31642d9859p-2, 0x1.bdf6e86469066p-2, 0x1.173972cddbf1dp+0}, /* i=39 */
   {0x1.bc0d01fb0cd7dp-2, 0x1.ca196b17304bap-2, 0x1.187313f18cffep+0}, /* i=40 */
   {0x1.c726f17a83ce1p-2, 0x1.d649ba6dbb27ep-2, 0x1.19b525a984998p+0}, /* i=41 */
   {0x1.d240d97da6a61p-2, 0x1.e28827e4d6589p-2, 0x1.1affb05b6fd51p+0}, /* i=42 */
   {0x1.dd5ac7514b785p-2, 0x1.eed520738a874p-2, 0x1.1c52bf836854cp+0}, /* i=43 */
   {0x1.e874b420cbcd6p-2, 0x1.fb30fb5073b1cp-2, 0x1.1dae5c89e91afp+0}, /* i=44 */
   {0x1.f38ea183274ecp-2, 0x1.03ce0cb8a6addp-1, 0x1.1f129215b3726p+0}, /* i=45 */
   {0x1.fea88f3f910d7p-2, 0x1.0a0b6d1fc6295p-1, 0x1.207f6ad86d868p+0}, /* i=46 */
   {0x1.04e13ea1a7edep-1, 0x1.1050cee381eeap-1, 0x1.21f4f1caa6668p+0}, /* i=47 */
   {0x1.0a6e32f58c349p-1, 0x1.169e5f1c5146cp-1, 0x1.23733165c24eep+0}, /* i=48 */
   {0x1.0ffb2a1b141bep-1, 0x1.1cf4549568f3ep-1, 0x1.24fa36a503b4p+0}, /* i=49 */
   {0x1.1588232379a08p-1, 0x1.2352df161c6bep-1, 0x1.268a0d15b6f38p+0}, /* i=50 */
   {0x1.1b151ef15d099p-1, 0x1.29ba30bcc7131p-1, 0x1.2822c10a8a623p+0}, /* i=51 */
   {0x1.20a216764fe1cp-1, 0x1.302a72b4b2256p-1, 0x1.29c45cc456b4ep+0}, /* i=52 */
   {0x1.262f058d564bdp-1, 0x1.36a3d1a79fb3cp-1, 0x1.2b6eeb7d307e6p+0}, /* i=53 */
   {0x1.2bbbf8e958f45p-1, 0x1.3d268e38d0b48p-1, 0x1.2d227dd06f692p+0}, /* i=54 */
   {0x1.3148f053f5c8ap-1, 0x1.43b2da637bfdep-1, 0x1.2edf20dca6acep+0}, /* i=55 */
   {0x1.36d5e750aae26p-1, 0x1.4a48e355b8387p-1, 0x1.30a4e0a422e0fp+0}, /* i=56 */
   {0x1.3c62da60fe4ep-1, 0x1.50e8d7992419cp-1, 0x1.3273c9a578e11p+0}, /* i=57 */
   {0x1.41efd290021cap-1, 0x1.5792f5033f9e1p-1, 0x1.344becc367264p+0}, /* i=58 */
   {0x1.477cc738cee97p-1, 0x1.5e47648c90cb8p-1, 0x1.362d5557355cep+0}, /* i=59 */
   {0x1.4d09be2649214p-1, 0x1.650660d5aee4p-1, 0x1.381813d02e948p+0}, /* i=60 */
   {0x1.5296b73d51261p-1, 0x1.6bd01dc073899p-1, 0x1.3a0c36f41ca63p+0}, /* i=61 */
   {0x1.5823ad23ec1bap-1, 0x1.72a4c907541f8p-1, 0x1.3c09cbea7a9bcp+0}, /* i=62 */
   {0x1.5db0a27c26dffp-1, 0x1.79849a79c296ap-1, 0x1.3e10e2f07872dp+0}, /* i=63 */
   {0x1.633d9f60b9246p-1, 0x1.806fd129acebfp-1, 0x1.40218eac5ecc7p+0}, /* i=64 */
   {0x1.68ca90c6740c3p-1, 0x1.87668a98839e9p-1, 0x1.423bd7e621fcdp+0}, /* i=65 */
   {0x1.6e57881e88846p-1, 0x1.8e69123a88cap-1, 0x1.445fd55d0379fp+0}, /* i=66 */
   {0x1.73e47bcf5c44ap-1, 0x1.957791f561e14p-1, 0x1.468d93e00cf8cp+0}, /* i=67 */
   {0x1.797173af5dce5p-1, 0x1.9c924a1c3cc89p-1, 0x1.48c5274c3e07dp+0}, /* i=68 */
   {0x1.7efe6928f0a28p-1, 0x1.a3b969000c7cp-1, 0x1.4b069e18bdc27p+0}, /* i=69 */
   {0x1.848b60250ca55p-1, 0x1.aaed2abfeee4p-1, 0x1.4d520b39de1bdp+0}, /* i=70 */
   {0x1.8a1853870b33p-1, 0x1.b22dc02bf5ebdp-1, 0x1.4fa77e398ca22p+0}, /* i=71 */
   {0x1.8fa5561b3d55p-1, 0x1.b97b79d937fbcp-1, 0x1.52071118b46a6p+0}, /* i=72 */
   {0x1.95324736e0b3p-1, 0x1.c0d665052cfa6p-1, 0x1.5470c8212789bp+0}, /* i=73 */
   {0x1.9abf3e640515dp-1, 0x1.c83ed95f84a98p-1, 0x1.56e4c0053af4p+0}, /* i=74 */
   {0x1.a04c2f886cf38p-1, 0x1.cfb4ffd04a95ep-1, 0x1.59630650917aap+0}, /* i=75 */
   {0x1.a5d927b80c96p-1, 0x1.d7392368e187ep-1, 0x1.5bebb41a334cp+0}, /* i=76 */
   {0x1.ab661e4dcd2fbp-1, 0x1.decb726d989fep-1, 0x1.5e7ed903870e3p+0}, /* i=77 */
   {0x1.b0f316638c665p-1, 0x1.e66c2b6ad87acp-1, 0x1.611c8a51997bfp+0}, /* i=78 */
   {0x1.b6800bbb04412p-1, 0x1.ee1b8347b1afep-1, 0x1.63c4da272bde9p+0}, /* i=79 */
   {0x1.bc0cffc5328abp-1, 0x1.f5d9b72441325p-1, 0x1.6677dda150a68p+0}, /* i=80 */
   {0x1.c199f6a9fd6efp-1, 0x1.fda70876c2f5p-1, 0x1.6935ab946ac92p+0}, /* i=81 */
   {0x1.c726ec50760fap-1, 0x1.02c1d6cdc89d6p+0, 0x1.6bfe571cb33b5p+0}, /* i=82 */
   {0x1.ccb3e72d566fdp-1, 0x1.06b7f62b5c876p+0, 0x1.6ed1f8f42c865p+0}, /* i=83 */
   {0x1.d240dabda1318p-1, 0x1.0ab5f7e4ed4a5p+0, 0x1.71b0a07b32f34p+0}, /* i=84 */
   {0x1.d7cdd152c608bp-1, 0x1.0ebc021d2f2a3p+0, 0x1.749a69155b5e6p+0}, /* i=85 */
   {0x1.dd5ac6dae4115p-1, 0x1.12ca30e2f16aap+0, 0x1.778f671611c83p+0}, /* i=86 */
   {0x1.e2e7c11db3b88p-1, 0x1.16e0a7b423818p+0, 0x1.7a8fb4639900fp+0}, /* i=87 */
   {0x1.e874b7399c00ap-1, 0x1.1aff7f8402f59p+0, 0x1.7d9b635032302p+0}, /* i=88 */
   {0x1.ee01ac297f1b9p-1, 0x1.1f26da39ab817p+0, 0x1.80b28ce387ea4p+0}, /* i=89 */
   {0x1.f38ea31f26d15p-1, 0x1.2356da398660dp+0, 0x1.83d54ab57d24ap+0}, /* i=90 */
   {0x1.f91b98fe077fbp-1, 0x1.278f9d6dd52f8p+0, 0x1.8703b32e0f2c7p+0}, /* i=91 */
   {0x1.fea88dd74d533p-1, 0x1.2bd1446320b03p+0, 0x1.8a3dded154b91p+0}, /* i=92 */
   {0x1.021ac24954c52p+0, 0x1.301bf21ecbcefp+0, 0x1.8d83e82dff514p+0}, /* i=93 */
   {0x1.04e140e0621bdp+0, 0x1.346fcb4b8d34fp+0, 0x1.90d5eb4727b61p+0}, /* i=94 */
   {0x1.07a7b30a468p+0, 0x1.38ccd8a900deap+0, 0x1.9433eeba9b341p+0}, /* i=95 */
   {0x1.0a6e358285578p+0, 0x1.3d336914add3p+0, 0x1.979e2f6e9904cp+0}, /* i=96 */
   {0x1.0d34a87abd7e2p+0, 0x1.41a36c08beff3p+0, 0x1.9b14a092bfdd6p+0}, /* i=97 */
   {0x1.0ffb27ce2d7aap+0, 0x1.461d302e83927p+0, 0x1.9e977f803bf53p+0}, /* i=98 */
   {0x1.12c1a8921b084p+0, 0x1.4aa0c6acd77e2p+0, 0x1.a226d9cf5cc34p+0}, /* i=99 */
   {0x1.158821891b2dbp+0, 0x1.4f2e431f7c5f9p+0, 0x1.a5c2bee6acd49p+0}, /* i=100 */
   {0x1.184e9d1d8762fp+0, 0x1.53c5d991cd7a1p+0, 0x1.a96b57fa98ee8p+0}, /* i=101 */
   {0x1.1b15124ec5c3p+0, 0x1.58679e684b039p+0, 0x1.ad20b540d7b87p+0}, /* i=102 */
   {0x1.1ddb927f0cadap+0, 0x1.5d13d26dc6a1ap+0, 0x1.b0e30aa4d7b91p+0}, /* i=103 */
   {0x1.20a2110d58dd2p+0, 0x1.61ca849173e17p+0, 0x1.b4b264417a8dep+0}, /* i=104 */
   {0x1.2368827501d3ep+0, 0x1.668bc53d197bdp+0, 0x1.b88ecf36bb19p+0}, /* i=105 */
   {0x1.262f062c6e8f5p+0, 0x1.6b57ef1760dccp+0, 0x1.bc78952d1d052p+0}, /* i=106 */
   {0x1.28f57e79e497dp+0, 0x1.702ef3da08fap+0, 0x1.c06faa9094fa5p+0}, /* i=107 */
   {0x1.2bbbfae83fc8ep+0, 0x1.751113bcb54fcp+0, 0x1.c47443ec26c4fp+0}, /* i=108 */
   {0x1.2e827522f72a5p+0, 0x1.79fe69446d989p+0, 0x1.c8867716604b1p+0}, /* i=109 */
   {0x1.3148f084a4656p+0, 0x1.7ef720528e15ep+0, 0x1.cca6684e8f50dp+0}, /* i=110 */
   {0x1.340f691eb2ceep+0, 0x1.83fb581cb81cap+0, 0x1.d0d43175bb3d8p+0}, /* i=111 */
   {0x1.36d5e81bfbe23p+0, 0x1.890b47ef80635p+0, 0x1.d51000a4fcd2ep+0}, /* i=112 */
   {0x1.399c5cf37de74p+0, 0x1.8e26f88767985p+0, 0x1.d959dd2229e57p+0}, /* i=113 */
   {0x1.3c62dd7d54b1p+0, 0x1.934eb97c7422dp+0, 0x1.ddb209c34dcbfp+0}, /* i=114 */
   {0x1.3f29598287d2bp+0, 0x1.98829493fdf03p+0, 0x1.e2188ed469cfep+0}, /* i=115 */
   {0x1.41efd294a1f5bp+0, 0x1.9dc2b4a033c1ap+0, 0x1.e68d907d04bc4p+0}, /* i=116 */
   {0x1.44b6506d73441p+0, 0x1.a30f50bc60988p+0, 0x1.eb113d8898745p+0}, /* i=117 */
   {0x1.477ccbe82e971p+0, 0x1.a868842604797p+0, 0x1.efa3ad3104bdcp+0}, /* i=118 */
   {0x1.4a4347ae24cdbp+0, 0x1.adce7d224f8c2p+0, 0x1.f44506ff72f82p+0}, /* i=119 */
   {0x1.4d09bdcd00bf2p+0, 0x1.b341598f346c7p+0, 0x1.f8f56486876d1p+0}, /* i=120 */
   {0x1.4fd04098b9336p+0, 0x1.b8c167925863dp+0, 0x1.fdb5091542e57p+0}, /* i=121 */
   {0x1.5296b7d08e393p+0, 0x1.be4ea18030985p+0, 0x1.0141f7e461c68p+1}, /* i=122 */
   {0x1.555d353cdb7ap+0, 0x1.c3e9558836927p+0, 0x1.03b12e2677b12p+1}, /* i=123 */
   {0x1.5823a916e9e1ep+0, 0x1.c9918ee05e92p+0, 0x1.06282c2c14456p+1}, /* i=124 */
   {0x1.5aea29f68cff1p+0, 0x1.cf47a743b9067p+0, 0x1.08a719114008bp+1}, /* i=125 */
   {0x1.5db0a487c2d9cp+0, 0x1.d50ba31b8900bp+0, 0x1.0b2df6ca40be9p+1}, /* i=126 */
   {0x1.607723b0cc5d9p+0, 0x1.daddc568cc0e7p+0, 0x1.0dbce2b2db7b7p+1}, /* i=127 */
   {0x1.633d9e439895bp+0, 0x1.e0be27c04def9p+0, 0x1.1053e8094c879p+1}, /* i=128 */
   {0x1.6604186659afep+0, 0x1.e6ad0007cd65ep+0, 0x1.12f31e89a8eeep+1}, /* i=129 */
   {0x1.68ca92bb07639p+0, 0x1.ecaa7d4a62a7fp+0, 0x1.159a9afdd3eb5p+1}, /* i=130 */
   {0x1.6b910b91e7a06p+0, 0x1.f2b6c9fd4e6ecp+0, 0x1.184a7034b286fp+1}, /* i=131 */
   {0x1.6e5789b818b5ap+0, 0x1.f8d223a7941ap+0, 0x1.1b02b9874fb0ep+1}, /* i=132 */
   {0x1.711e02af0b799p+0, 0x1.fefca23c1235bp+0, 0x1.1dc381a05c484p+1}, /* i=133 */
   {0x1.73e47d49fcc75p+0, 0x1.029b4222f757p+1, 0x1.208ce46956a66p+1}, /* i=134 */
   {0x1.76aaf67d22c7ap+0, 0x1.05bff973236a6p+1, 0x1.235ef44b0e72ep+1}, /* i=135 */
   {0x1.79717230ef3adp+0, 0x1.08ec93b6a5d1ap+1, 0x1.2639cafd1cdfdp+1}, /* i=136 */
   {0x1.7c37f15f629d9p+0, 0x1.0c212a98283d3p+1, 0x1.291d7f9a7985p+1}, /* i=137 */
   {0x1.7efe68030af6cp+0, 0x1.0f5dc8cfcdf8cp+1, 0x1.2c0a1bc83bb3bp+1}, /* i=138 */
   {0x1.81c4e6c6db2c7p+0, 0x1.12a29abb0b75bp+1, 0x1.2effc79ad798cp+1}, /* i=139 */
   {0x1.848b633655629p+0, 0x1.15efad50ab23ap+1, 0x1.31fe8ed4504a6p+1}, /* i=140 */
   {0x1.8751de1804e05p+0, 0x1.19451ad8a3dcfp+1, 0x1.35068949e1721p+1}, /* i=141 */
   {0x1.8a185944089c7p+0, 0x1.1ca2ff31decf9p+1, 0x1.3817d051dbc02p+1}, /* i=142 */
   {0x1.8cded2771f5abp+0, 0x1.2009718484b8fp+1, 0x1.3b327904009d2p+1}, /* i=143 */
   {0x1.8fa54d6596351p+0, 0x1.23789089df1c8p+1, 0x1.3e569f6c95ad6p+1}, /* i=144 */
   {0x1.926bc73ee7fcfp+0, 0x1.26f0733cb579p+1, 0x1.418458913fd67p+1}, /* i=145 */
   {0x1.953241613d998p+0, 0x1.2a713605af0b1p+1, 0x1.44bbbe7a8602bp+1}, /* i=146 */
   {0x1.97f8c030d80fap+0, 0x1.2dfaf98321b5ap+1, 0x1.47fcef1fb45f8p+1}, /* i=147 */
   {0x1.9abf411d90f7fp+0, 0x1.318dd5d280dap+1, 0x1.4b4800b524f91p+1}, /* i=148 */
   {0x1.9d85b616d38ap+0, 0x1.3529d42d8bap+1, 0x1.4e9cfbb34dfddp+1}, /* i=149 */
   {0x1.a04c32883458ap+0, 0x1.38cf298d11838p+1, 0x1.51fc10fd32e11p+1}, /* i=150 */
   {0x1.a312ad4e38e56p+0, 0x1.3c7de621cf1cp+1, 0x1.55654f90310ebp+1}, /* i=151 */
   {0x1.a5d927b0d2725p+0, 0x1.403627f1cccdfp+1, 0x1.58d8d33a23ab1p+1}, /* i=152 */
   {0x1.a89fa19876f8dp+0, 0x1.43f80b7fa6215p+1, 0x1.5c56b66e17f05p+1}, /* i=153 */
   {0x1.ab661fd1c5eacp+0, 0x1.47c3b44fdeeecp+1, 0x1.5fdf1a30b807dp+1}, /* i=154 */
   {0x1.ae2c98bdad4c6p+0, 0x1.4b993270c989bp+1, 0x1.63720d73060bep+1}, /* i=155 */
   {0x1.b0f315718efd1p+0, 0x1.4f78afeece7fdp+1, 0x1.670fb766792bp+1}, /* i=156 */
   {0x1.b3b98ed8efe2p+0, 0x1.536240bfef161p+1, 0x1.6ab82aae2c36cp+1}, /* i=157 */
   {0x1.b6800f7b22ae1p+0, 0x1.575611f43a5e1p+1, 0x1.6e6b916f1514ap+1}, /* i=158 */
   {0x1.b946843b07109p+0, 0x1.5b5426aaad79dp+1, 0x1.7229ee91ec4d7p+1}, /* i=159 */
   {0x1.bc0d00ca93839p+0, 0x1.5f5cb9e2d0f4bp+1, 0x1.75f37965b80bbp+1}, /* i=160 */
   {0x1.bed37ebe3464bp+0, 0x1.636fe18d705e7p+1, 0x1.79c8468afb649p+1}, /* i=161 */
   {0x1.c199fa1497cbfp+0, 0x1.678db72509cap+1, 0x1.7da86df576378p+1}, /* i=162 */
   {0x1.c46070cf0661fp+0, 0x1.6bb6573f2050cp+1, 0x1.81940a894fc56p+1}, /* i=163 */
   {0x1.c726ec8d02f02p+0, 0x1.6fe9f051f103cp+1, 0x1.858b4813642d4p+1}, /* i=164 */
   {0x1.c9ed671d2fc9bp+0, 0x1.74289970b6dbap+1, 0x1.898e3c5e14d99p+1}, /* i=165 */
   {0x1.ccb3e56995d8p+0, 0x1.78727ae11c5f9p+1, 0x1.8d9d0d7a5334cp+1}, /* i=166 */
   {0x1.cf7a5c2680e96p+0, 0x1.7cc7a41f95149p+1, 0x1.91b7ca0cb79f3p+1}, /* i=167 */
   {0x1.d240dc37e3352p+0, 0x1.812850f5304dp+1, 0x1.95deaab94fb93p+1}, /* i=168 */
   {0x1.d50757371bf46p+0, 0x1.85948c7d3308p+1, 0x1.9a11ba0b8c8d7p+1}, /* i=169 */
   {0x1.d7cdd1f3a32d6p+0, 0x1.8a0c8051b99d4p+1, 0x1.9e511f804693ep+1}, /* i=170 */
   {0x1.da944cf36e9a8p+0, 0x1.8e904fb5b41f1p+1, 0x1.a29cfc9c2546fp+1}, /* i=171 */
   {0x1.dd5ac42a5199fp+0, 0x1.932016c322f47p+1, 0x1.a6f56c1aaf5bbp+1}, /* i=172 */
   {0x1.e02142cd43d3cp+0, 0x1.97bc0b205d2bcp+1, 0x1.ab5aa119a7868p+1}, /* i=173 */
   {0x1.e2e7bc4344e63p+0, 0x1.9c643b707d979p+1, 0x1.afcca98fb555dp+1}, /* i=174 */
   {0x1.e5ae3b7f31156p+0, 0x1.a118ddf5f65eap+1, 0x1.b44bb948af1eap+1}, /* i=175 */
   {0x1.e874b50c1d9fep+0, 0x1.a5da03831ba06p+1, 0x1.b8d7e0577368bp+1}, /* i=176 */
   {0x1.eb3b3040e3a91p+0, 0x1.aaa7dd2e74e1cp+1, 0x1.bd714dacf2e67p+1}, /* i=177 */
   {0x1.ee01a9ce1e7a3p+0, 0x1.af828a3ddbdfbp+1, 0x1.c2181f39c64ffp+1}, /* i=178 */
   {0x1.f0c822d68c8efp+0, 0x1.b46a32023dac2p+1, 0x1.c6cc7aa98c0a7p+1}, /* i=179 */
   {0x1.f38ea58294d07p+0, 0x1.b95f0c77c8a94p+1, 0x1.cb8e95b44ba76p+1}, /* i=180 */
   {0x1.f6551ce3d4d4bp+0, 0x1.be611a5e62d83p+1, 0x1.d05e711147bb3p+1}, /* i=181 */
   {0x1.f91b95fd9b4fcp+0, 0x1.c370997fb61e5p+1, 0x1.d53c4819523acp+1}, /* i=182 */
   {0x1.fbe2106d41d3dp+0, 0x1.c88db038730e2p+1, 0x1.da283fb26c7e4p+1}, /* i=183 */
   {0x1.fea8943d2ca97p+0, 0x1.cdb8950bf8154p+1, 0x1.df228c5a5a75ep+1}, /* i=184 */
   {0x1.00b785dd10f72p+1, 0x1.d2f1474a47642p+1, 0x1.e42b2d65ef6c6p+1}, /* i=185 */
   {0x1.021ac25d1ef66p+1, 0x1.d83808e24427ap+1, 0x1.e942625f1dc82p+1}, /* i=186 */
   {0x1.037e06f3d02e4p+1, 0x1.dd8d1ee13097ep+1, 0x1.ee686dee908ap+1}, /* i=187 */
   {0x1.04e13df1efa17p+1, 0x1.e2f05eb161312p+1, 0x1.f39d26f050449p+1}, /* i=188 */
   {0x1.06447913ca60cp+1, 0x1.e862361ccca9ap+1, 0x1.f8e0f77014f12p+1}, /* i=189 */
   {0x1.07a7b9b08551dp+1, 0x1.ede2d4ed68d18p+1, 0x1.fe340dae0cb02p+1}, /* i=190 */
   {0x1.090af79bf7b44p+1, 0x1.f372454991c9bp+1, 0x1.01cb39bfc1807p+2}, /* i=191 */
   {0x1.0a6e30c5ef91ep+1, 0x1.f910a9484c9d3p+1, 0x1.048424efcd623p+2}, /* i=192 */
   {0x1.0bd16cb65a15p+1, 0x1.febe4a8efa729p+1, 0x1.0744ec14d9a1fp+2}, /* i=193 */
   {0x1.0d34aeae2637bp+1, 0x1.023db15dfe647p+2, 0x1.0a0dab22f3b8fp+2}, /* i=194 */
   {0x1.0e97ecb6dcf1cp+1, 0x1.0523fa79655b5p+2, 0x1.0cde63a0f8d3p+2}, /* i=195 */
   {0x1.0ffb2589abef8p+1, 0x1.081213f6f414p+2, 0x1.0fb7285798238p+2}, /* i=196 */
   {0x1.115e675cf4052p+1, 0x1.0b08328bf6b61p+2, 0x1.12982c7e440dep+2}, /* i=197 */
   {0x1.12c1a6858042fp+1, 0x1.0e0654578d6dcp+2, 0x1.15816e449674fp+2}, /* i=198 */
   {0x1.1424e19bccc4ap+1, 0x1.110c8d24bd03ap+2, 0x1.187300e843871p+2}, /* i=199 */
   {0x1.1588237f32c99p+1, 0x1.141b0c09c48a4p+2, 0x1.1b6d123f592dfp+2}, /* i=200 */
   {0x1.16eb624274e9cp+1, 0x1.1731d2d516a64p+2, 0x1.1e6fa40df411dp+2}, /* i=201 */
   {0x1.184e9efa40a35p+1, 0x1.1a50fb8b0e3a8p+2, 0x1.217acfabd276ap+2}, /* i=202 */
   {0x1.19b1d3c9b72c7p+1, 0x1.1d7890ab4e11ap+2, 0x1.248e9f50c3427p+2}, /* i=203 */
   {0x1.1b15162cec7a9p+1, 0x1.20a8db8fb38bfp+2, 0x1.27ab5a8d82f52p+2}, /* i=204 */
   {0x1.1c785207aed8dp+1, 0x1.23e1c6d36ebf5p+2, 0x1.2ad0ec83b0afcp+2}, /* i=205 */
   {0x1.1ddb909a1ae42p+1, 0x1.272380a0092a2p+2, 0x1.2dff8245f0a3fp+2}, /* i=206 */
   {0x1.1f3ed10edf94bp+1, 0x1.2a6e20429e98ep+2, 0x1.3137329a766cp+2}, /* i=207 */
   {0x1.20a20e4c262b1p+1, 0x1.2dc1b2ef48a19p+2, 0x1.34780a68d1713p+2}, /* i=208 */
   {0x1.22054d61213fdp+1, 0x1.311e5e58fdd55p+2, 0x1.37c22e8f1acbbp+2}, /* i=209 */
   {0x1.236887f3161fep+1, 0x1.34842ce53f769p+2, 0x1.3b15a938624a4p+2}, /* i=210 */
   {0x1.24cbc5d6d5733p+1, 0x1.37f34be0f1061p+2, 0x1.3e72a6bcae166p+2}, /* i=211 */
   {0x1.262f032649238p+1, 0x1.3b6bcc34deb6p+2, 0x1.41d937acaaccap+2}, /* i=212 */
   {0x1.27923f2ffa89cp+1, 0x1.3eedc6cdd40e3p+2, 0x1.45497473918p+2}, /* i=213 */
   {0x1.28f5809d451a2p+1, 0x1.427967b2205b4p+2, 0x1.48c3883ab8219p+2}, /* i=214 */
   {0x1.2a58bd73dde5dp+1, 0x1.460eb0bafb787p+2, 0x1.4c4774d28368bp+2}, /* i=215 */
   {0x1.2bbbfb166f167p+1, 0x1.49adcb420bd0dp+2, 0x1.4fd562cab8c29p+2}, /* i=216 */
   {0x1.2d1f382ea6d2ep+1, 0x1.4d56cfb18a874p+2, 0x1.536d6a19e35d1p+2}, /* i=217 */
   {0x1.2e82700750344p+1, 0x1.5109cd90ebf2dp+2, 0x1.570f99fd1aef7p+2}, /* i=218 */
   {0x1.2fe5b50590476p+1, 0x1.54c712bf7491ap+2, 0x1.5abc3ef4f7a22p+2}, /* i=219 */
   {0x1.3148efa73a1d9p+1, 0x1.588e7cd065e99p+2, 0x1.5e73372edd994p+2}, /* i=220 */
   {0x1.32ac2dde293e1p+1, 0x1.5c604e8e68ef6p+2, 0x1.6234c45037596p+2}, /* i=221 */
   {0x1.340f6bac4f3f5p+1, 0x1.603c9a8b5ded8p+2, 0x1.6600f89d013f2p+2}, /* i=222 */
   {0x1.3572a48cf1aacp+1, 0x1.642371b4566e1p+2, 0x1.69d7e4b9bf698p+2}, /* i=223 */
   {0x1.36d5e5aef27e4p+1, 0x1.6815174f58022p+2, 0x1.6db9cadb3cafdp+2}, /* i=224 */
   {0x1.38391e0b09ab2p+1, 0x1.6c117950e6b87p+2, 0x1.71a6993dd9fb8p+2}, /* i=225 */
   {0x1.399c621b78e6dp+1, 0x1.7018f162936fap+2, 0x1.759ea82dd254fp+2}, /* i=226 */
   {0x1.3aff9f47f6921p+1, 0x1.742b68c9057d8p+2, 0x1.79a1e14878854p+2}, /* i=227 */
   {0x1.3c62e55c36e4bp+1, 0x1.78492d55f07a2p+2, 0x1.7db0913ce9968p+2}, /* i=228 */
   {0x1.3dc619041bb9p+1, 0x1.7c720d39f7679p+2, 0x1.81ca86eeed95ep+2}, /* i=229 */
   {0x1.3f295876fc34cp+1, 0x1.80a682c4ee3ffp+2, 0x1.85f03af57c13fp+2}, /* i=230 */
   {0x1.408c989b0304bp+1, 0x1.84e68d5192b0cp+2, 0x1.8a21acb1e0735p+2}, /* i=231 */
   {0x1.41efd921048a3p+1, 0x1.89324cb7f7e1ep+2, 0x1.8e5efb8f2d893p+2}, /* i=232 */
   {0x1.43531022de11cp+1, 0x1.8d89c2f561ed8p+2, 0x1.92a82980e58cep+2}, /* i=233 */
   {0x1.44b65018c33d4p+1, 0x1.91ed4b3eaec19p+2, 0x1.96fd908fad69p+2}, /* i=234 */
   {0x1.46198afe18c5dp+1, 0x1.965cdb3ddf20ap+2, 0x1.9b5f2688144bbp+2}, /* i=235 */
   {0x1.477cca027c6eep+1, 0x1.9ad8b26f03da7p+2, 0x1.9fcd2a1dd4007p+2}, /* i=236 */
   {0x1.48e0042186657p+1, 0x1.9f60d62d9fafdp+2, 0x1.a447a09baa115p+2}, /* i=237 */
   {0x1.4a43451bdb06fp+1, 0x1.a3f58fde7effp+2, 0x1.a8ced28983037p+2}, /* i=238 */
   {0x1.4ba686151e7a1p+1, 0x1.a896ec819f8e8p+2, 0x1.ad62ccc32026dp+2}, /* i=239 */
   {0x1.4d09bc035c261p+1, 0x1.ad44ea53b484ap+2, 0x1.b2038d86db01fp+2}, /* i=240 */
   {0x1.4e6cfe6d26714p+1, 0x1.b1fffccd1b17fp+2, 0x1.b6b18706d9c75p+2}, /* i=241 */
   {0x1.4fd0379266fa7p+1, 0x1.b6c7fe6444836p+2, 0x1.bb6c942275fb5p+2}, /* i=242 */
   {0x1.513379a1812abp+1, 0x1.bb9d528ef212p+2, 0x1.c0351741d81b1p+2}, /* i=243 */
   {0x1.5296ba9ce6fb2p+1, 0x1.c07ffc2c28a3ep+2, 0x1.c50b133ebf77dp+2}, /* i=244 */
   {0x1.53f9f922f8c36p+1, 0x1.c5701bcd6e995p+2, 0x1.c9eea85521f2fp+2}, /* i=245 */
   {0x1.555d326c961cbp+1, 0x1.ca6dcd28210bcp+2, 0x1.cedff1f1bd3b8p+2}, /* i=246 */
   {0x1.56c06f34907ecp+1, 0x1.cf7955fa39a6ap+2, 0x1.d3df3524b1559p+2}, /* i=247 */
   {0x1.5823afb80f703p+1, 0x1.d492de6a3c95dp+2, 0x1.d8ec99b4e9c8dp+2}, /* i=248 */
   {0x1.5986e95a4850ap+1, 0x1.d9ba668f94169p+2, 0x1.de081fb8b4e1fp+2}, /* i=249 */
   {0x1.5aea27b72ad74p+1, 0x1.def0411d1d4d4p+2, 0x1.e3321921ab733p+2}, /* i=250 */
   {0x1.5c4d68dd59eadp+1, 0x1.e4348f5db078ep+2, 0x1.e86aa6f0d33c7p+2}, /* i=251 */
   {0x1.5db0a31f017f6p+1, 0x1.e98754e5b66c1p+2, 0x1.edb1ccb178051p+2}, /* i=252 */
   {0x1.5f13e7c503116p+1, 0x1.eee8fd3e010c3p+2, 0x1.f307f5028e2ffp+2}, /* i=253 */
   {0x1.60771eeb11416p+1, 0x1.f459550104384p+2, 0x1.f86ceceb8980dp+2}, /* i=254 */
   {0x1.61da57ecb5a17p+1, 0x1.f9d8c189bc544p+2, 0x1.fde118f05739dp+2}, /* i=255 */
};

/* The following table is used in the accurate path only.
   For each i, 0 <= i < 256, let {xi, si, ei} be the T[i] values,
   such that si and si+ei approximate sinh(xi) and cosh(xi) with a few
   extra bits, then si + Tl[i][0] and si + ei + Tl[j][1] approximate
   sinh(xi) and cosh(xi) with at least 107 bits.
   Generated with build_table_Tl(T0,T1,T2) from the file sinh.sage,
   where T0,T1,T2 are printed using the printT() routine below. */
static const double Tl[256][2] = {
   {0x0p+0, 0x0p+0}, /* i=0 */
   {0x1.06aceafcaf699p-68, 0x1.89951332da0f8p-60}, /* i=1 */
   {-0x1.8908874e7cf5fp-63, -0x1.0a4b871342323p-63}, /* i=2 */
   {-0x1.4ffcceae426f7p-60, -0x1.1ecf389f4e1cbp-68}, /* i=3 */
   {0x1.4ebc34cc22f4fp-55, -0x1.d328ffac3f7adp-71}, /* i=4 */
   {-0x1.5d23181c86ca5p-52, 0x1.98a48963b343bp-77}, /* i=5 */
   {0x1.a54ae33997becp-47, 0x1.2b27d74d5e2bp-79}, /* i=6 */
   {0x1.ed3724ec15057p-43, -0x1.723b02b47a002p-84}, /* i=7 */
   {-0x1.ef3d2b27d318bp-42, -0x1.63098ba71cd3fp-87}, /* i=8 */
   {-0x1.dafca840f6054p-38, 0x1.2a18eedcecd0ep-92}, /* i=9 */
   {0x1.5ac56fdf3db22p-33, 0x1.863ead7dc0445p-96}, /* i=10 */
   {0x1.0fde24da43facp-27, -0x1.42fb46ff1f236p-100}, /* i=11 */
   {0x1.90bd1d0f98439p-23, 0x1.1090ef61fcdp-106}, /* i=12 */
   {-0x1.b9f6238f23666p-23, 0x1.f41b9508bcp-109}, /* i=13 */
   {-0x1.f355b53ba7163p-16, 0x1.c1d0fc3ep-111}, /* i=14 */
   {-0x1.3edecec95bd27p-13, -0x1.82e5cp-119}, /* i=15 */
   {0x1.1c6067dbca196p-7, 0x1.bfc4p-119}, /* i=16 */
   {-0x1.a2585e249841ep-3, 0x1.dap-125}, /* i=17 */
   {-0x1.1ea4f401a237ap+1, 0x0p+0}, /* i=18 */
   {0x1.a5dc4fa9d2411p+2, 0x0p+0}, /* i=19 */
   {-0x1.baefb695ffbd9p+9, 0x0p+0}, /* i=20 */
   {0x1.ca3c0a142ce4p+10, -0x1p-116}, /* i=21 */
   {-0x1.9299545111483p+9, 0x1p-112}, /* i=22 */
   {0x1.64520f059e7d6p+21, 0x0p+0}, /* i=23 */
   {-0x1.04c70c6690d06p+25, 0x0p+0}, /* i=24 */
   {0x1.7b97fcc17e8d4p+28, 0x0p+0}, /* i=25 */
   {0x1.58c3230f19effp+31, 0x0p+0}, /* i=26 */
   {0x1.a939ec6223421p+37, 0x0p+0}, /* i=27 */
   {-0x1.2e45eb231611p+41, 0x0p+0}, /* i=28 */
   {-0x1.88e6636b240f3p+45, 0x0p+0}, /* i=29 */
   {0x1.b4cef1da1d5d2p+49, 0x0p+0}, /* i=30 */
   {-0x1.7e339a620344ap+51, 0x0p+0}, /* i=31 */
   {0x1.091c5336824f9p+57, 0x0p+0}, /* i=32 */
   {0x1.4df2d9154c47ap+61, 0x0p+0}, /* i=33 */
   {0x1.758ef6205a121p+62, 0x0p+0}, /* i=34 */
   {-0x1.418a1efe98babp+69, 0x0p+0}, /* i=35 */
   {0x1.b4fa83a196d3ep+71, 0x0p+0}, /* i=36 */
   {0x1.51f15b1c1d96fp+77, 0x0p+0}, /* i=37 */
   {0x1.9c030768f6835p+79, 0x0p+0}, /* i=38 */
   {0x1.ee7c914f81c64p+83, 0x0p+0}, /* i=39 */
   {0x1.77730570a3fd8p+89, 0x0p+0}, /* i=40 */
   {0x1.9b712e3a55c38p+93, 0x0p+0}, /* i=41 */
   {0x1.f6e7d065cae93p+96, 0x0p+0}, /* i=42 */
   {-0x1.8c16e27b3f96cp+100, 0x0p+0}, /* i=43 */
   {0x1.1d2d00d1352a5p+105, 0x0p+0}, /* i=44 */
   {-0x1.161c8d5cabaeap+107, 0x0p+0}, /* i=45 */
   {0x1.768f5b026ea75p+113, 0x0p+0}, /* i=46 */
   {0x1.e30c8b9543f51p+116, 0x0p+0}, /* i=47 */
   {-0x1.6e6b658417995p+118, 0x0p+0}, /* i=48 */
   {-0x1.e6dc7bdf696a5p+116, 0x0p+0}, /* i=49 */
   {0x1.38cdc947ca647p+129, 0x0p+0}, /* i=50 */
   {0x1.bd7257c83a429p+133, 0x0p+0}, /* i=51 */
   {0x1.74c3a2c6d3b7ap+137, 0x0p+0}, /* i=52 */
   {-0x1.407d1b51ac29bp+140, 0x0p+0}, /* i=53 */
   {0x1.39184547dd6a3p+143, 0x0p+0}, /* i=54 */
   {-0x1.a85a84babd68ap+149, 0x0p+0}, /* i=55 */
   {-0x1.ce237fd8c7bd9p+153, 0x0p+0}, /* i=56 */
   {-0x1.fd780c721f974p+157, 0x0p+0}, /* i=57 */
   {-0x1.8ee943fe05ae5p+159, 0x0p+0}, /* i=58 */
   {-0x1.3185f017f792cp+163, 0x0p+0}, /* i=59 */
   {-0x1.2e6fff5246abep+168, 0x0p+0}, /* i=60 */
   {0x1.84fadfe9c8801p+173, 0x0p+0}, /* i=61 */
   {-0x1.a3115b0a299a3p+177, 0x0p+0}, /* i=62 */
   {0x1.78c7132c77e36p+179, 0x0p+0}, /* i=63 */
   {0x1.6335df1d6ff7p+185, 0x0p+0}, /* i=64 */
   {-0x1.efe70270275abp+188, 0x0p+0}, /* i=65 */
   {-0x1.601d7ae443ac8p+193, 0x0p+0}, /* i=66 */
   {-0x1.848648afb8f7cp+197, 0x0p+0}, /* i=67 */
   {0x1.70a36cb771e8ap+200, 0x0p+0}, /* i=68 */
   {-0x1.5fe13bc5c698bp+205, 0x0p+0}, /* i=69 */
   {0x1.4cc2f6c46e884p+209, 0x0p+0}, /* i=70 */
   {0x1.85143897bfcb3p+213, 0x0p+0}, /* i=71 */
   {-0x1.574380474c2c3p+217, 0x0p+0}, /* i=72 */
   {0x1.3a26181c01c5p+220, 0x0p+0}, /* i=73 */
   {0x1.2135695fef858p+224, 0x0p+0}, /* i=74 */
   {0x1.a98f09bf039e5p+226, 0x0p+0}, /* i=75 */
   {-0x1.8e5f124a8b5cap+233, 0x0p+0}, /* i=76 */
   {0x1.c099cec491868p+237, 0x0p+0}, /* i=77 */
   {0x1.393238ff0cf95p+240, 0x0p+0}, /* i=78 */
   {-0x1.27f7d0129f616p+245, 0x0p+0}, /* i=79 */
   {-0x1.7e862a141287ep+249, 0x0p+0}, /* i=80 */
   {0x1.0a6aabd7ef2d3p+252, 0x0p+0}, /* i=81 */
   {0x1.9a192cefed2c8p+254, 0x0p+0}, /* i=82 */
   {0x1.04ba9e7be17a2p+260, 0x0p+0}, /* i=83 */
   {-0x1.ac9d5056e072fp+264, 0x0p+0}, /* i=84 */
   {0x1.ebed498c90296p+267, 0x0p+0}, /* i=85 */
   {-0x1.073cf09f94b4p+273, 0x0p+0}, /* i=86 */
   {-0x1.0e7a8d8dfb6f2p+277, 0x0p+0}, /* i=87 */
   {0x1.60457faaff744p+281, 0x0p+0}, /* i=88 */
   {0x1.8c7b06f501a21p+285, 0x0p+0}, /* i=89 */
   {0x1.9f4cd6bbc700dp+288, 0x0p+0}, /* i=90 */
   {0x1.94a17df9c7111p+293, 0x0p+0}, /* i=91 */
   {0x1.df796a889cad7p+297, 0x0p+0}, /* i=92 */
   {0x1.3186799828a96p+300, 0x0p+0}, /* i=93 */
   {0x1.39f0b43095221p+305, 0x0p+0}, /* i=94 */
   {0x1.be328b9869745p+309, 0x0p+0}, /* i=95 */
   {-0x1.2bbd18bec103bp+312, 0x0p+0}, /* i=96 */
   {-0x1.371bb1c77f8d9p+317, 0x0p+0}, /* i=97 */
   {0x1.9187e473e7d61p+321, 0x0p+0}, /* i=98 */
   {-0x1.660f94553432bp+325, 0x0p+0}, /* i=99 */
   {0x1.b6147c209db17p+329, 0x0p+0}, /* i=100 */
   {0x1.ecb0b92979a72p+331, 0x0p+0}, /* i=101 */
   {0x1.6f168b702c752p+337, 0x0p+0}, /* i=102 */
   {-0x1.cef356609e9ebp+339, 0x0p+0}, /* i=103 */
   {-0x1.b54ab55486df6p+341, 0x0p+0}, /* i=104 */
   {-0x1.af8f935f7d7bfp+348, 0x0p+0}, /* i=105 */
   {0x1.0c00faf62e0ep+353, 0x0p+0}, /* i=106 */
   {-0x1.e15d89f2c724dp+357, 0x0p+0}, /* i=107 */
   {-0x1.c9ace414e8aedp+360, 0x0p+0}, /* i=108 */
   {0x1.790056ad98112p+362, 0x0p+0}, /* i=109 */
   {-0x1.8ae45b93f311fp+369, 0x0p+0}, /* i=110 */
   {-0x1.37242e076a43p+372, 0x0p+0}, /* i=111 */
   {-0x1.ca66feca54132p+376, 0x0p+0}, /* i=112 */
   {0x1.7604774d52754p+381, 0x0p+0}, /* i=113 */
   {-0x1.6a7fdc4304ac1p+384, 0x0p+0}, /* i=114 */
   {0x1.b8637cc5db9b8p+388, 0x0p+0}, /* i=115 */
   {-0x1.a3235cf1f5a78p+393, 0x0p+0}, /* i=116 */
   {-0x1.0363aee3502d1p+397, 0x0p+0}, /* i=117 */
   {0x1.85a3b78f500e1p+399, 0x0p+0}, /* i=118 */
   {-0x1.406fffe30718ep+401, 0x0p+0}, /* i=119 */
   {-0x1.de6174e7029dp+408, 0x0p+0}, /* i=120 */
   {-0x1.e743484ceb813p+413, 0x0p+0}, /* i=121 */
   {-0x1.265014585e72fp+417, 0x0p+0}, /* i=122 */
   {0x1.8c853c28cb54fp+417, 0x0p+0}, /* i=123 */
   {0x1.761fcfdcef49ap+425, 0x0p+0}, /* i=124 */
   {-0x1.ad4c92d4da216p+429, 0x0p+0}, /* i=125 */
   {0x1.8862be6323cfp+430, 0x0p+0}, /* i=126 */
   {-0x1.40ad5db43b071p+437, 0x0p+0}, /* i=127 */
   {0x1.183b55c3e537p+441, 0x0p+0}, /* i=128 */
   {-0x1.7ecc50aa1b5e4p+445, 0x0p+0}, /* i=129 */
   {0x1.b52473193b2bep+449, 0x0p+0}, /* i=130 */
   {0x1.14051eda2b5a5p+452, 0x0p+0}, /* i=131 */
   {0x1.dcf7d3fe8329cp+456, 0x0p+0}, /* i=132 */
   {0x1.713d6ae750bebp+459, 0x0p+0}, /* i=133 */
   {0x1.965fb0d5960bdp+463, 0x0p+0}, /* i=134 */
   {-0x1.e91505912444bp+467, 0x0p+0}, /* i=135 */
   {0x1.88df51e9006a1p+473, 0x0p+0}, /* i=136 */
   {0x1.56422c18b1842p+474, 0x0p+0}, /* i=137 */
   {-0x1.f383324f3d451p+481, 0x0p+0}, /* i=138 */
   {-0x1.a21db13716722p+483, 0x0p+0}, /* i=139 */
   {-0x1.2c4b4ef450de7p+487, 0x0p+0}, /* i=140 */
   {-0x1.3ab80887843fcp+488, 0x0p+0}, /* i=141 */
   {0x1.cc5a7cb69910cp+497, 0x0p+0}, /* i=142 */
   {-0x1.e45df8e510fb6p+500, 0x0p+0}, /* i=143 */
   {-0x1.ee490a7742a19p+504, 0x0p+0}, /* i=144 */
   {-0x1.18b91bd8eec7ap+508, 0x0p+0}, /* i=145 */
   {0x1.8836bfd42bc1p+511, 0x0p+0}, /* i=146 */
   {0x1.4360e3d1d609dp+515, 0x0p+0}, /* i=147 */
   {-0x1.eeb5b22d6405fp+521, 0x0p+0}, /* i=148 */
   {-0x1.0fc1574af0f85p+525, 0x0p+0}, /* i=149 */
   {-0x1.1d925d1ea1f7ep+529, 0x0p+0}, /* i=150 */
   {0x1.6906f51679338p+529, 0x0p+0}, /* i=151 */
   {-0x1.323927ed9b402p+536, 0x0p+0}, /* i=152 */
   {0x1.413f32cc6add8p+541, 0x0p+0}, /* i=153 */
   {0x1.f4b84588f4b9bp+543, 0x0p+0}, /* i=154 */
   {0x1.00ca2bd8d771bp+549, 0x0p+0}, /* i=155 */
   {-0x1.75e41da870bfep+551, 0x0p+0}, /* i=156 */
   {-0x1.71abdf57f3e44p+556, 0x0p+0}, /* i=157 */
   {-0x1.08aa0167dafa3p+561, 0x0p+0}, /* i=158 */
   {0x1.5ccdd597524a7p+565, 0x0p+0}, /* i=159 */
   {-0x1.257c128b16a3bp+569, 0x0p+0}, /* i=160 */
   {0x1.99f6ec0510a89p+571, 0x0p+0}, /* i=161 */
   {-0x1.9cbc7c8bb05ddp+574, 0x0p+0}, /* i=162 */
   {-0x1.5247bbaa9ad87p+580, 0x0p+0}, /* i=163 */
   {-0x1.3041f9ad336afp+584, 0x0p+0}, /* i=164 */
   {0x1.d1d7e443e6455p+588, 0x0p+0}, /* i=165 */
   {0x1.01bf5ae72de35p+593, 0x0p+0}, /* i=166 */
   {-0x1.0221ccceecd2cp+596, 0x0p+0}, /* i=167 */
   {-0x1.2c667b5bcf6e5p+599, 0x0p+0}, /* i=168 */
   {0x1.a7b4e60920109p+603, 0x0p+0}, /* i=169 */
   {0x1.e9491f1f0f2cdp+606, 0x0p+0}, /* i=170 */
   {0x1.f90ee22052cf4p+610, 0x0p+0}, /* i=171 */
   {0x1.a6ca24fe9290dp+617, 0x0p+0}, /* i=172 */
   {-0x1.7265b2e1e0024p+620, 0x0p+0}, /* i=173 */
   {0x1.8d9bfeac5d692p+625, 0x0p+0}, /* i=174 */
   {-0x1.4a7b3e833c947p+628, 0x0p+0}, /* i=175 */
   {0x1.571055c21341dp+632, 0x0p+0}, /* i=176 */
   {-0x1.9e1fcda7c36fcp+637, 0x0p+0}, /* i=177 */
   {0x1.f9d708bb93e7cp+641, 0x0p+0}, /* i=178 */
   {-0x1.e5258660c8389p+645, 0x0p+0}, /* i=179 */
   {-0x1.e406c88e3a56bp+649, 0x0p+0}, /* i=180 */
   {-0x1.3b29f64d9fc24p+653, 0x0p+0}, /* i=181 */
   {0x1.23b4636765f1cp+654, 0x0p+0}, /* i=182 */
   {0x1.056a65c464405p+660, 0x0p+0}, /* i=183 */
   {0x1.c9c4b633bda9bp+665, 0x0p+0}, /* i=184 */
   {0x1.1cd75ebe86c9p+667, 0x0p+0}, /* i=185 */
   {-0x1.ee0a21be5723dp+672, 0x0p+0}, /* i=186 */
   {-0x1.96a2eedd5af78p+669, 0x0p+0}, /* i=187 */
   {-0x1.214900d86efep+681, 0x0p+0}, /* i=188 */
   {-0x1.35c5cbc382dbap+684, 0x0p+0}, /* i=189 */
   {-0x1.78dc634f77139p+688, 0x0p+0}, /* i=190 */
   {-0x1.638236a772c9fp+691, 0x0p+0}, /* i=191 */
   {-0x1.7421a8ece234dp+696, 0x0p+0}, /* i=192 */
   {0x1.8cc2f1b6c340ep+700, 0x0p+0}, /* i=193 */
   {-0x1.7b0b22f29fc26p+705, 0x0p+0}, /* i=194 */
   {0x1.cecc491832d13p+704, 0x0p+0}, /* i=195 */
   {0x1.f67ce173b0d05p+711, 0x0p+0}, /* i=196 */
   {0x1.1c7f01f22bee3p+714, 0x0p+0}, /* i=197 */
   {-0x1.0bbf8b14a2cf3p+718, 0x0p+0}, /* i=198 */
   {0x1.11c0bdd3625efp+721, 0x0p+0}, /* i=199 */
   {0x1.cd9a3fcbca22ap+729, 0x0p+0}, /* i=200 */
   {0x1.2a8fa60062e02p+732, 0x0p+0}, /* i=201 */
   {-0x1.4c459f92f4782p+736, 0x0p+0}, /* i=202 */
   {-0x1.a5f12240609e8p+740, 0x0p+0}, /* i=203 */
   {0x1.fac2281dcd501p+744, 0x0p+0}, /* i=204 */
   {0x1.9bd7e4cf7c3c9p+747, 0x0p+0}, /* i=205 */
   {0x1.2c16a1542ba2fp+753, 0x0p+0}, /* i=206 */
   {-0x1.328a72b57ddedp+756, 0x0p+0}, /* i=207 */
   {-0x1.8ccf6f5340f04p+760, 0x0p+0}, /* i=208 */
   {0x1.b01921d47da4dp+763, 0x0p+0}, /* i=209 */
   {-0x1.dc1bbfe6fa924p+768, 0x0p+0}, /* i=210 */
   {-0x1.44d29a0adf1b2p+773, 0x0p+0}, /* i=211 */
   {0x1.2fa8c6902b4e4p+777, 0x0p+0}, /* i=212 */
   {0x1.e0970bcbcfa3dp+781, 0x0p+0}, /* i=213 */
   {0x1.9aa21255fdf03p+784, 0x0p+0}, /* i=214 */
   {-0x1.8eb87c961f527p+788, 0x0p+0}, /* i=215 */
   {-0x1.4288c548a39d1p+792, 0x0p+0}, /* i=216 */
   {0x1.8b8b3b9842aedp+797, 0x0p+0}, /* i=217 */
   {0x1.7140767725542p+801, 0x0p+0}, /* i=218 */
   {-0x1.73aac8e435bfap+803, 0x0p+0}, /* i=219 */
   {0x1.bfdbc5451fd4ep+809, 0x0p+0}, /* i=220 */
   {-0x1.7782ca8657155p+812, 0x0p+0}, /* i=221 */
   {-0x1.904f1fbd3be86p+817, 0x0p+0}, /* i=222 */
   {-0x1.b1cc912730d28p+821, 0x0p+0}, /* i=223 */
   {0x1.ea7f86a0e9014p+825, 0x0p+0}, /* i=224 */
   {0x1.5f00241477d19p+829, 0x0p+0}, /* i=225 */
   {0x1.40d4bfddd6e0dp+833, 0x0p+0}, /* i=226 */
   {0x1.ffd47987dfc74p+837, 0x0p+0}, /* i=227 */
   {0x1.9f24294062598p+839, 0x0p+0}, /* i=228 */
   {0x1.e57e3109b9fa4p+844, 0x0p+0}, /* i=229 */
   {-0x1.ad8e9c94e3977p+849, 0x0p+0}, /* i=230 */
   {-0x1.8286a2957e77bp+853, 0x0p+0}, /* i=231 */
   {-0x1.24364d6e33b21p+856, 0x0p+0}, /* i=232 */
   {0x1.0cd24e0be255cp+861, 0x0p+0}, /* i=233 */
   {0x1.11a1d20168abep+865, 0x0p+0}, /* i=234 */
   {-0x1.37111513f699fp+869, 0x0p+0}, /* i=235 */
   {0x1.81482d1952c3ap+869, 0x0p+0}, /* i=236 */
   {0x1.d512d5dd7e3cbp+877, 0x0p+0}, /* i=237 */
   {0x1.a0d6f12d9f5a3p+881, 0x0p+0}, /* i=238 */
   {0x1.14ac9b109e7d4p+876, 0x0p+0}, /* i=239 */
   {-0x1.e65eb1e11b6b6p+889, 0x0p+0}, /* i=240 */
   {0x1.5cf8a97a38767p+887, 0x0p+0}, /* i=241 */
   {0x1.cd78e41263d53p+897, 0x0p+0}, /* i=242 */
   {0x1.5520738a9eaadp+900, 0x0p+0}, /* i=243 */
   {0x1.0270ea47ee714p+897, 0x0p+0}, /* i=244 */
   {-0x1.103da3ae27153p+909, 0x0p+0}, /* i=245 */
   {0x1.0becb01dbdae2p+911, 0x0p+0}, /* i=246 */
   {0x1.a4cfe7792b76cp+916, 0x0p+0}, /* i=247 */
   {-0x1.f18ce502a948ep+919, 0x0p+0}, /* i=248 */
   {-0x1.6f785ea892c19p+924, 0x0p+0}, /* i=249 */
   {0x1.5df0ca797ea65p+929, 0x0p+0}, /* i=250 */
   {-0x1.591b32bdcaa3cp+933, 0x0p+0}, /* i=251 */
   {0x1.c833cdff5a68cp+937, 0x0p+0}, /* i=252 */
   {0x1.03f99b046bf52p+938, 0x0p+0}, /* i=253 */
   {-0x1.d1b74f3418125p+944, 0x0p+0}, /* i=254 */
   {0x1.3300ca4eeb124p+950, 0x0p+0}, /* i=255 */
};

/* The following table is used in the accurate path only.
   For each j, 0 <= j < 256, let {xj, sj, cj} be the U[j] values,
   such that sj and cj approximate sinh(xj) and cosh(xj) with a few
   extra bits, then sj + Ul[j][0] and cj + Ul[j][1] approximate
   sinh(xj) and cosh(xj) with at least 107 bits.
   Generated with build_table_Ul(U0,U1,U2) from the file sinh.sage,
   where U0,U1,U2 are printed using the printU() routine below. */
static const double Ul[256][2] = {
   {0x0p+0, 0x0p+0}, /* i=0 */
   {-0x1.cc125d97df011p-76, -0x1.e6bae12de82cep-70}, /* i=1 */
   {0x1.30f98f760b88dp-75, -0x1.6ca2976e53d76p-71}, /* i=2 */
   {0x1.e23212f40caa6p-74, -0x1.d64a01a2ef735p-69}, /* i=3 */
   {-0x1.74d7fe5518d25p-75, -0x1.11f607ae2cde9p-74}, /* i=4 */
   {-0x1.f1e7a8971f37cp-74, 0x1.0cec70d2c9d5fp-69}, /* i=5 */
   {-0x1.ede7326aa5f09p-74, 0x1.a611a3ec0cf22p-69}, /* i=6 */
   {-0x1.525d4563052b1p-73, 0x1.93c5ed3caff24p-70}, /* i=7 */
   {-0x1.3efc8cb924729p-73, 0x1.9b2147bfdfc2p-69}, /* i=8 */
   {-0x1.041d3153ed724p-73, -0x1.034c82669fef1p-71}, /* i=9 */
   {-0x1.8e2f60ccf0ep-74, 0x1.eef11e1433e5ep-69}, /* i=10 */
   {-0x1.8f985ba67a37p-73, -0x1.b9c803125fd4p-73}, /* i=11 */
   {-0x1.740ee01ad7bc3p-72, -0x1.338cc7288b73bp-71}, /* i=12 */
   {0x1.31b3256e78c3cp-74, -0x1.eb3ae56b94f65p-69}, /* i=13 */
   {0x1.3c4997ff0424bp-72, -0x1.c4c9bf012a0ddp-69}, /* i=14 */
   {0x1.79b38681f84c3p-72, 0x1.6eb9df80ee3cep-69}, /* i=15 */
   {0x1.8fed91ae794aap-72, -0x1.1a6a8a28f04c4p-69}, /* i=16 */
   {-0x1.13d88555cac7ep-72, -0x1.0c000fa79a59p-71}, /* i=17 */
   {0x1.ffb7efc57f4cp-72, -0x1.868a414ac91ap-69}, /* i=18 */
   {0x1.338d9ce94afcap-72, -0x1.e69c8f1ed4087p-69}, /* i=19 */
   {-0x1.b92cde07ec9ffp-76, -0x1.c2dad21d7af5p-69}, /* i=20 */
   {-0x1.6740a2aebe265p-72, 0x1.305193ecab1dbp-69}, /* i=21 */
   {0x1.085e8eb198df9p-73, 0x1.8c42c2affc91ep-69}, /* i=22 */
   {0x1.28e352f80f835p-74, 0x1.6560d765ba457p-69}, /* i=23 */
   {-0x1.fc14060bbf636p-72, -0x1.268ccf3bb4962p-74}, /* i=24 */
   {0x1.f247008ac3e3ep-71, 0x1.76b57d1aa7777p-70}, /* i=25 */
   {0x1.57a6ceb3b4336p-73, -0x1.4e4c0afd5ac3fp-72}, /* i=26 */
   {0x1.d0b3fd87e7afap-71, 0x1.367e3ccd75e3cp-73}, /* i=27 */
   {0x1.d3db0418c307ap-72, -0x1.15522378229b4p-69}, /* i=28 */
   {-0x1.f7bfa9d007548p-72, 0x1.eb09489ebc265p-75}, /* i=29 */
   {0x1.95e8360490ep-71, -0x1.588a215b8d24dp-75}, /* i=30 */
   {-0x1.bfd369754dbc4p-71, -0x1.a67fccb65954bp-75}, /* i=31 */
   {-0x1.8036da78049a1p-71, 0x1.dfa388c3a803p-69}, /* i=32 */
   {-0x1.60aa699d644a7p-73, 0x1.132c1cabf36ecp-70}, /* i=33 */
   {0x1.a7b6efe96297cp-71, -0x1.0d4cad19224a7p-69}, /* i=34 */
   {-0x1.21f632c7643dep-72, -0x1.cf140022a1fe7p-69}, /* i=35 */
   {0x1.cbe0d9b08abd7p-71, -0x1.66b4ec4fdd4b1p-71}, /* i=36 */
   {-0x1.bf3acf5f270dcp-71, 0x1.685bbecc9be3bp-70}, /* i=37 */
   {-0x1.81914a780f0c2p-71, 0x1.9b4cc5b8bbd6p-69}, /* i=38 */
   {0x1.64698d6aa32fbp-72, -0x1.414d7f0b5db0bp-69}, /* i=39 */
   {0x1.b5604cbda5badp-72, -0x1.b9b20d4f6783cp-69}, /* i=40 */
   {0x1.bfa247fa3db9dp-71, 0x1.4736ed609deb7p-69}, /* i=41 */
   {0x1.e59c3dfe07ef4p-71, -0x1.a2c53292eb267p-69}, /* i=42 */
   {-0x1.bc0e71dbff684p-73, 0x1.6006480d88602p-69}, /* i=43 */
   {0x1.061fd02b0d1ecp-71, 0x1.97f10247dd219p-72}, /* i=44 */
   {0x1.52ca9d90f7528p-70, -0x1.75f26f51d468cp-69}, /* i=45 */
   {-0x1.6cafe82602e85p-70, -0x1.4f5036669f078p-75}, /* i=46 */
   {-0x1.08cbf49cc11dfp-73, -0x1.2d084a6813884p-69}, /* i=47 */
   {-0x1.65273fdbc7039p-70, -0x1.0284835c62871p-71}, /* i=48 */
   {0x1.49cd697517a38p-71, -0x1.75d1f59cf8abfp-69}, /* i=49 */
   {0x1.b484808d541e9p-70, -0x1.fadb289ea2362p-69}, /* i=50 */
   {0x1.ca3a71ee25d48p-71, -0x1.0668584094654p-70}, /* i=51 */
   {-0x1.ebece79c0cdc4p-72, 0x1.62aa793165d04p-70}, /* i=52 */
   {0x1.ed3c708fb0577p-71, 0x1.4122a761940fap-70}, /* i=53 */
   {-0x1.c56e49d6ce5e1p-70, 0x1.b938266321177p-69}, /* i=54 */
   {-0x1.c803b3cac637bp-71, -0x1.c0a480483be05p-69}, /* i=55 */
   {-0x1.c393559e63409p-70, -0x1.fd6cc2fbaedbfp-69}, /* i=56 */
   {0x1.3955e4a058217p-71, 0x1.ab4ba53b705c1p-71}, /* i=57 */
   {0x1.0d34952a2e57fp-73, 0x1.34feb237c3eafp-72}, /* i=58 */
   {-0x1.db541bc3d687ap-70, -0x1.7250e04efb13dp-69}, /* i=59 */
   {-0x1.8c73dfc89d19bp-70, 0x1.e0333989efdb6p-70}, /* i=60 */
   {-0x1.e1530f3507293p-73, -0x1.dd9736fc563d9p-70}, /* i=61 */
   {-0x1.773970fbd9fe5p-71, -0x1.5f4b88b364a8cp-69}, /* i=62 */
   {-0x1.511d2f35062acp-70, 0x1.a593df54a8efep-69}, /* i=63 */
   {0x1.a2c80dc17d90bp-70, -0x1.55563b510d2c3p-69}, /* i=64 */
   {0x1.0dce78b1dac7p-71, -0x1.9d8e94061981ap-71}, /* i=65 */
   {0x1.3d4497fd35f14p-70, 0x1.b81a8b8c3fe0ap-69}, /* i=66 */
   {-0x1.ef876397ac2cep-72, -0x1.3cf3e7dfb8492p-69}, /* i=67 */
   {0x1.e32bf03a7a561p-70, 0x1.05eb63c7453cdp-70}, /* i=68 */
   {0x1.7330935fab8e1p-71, 0x1.80981552d788ap-70}, /* i=69 */
   {-0x1.43d07261cd239p-75, -0x1.5264370bc257cp-72}, /* i=70 */
   {0x1.ee206961285fdp-71, 0x1.3cc5d4d31f75bp-69}, /* i=71 */
   {-0x1.129182ed61e77p-71, -0x1.4c1bb9ad2522dp-69}, /* i=72 */
   {-0x1.0734657db52adp-72, -0x1.c9c73971ebb56p-69}, /* i=73 */
   {-0x1.6cc4ac7ab6453p-72, -0x1.21caf584f8cfap-69}, /* i=74 */
   {-0x1.0083f9881ca35p-71, -0x1.065a9cf08ef22p-70}, /* i=75 */
   {0x1.9d3fe6fbd6c88p-70, -0x1.5b37df08c8febp-70}, /* i=76 */
   {0x1.0fd94868678b4p-70, 0x1.36ac7d063279cp-71}, /* i=77 */
   {-0x1.1b99aedc6cdd3p-70, 0x1.ef1ce2fe45f61p-69}, /* i=78 */
   {0x1.913f2ad422e09p-70, -0x1.cb40ad3af1931p-70}, /* i=79 */
   {-0x1.ed6f05ce3b18dp-70, -0x1.ef1f6fa92997ap-70}, /* i=80 */
   {0x1.d714984430cb5p-72, -0x1.e30e3cac51392p-69}, /* i=81 */
   {0x1.6ea6ed11d5593p-69, 0x1.258ccd243875dp-69}, /* i=82 */
   {0x1.3fe9878c90f24p-69, -0x1.e404beb1985c7p-71}, /* i=83 */
   {0x1.804dbffe81dffp-69, 0x1.07ffef71a553fp-69}, /* i=84 */
   {0x1.94edb6860ae1dp-69, -0x1.1f9af00117f2ep-70}, /* i=85 */
   {0x1.f88f3346c57bfp-69, -0x1.88214c48fa801p-72}, /* i=86 */
   {0x1.a978418cc69d5p-69, -0x1.4bd55b47b6ad9p-69}, /* i=87 */
   {-0x1.dcc62caeccf26p-69, 0x1.6301dbc0f4ee5p-70}, /* i=88 */
   {-0x1.113b6de9cba6ap-70, -0x1.f9c7bef394e72p-71}, /* i=89 */
   {-0x1.2fe06d0c871fap-70, 0x1.c8d03d7bccc19p-70}, /* i=90 */
   {-0x1.9a56e5024be36p-69, 0x1.a1b9fce65fda2p-69}, /* i=91 */
   {-0x1.e85f8d5f88b44p-69, 0x1.f794517b9773bp-70}, /* i=92 */
   {0x1.a46da87cd2481p-70, -0x1.f9c6cd1a01069p-69}, /* i=93 */
   {-0x1.cdc9ac7ae3979p-69, 0x1.22ebca44ef6a9p-70}, /* i=94 */
   {0x1.ad52403d4b414p-69, -0x1.37d8dc42157e8p-69}, /* i=95 */
   {0x1.cf37c155fc817p-70, 0x1.c30b25c1af414p-69}, /* i=96 */
   {0x1.e4e00199814ecp-69, -0x1.c191afc07501fp-71}, /* i=97 */
   {0x1.9ce25d7888fbbp-70, -0x1.572a889b3fb69p-73}, /* i=98 */
   {-0x1.fd1a731689ee9p-71, 0x1.c4df9e0b184aep-69}, /* i=99 */
   {-0x1.e0837b42a9db5p-69, 0x1.de3c949b6cb0fp-69}, /* i=100 */
   {0x1.86a3e9a685316p-70, -0x1.30c5583616607p-69}, /* i=101 */
   {-0x1.8a72cb1d8fd67p-75, -0x1.5cae4be6132ap-70}, /* i=102 */
   {0x1.70b78cfcdc386p-69, -0x1.6bfc485a682b7p-69}, /* i=103 */
   {-0x1.20250cf6d61efp-69, -0x1.ce24a83904d33p-69}, /* i=104 */
   {0x1.f76751ddd36eap-69, -0x1.4a9187b2e93f2p-69}, /* i=105 */
   {-0x1.09e8d0c89f73p-71, 0x1.e4e500b332b26p-72}, /* i=106 */
   {0x1.c733e1f11ab31p-70, 0x1.ff8f78984ca64p-69}, /* i=107 */
   {0x1.018ffc3922ca2p-70, 0x1.c7d4555d4515fp-70}, /* i=108 */
   {0x1.9abeb0bcf3981p-72, 0x1.bda250cfd5432p-70}, /* i=109 */
   {0x1.d36e41eaf4272p-71, 0x1.4b9700643c389p-71}, /* i=110 */
   {-0x1.a10095d7ea98cp-70, -0x1.7657a2cdd2681p-70}, /* i=111 */
   {-0x1.54d5e72a63e9cp-69, -0x1.a9bfb2a60f6afp-70}, /* i=112 */
   {-0x1.33d29f0415ca2p-69, -0x1.a8d74d1a43e6p-74}, /* i=113 */
   {0x1.c450bf26fdbp-72, 0x1.b393f363cfe96p-70}, /* i=114 */
   {-0x1.da428de1451adp-69, -0x1.25ec2dc392dp-69}, /* i=115 */
   {0x1.a2e96355e30b6p-69, -0x1.d86381e4d8e87p-69}, /* i=116 */
   {-0x1.28a852d56a4edp-69, -0x1.83277cc96eafp-69}, /* i=117 */
   {0x1.67ae00f9bae17p-70, -0x1.2932116d6533fp-69}, /* i=118 */
   {-0x1.ee059d97ec689p-69, 0x1.fadab4f330b1fp-69}, /* i=119 */
   {-0x1.ff8b8f145b91fp-69, -0x1.bd7b9743c8963p-69}, /* i=120 */
   {-0x1.b4f57161ff472p-70, 0x1.ec1476d6563c3p-69}, /* i=121 */
   {-0x1.5b0dd5b349046p-74, -0x1.c89666d7248fp-68}, /* i=122 */
   {-0x1.4d2a47d3fc853p-69, 0x1.7d228c291c71fp-68}, /* i=123 */
   {-0x1.071f40bbfcc6p-69, 0x1.9038cac5e7c17p-69}, /* i=124 */
   {0x1.f96117f7a8984p-69, -0x1.b093cfada9266p-68}, /* i=125 */
   {-0x1.d9f29cd513256p-70, -0x1.3dfec2098ef11p-68}, /* i=126 */
   {0x1.52c9375ce8e42p-70, 0x1.f56a9c61471adp-68}, /* i=127 */
   {0x1.7fed70054f934p-70, -0x1.711c554acccadp-68}, /* i=128 */
   {-0x1.6b22755c76f47p-69, -0x1.c4085d0b7c55cp-68}, /* i=129 */
   {-0x1.ca47f3f3650adp-69, 0x1.85be32ff09cd7p-71}, /* i=130 */
   {-0x1.253318bf58a74p-69, -0x1.9b22f6dbb1975p-68}, /* i=131 */
   {-0x1.fac6ae590ca9fp-69, 0x1.0720f36f66033p-68}, /* i=132 */
   {-0x1.fae47c181ae8ap-73, -0x1.8a8e55cd9ea65p-69}, /* i=133 */
   {-0x1.300abe9518c85p-72, -0x1.c425810950019p-69}, /* i=134 */
   {0x1.077ee5a85a09fp-71, 0x1.1bf8cd24a26d1p-68}, /* i=135 */
   {0x1.4937072528aa7p-69, 0x1.bcd61e769a7f7p-69}, /* i=136 */
   {0x1.c560f914dcfa8p-68, -0x1.6a5546744d556p-71}, /* i=137 */
   {0x1.cc39b5edc4e35p-68, 0x1.237e9c4df1c9ap-69}, /* i=138 */
   {-0x1.c2415a40346fcp-70, -0x1.190d55cedda2dp-70}, /* i=139 */
   {-0x1.53881dac89e96p-68, 0x1.b9ff2df504af6p-68}, /* i=140 */
   {-0x1.50890a9134accp-68, -0x1.358e8608723dcp-68}, /* i=141 */
   {-0x1.031cfc55a6ef4p-69, 0x1.08a18cbfb1a86p-70}, /* i=142 */
   {-0x1.1acd45619b29cp-68, 0x1.f264e8ef6843ep-68}, /* i=143 */
   {-0x1.105967a61e57cp-69, 0x1.1dd1bbe1c4608p-68}, /* i=144 */
   {-0x1.1a4b5366b3339p-69, -0x1.a7f4c61818459p-69}, /* i=145 */
   {-0x1.84c4f5a0d234dp-68, -0x1.a588bb67a659dp-68}, /* i=146 */
   {0x1.5ec829c41319p-70, -0x1.b757f020a361ap-68}, /* i=147 */
   {-0x1.3d3d8f0bfd99p-72, 0x1.540903e9d67aep-68}, /* i=148 */
   {-0x1.f8e48fb8b7098p-68, -0x1.a3a6adbcdf3c6p-68}, /* i=149 */
   {-0x1.f9e60e1a440b5p-68, 0x1.6972705258e14p-70}, /* i=150 */
   {-0x1.5fd0345b73233p-68, 0x1.14718d528da85p-68}, /* i=151 */
   {-0x1.f514e16df0e6p-71, 0x1.63610d3acf335p-69}, /* i=152 */
   {0x1.c373fb8185faep-70, 0x1.7a1de2017abb2p-68}, /* i=153 */
   {-0x1.2316580237529p-71, 0x1.0cf46c1122715p-68}, /* i=154 */
   {0x1.f4c87d77a1287p-70, 0x1.29e0a9a2316c4p-70}, /* i=155 */
   {-0x1.7a0b0e3cdaa05p-69, 0x1.b173c43fca317p-70}, /* i=156 */
   {0x1.c98b6f73fd8b4p-68, -0x1.93cf31fd7c2d8p-70}, /* i=157 */
   {-0x1.18c79f8e2e97ap-68, 0x1.6d21fb2bb8525p-69}, /* i=158 */
   {-0x1.b0787601298f3p-68, -0x1.f6af77ad23345p-71}, /* i=159 */
   {0x1.20f079c66d7dep-69, -0x1.38d89fe404aebp-68}, /* i=160 */
   {0x1.0b24e12370c9cp-73, 0x1.397e1ae2045f4p-68}, /* i=161 */
   {0x1.0e4aa92e73dfcp-69, -0x1.9bf5e570d55ddp-69}, /* i=162 */
   {0x1.007b0a5b67871p-72, -0x1.d0817c2c0dbafp-69}, /* i=163 */
   {-0x1.f05540451d851p-69, -0x1.5a3e4bc2025f3p-68}, /* i=164 */
   {0x1.5abe5c6ed3e4ap-68, 0x1.1486ad7d84fe5p-68}, /* i=165 */
   {0x1.7be62ab7a2101p-68, 0x1.ceef7711380e9p-69}, /* i=166 */
   {-0x1.871c36e9d4139p-68, -0x1.fb83afe442bbap-70}, /* i=167 */
   {0x1.213d9e4081dc9p-70, 0x1.25504bd141d71p-68}, /* i=168 */
   {0x1.80996ffffe784p-69, -0x1.b373de91f6a75p-69}, /* i=169 */
   {0x1.5a65c38510e4fp-68, 0x1.7603d7b8eac88p-68}, /* i=170 */
   {0x1.bed20d0bf1d87p-68, -0x1.6789c936540d3p-69}, /* i=171 */
   {-0x1.716dcb1b77895p-70, -0x1.81e8bc196b514p-69}, /* i=172 */
   {0x1.21f4e3d09c21ep-68, -0x1.005d99a56543ap-70}, /* i=173 */
   {0x1.dfbf08609e775p-68, -0x1.ed92824d605bp-68}, /* i=174 */
   {0x1.f7e63692b53fep-69, -0x1.372de32217dabp-71}, /* i=175 */
   {0x1.c0c7104756a83p-71, 0x1.bd3b4b54a33ddp-69}, /* i=176 */
   {0x1.41cbd1d7a7bf9p-68, -0x1.a6d4073d2c183p-68}, /* i=177 */
   {-0x1.385e66f085091p-69, 0x1.ce3c07dd5a267p-69}, /* i=178 */
   {-0x1.21e2f31939c01p-69, 0x1.e5195e837065p-70}, /* i=179 */
   {0x1.710ddc70f8a83p-69, -0x1.c63451b3b66d2p-68}, /* i=180 */
   {-0x1.7f14f4fc4452ep-70, -0x1.1244ea576aa7ap-68}, /* i=181 */
   {-0x1.b66c9243b00e3p-69, 0x1.21e8cac31f56bp-68}, /* i=182 */
   {-0x1.d05a86edcd3b9p-70, -0x1.17c0c973eb37p-70}, /* i=183 */
   {0x1.1324a13f64769p-70, -0x1.1d250a88c328ap-68}, /* i=184 */
   {0x1.6b9558896a52ap-68, 0x1.50c3d165671f9p-68}, /* i=185 */
   {-0x1.5461c2e8f209p-68, -0x1.ce9cda41e4b78p-69}, /* i=186 */
   {-0x1.8f32f2e702724p-68, 0x1.07612d70a0702p-70}, /* i=187 */
   {-0x1.fb25cc88bb058p-72, 0x1.916d8f16c9f57p-71}, /* i=188 */
   {0x1.3fedacb479337p-69, -0x1.572677d495ea4p-69}, /* i=189 */
   {0x1.4787b4e497b07p-68, 0x1.f9e2f5611e979p-68}, /* i=190 */
   {0x1.c68d2ad570c07p-69, -0x1.e6abae4538b3dp-67}, /* i=191 */
   {0x1.51cbc86efb9e8p-69, -0x1.8057294f6d0e6p-68}, /* i=192 */
   {0x1.c020de7181543p-68, -0x1.94ca767ba032ep-68}, /* i=193 */
   {0x1.3fe5830ff1a57p-67, 0x1.24f892180b0b2p-67}, /* i=194 */
   {0x1.336d646a0b267p-68, 0x1.c5f2047eb5bacp-68}, /* i=195 */
   {-0x1.3604f7fbad5bcp-68, 0x1.189093e059c56p-67}, /* i=196 */
   {-0x1.b6bff63be9f98p-67, 0x1.046e80d9c4e13p-69}, /* i=197 */
   {0x1.d853de738314p-68, -0x1.dfbba064c59cp-67}, /* i=198 */
   {0x1.db970f23f0bcp-67, -0x1.dd22cb754caaep-70}, /* i=199 */
   {0x1.470924218186ep-67, 0x1.8985d90a78c99p-67}, /* i=200 */
   {0x1.ac198a955c05p-72, 0x1.bd113cde03a7fp-68}, /* i=201 */
   {-0x1.152aecfd663b3p-67, -0x1.37d344e58a72cp-67}, /* i=202 */
   {0x1.3ff8910dc55a8p-67, -0x1.4ef8946b1a2f1p-67}, /* i=203 */
   {0x1.35466cf45fddfp-67, -0x1.e942b73993db2p-69}, /* i=204 */
   {-0x1.86ba464827a75p-68, -0x1.374b5238fe75bp-67}, /* i=205 */
   {-0x1.b662fede103aap-67, -0x1.0b639090694bbp-69}, /* i=206 */
   {0x1.70c6e4d5b0bdp-67, 0x1.6d71bd1d386bp-67}, /* i=207 */
   {0x1.ba9d154d0e985p-68, -0x1.ac141b0fc1ad3p-67}, /* i=208 */
   {-0x1.17007ba8c51b1p-67, -0x1.17f3ea450a94fp-67}, /* i=209 */
   {-0x1.af9cfcff75aa2p-69, 0x1.da7861b70f94cp-67}, /* i=210 */
   {0x1.535f439324327p-74, -0x1.ac01f0f30f4f7p-68}, /* i=211 */
   {-0x1.bdc60333aba51p-67, 0x1.87828ee70336fp-69}, /* i=212 */
   {-0x1.f200e9d3eb15bp-67, -0x1.c16f5d84ba9ffp-68}, /* i=213 */
   {-0x1.97af320485e52p-68, -0x1.5f78a57458157p-67}, /* i=214 */
   {0x1.47cd0fa68ba97p-67, -0x1.c3c83b4bfc515p-70}, /* i=215 */
   {0x1.d4841e5c03886p-67, -0x1.d9b4e24f0a52cp-67}, /* i=216 */
   {0x1.ed67e6a226b22p-67, -0x1.df43e8c66e09bp-67}, /* i=217 */
   {-0x1.50fde14f403c9p-68, 0x1.f480f535c75a7p-67}, /* i=218 */
   {0x1.2ac03f80ad9fbp-69, -0x1.568dace2df407p-67}, /* i=219 */
   {-0x1.c0f74f5c2cf96p-73, 0x1.f0c07bfabea77p-69}, /* i=220 */
   {0x1.aaebc950f0c18p-67, 0x1.1bf195319ec76p-67}, /* i=221 */
   {-0x1.de896af14f8b6p-67, -0x1.59a987f1bfb7cp-67}, /* i=222 */
   {-0x1.7a2a9ecd9f0afp-72, 0x1.5dda3046758b4p-68}, /* i=223 */
   {-0x1.a8325ca5c74b9p-69, 0x1.72cbbc2231025p-68}, /* i=224 */
   {-0x1.7acca2ad24d15p-69, -0x1.00f9765848bdcp-69}, /* i=225 */
   {0x1.e9bb7f553c2e3p-69, 0x1.9fa1b1cedc899p-68}, /* i=226 */
   {-0x1.c935e4fe3fa3ap-71, -0x1.c32dd48421613p-67}, /* i=227 */
   {-0x1.63a662cc7e483p-67, -0x1.f9950532d1e0ap-67}, /* i=228 */
   {0x1.da4284a2b3b04p-68, -0x1.4f95299cdeecep-68}, /* i=229 */
   {-0x1.d6d2999ac6e88p-70, 0x1.df2a03184427fp-68}, /* i=230 */
   {-0x1.e0476ffa51769p-71, 0x1.54f9ef4ce5082p-67}, /* i=231 */
   {0x1.1f030651cf7c8p-67, -0x1.344b94effe4ap-67}, /* i=232 */
   {-0x1.63c8bf00ca55bp-67, -0x1.589e97076e239p-68}, /* i=233 */
   {-0x1.748f95d9f2c24p-73, -0x1.a5d2a9a099ec3p-67}, /* i=234 */
   {-0x1.43b9f8ccd861ep-68, 0x1.d857672a61cf4p-67}, /* i=235 */
   {0x1.396c881d50c13p-67, -0x1.12a0a0ff4f684p-71}, /* i=236 */
   {-0x1.d243bf2b153f6p-70, 0x1.b9012df183604p-69}, /* i=237 */
   {-0x1.6041337ff51b3p-67, -0x1.415b42c7ee528p-67}, /* i=238 */
   {-0x1.36cd0bba43baep-71, -0x1.c64223aac9278p-67}, /* i=239 */
   {-0x1.b90d7ba8173d8p-67, 0x1.2a4c96a11fc4ap-68}, /* i=240 */
   {-0x1.3f60f0cf65b96p-68, -0x1.af54d96ae2d34p-67}, /* i=241 */
   {0x1.642fadcaed496p-69, -0x1.d22c1713be6bep-69}, /* i=242 */
   {0x1.e6fdf5e75db62p-68, 0x1.3879e04c3badp-67}, /* i=243 */
   {0x1.b8606c1c3bbe4p-71, -0x1.5d41415294b67p-68}, /* i=244 */
   {-0x1.3f57fe2e23d43p-69, -0x1.9a5c677640eabp-69}, /* i=245 */
   {0x1.82129a23715p-67, -0x1.24528efb7fc1cp-69}, /* i=246 */
   {0x1.8bbbd5809ba8ap-69, -0x1.23c5767bd1537p-68}, /* i=247 */
   {-0x1.3bb0b6b194422p-67, 0x1.6c8f89d76c4e1p-69}, /* i=248 */
   {0x1.02a35bd49318p-67, -0x1.f41a2d4dd4605p-68}, /* i=249 */
   {-0x1.41080bc3154c1p-67, -0x1.64b13bcb4353cp-67}, /* i=250 */
   {0x1.5f239966ad308p-67, 0x1.24b8d351ef35bp-70}, /* i=251 */
   {-0x1.5e4b62c9bb051p-67, 0x1.02c16542fcd2bp-67}, /* i=252 */
   {-0x1.f89aa76ab1e5p-67, -0x1.b42b78b3c565p-69}, /* i=253 */
   {-0x1.18fb74c82003p-67, 0x1.3a454097b1dfep-70}, /* i=254 */
   {0x1.78f85fb02deb6p-67, 0x1.96f01bda277cap-67}, /* i=255 */
};

typedef union { double f; uint64_t u; } d64u64;

/* Add a + b, such that *hi + *lo approximates a + b.
   Assumes |a| >= |b|.  */
static void
fast_two_sum (double *hi, double *lo, double a, double b)
{
  double e;

  //  assert (fabs (a) >= fabs (b));

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
  /* We neglect at input S2[4][0]*zl*w^7 which relatively to sinh(w) ~ w
     is less than S2[4][0]*zl*w^7/w < 2^-18.46*2^-68*2^-45.14 = 2^-131.6.
     We have |h| < 2^-12.29, we neglect in output ulp(h)*w^7, which relatively
     to w gives ulp(2^-12.29)*w^6 < 2^-110.14 */
  s_mul (h, l, *h, zh, zl);                     /* multiply by w^2 */
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

  /* since |T[i][0] - i*2^8/magic| < 2.36e-8 and
           |U[j][0] - j/magic| < 9.14e-07, we have:
     |x - T[i][0] - U[j][0]| < 0.00542055 + 2.36e-8 + 9.14e-07 < 0.00543 */

  /* we have |w| < 0.00543 */
  double swh, swl, cwh, cwl;
  eval_S (h, l, w);
  /* |h + l - sinh(w)| < 2^-67.99*|h| */

  if (k == 0)
    return __builtin_fma (0x1.02p-68, *h, 0x1p-1074);
  /* 2^-67.99 < 0x1.02p-68, and we add 2^-1074 to workaround cases
     when 0x1.02p-68 * h is rounded to zero */

  eval_C (&cwh, &cwl, w);
  /* |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| */
  
  swh = *h;
  swl = *l;
  double svh, svl, cvh, cvl, h1, l1, h2, l2;
  s_mul (&h1, &l1, U[j][1], cwh, cwl); /* U[j][1]*cosh(w) */
  /* |U[j][1] - sinh(U[j][0])| < 2^-16 ulp(U[j][1]) <= 2^-68 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| thus
     |h1+l1-sinh(U[j][0])*cosh(w)| < 2^-67.01*|h1+l1| */
  s_mul (&h2, &l2, U[j][2], swh, swl); /* U[j][2]*sinh(w) */
  /* |U[j][2] - cosh(U[j][0])| < 2^-16 ulp(U[j][1]) <= 2^-68 |U[j][1]|
     and |swh + swl - sinh(w)| < 2^-67.99*|swh+swl| thus
     |h2+l2-cosh(U[j][0])*sinh(w)| < 2^-66.99*|h2+l2| */

  fast_sum2 (h, l, h1, l1, h2, l2); /* h+l approximates sinh(v) */
  /* since h1+l1 and h2+l2 have a relative error bound < 2^-66.99, that bound
     holds for the sum of their absolute values, but we might have
     cancellation, the worst case being for j=1 and w=-0.00543,
     where h1+l1 >= 0.0108414. and h2+l2 >= -0.0054304,
     thus (|h1+l1| + |h2+l2|)/((|h1+l1| - |h2+l2|) < 3.008,
     thus |h + l - sinh(v)| < 3.008*2^-66.99 < 2^-65.40.
     Note: the rounding error in fast_sum2() is absorbed in the above
     error bound (which is over-estimated). */

  if (i == 0)
    return 0x1.85p-66 * *h; /* 2^-65.40 < 0x1.85p-66 */

  svh = *h;
  svl = *l;
  s_mul (&h1, &l1, U[j][1], swh, swl); /* U[j][1]*sinh(w) */
  /* |U[j][1] - sinh(U[j][0])| < 2^-16 ulp(U[j][1]) <= 2^-68 |U[j][1]|
     and |swh + swl - sinh(w)| < 2^-67.99*|swh+swl| thus
     |h1+l1-sinh(U[j][0])*sinh(w)| < 2^-66.99*|h1+l1| */
  s_mul (&h2, &l2, U[j][2], cwh, cwl); /* U[j][2]*cosh(w) */
  /* |U[j][2] - cosh(U[j][0])| < 2^-16 ulp(U[j][1]) <= 2^-68 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-68.04*|cwh+cwl| thus
     |h2+l2-cosh(U[j][0])*cosh(w)| < 2^-67.01*|h2+l2| */
  fast_sum2 (&cvh, &cvl, h2, l2, h1, l1); /* cvh+cvl approximates cosh(v) */
  /* since h1+l1 and h2+l2 have a relative error bound < 2^-66.99, that bound
     holds for the sum of their absolute values, but we might have
     cancellation, the worst case being for j=1 and w=-0.00543,
     where h2+l2 >= 1.0000735. and h1+l1 >= -0.0000589
     thus (|h1+l1| + |h2+l2|)/((|h1+l1| - |h2+l2|) < 1.000118,
     and |cvh + cvl - cosh(v)| < 1.000118*2^-66.99 < 2^-66.98.
     Note: he rounding errors in fast_sum2 are absorbed in the above error
     bound (which is over-estimated) */

  /* At this point svh+svl approximates sinh(v) with relative error bounded by
     2^-65.40, cvh+cvl approximates cosh(v) with relative error bounded
     by 2^-66.98, T[i][1] approximates sinh(T[i][0]) with relative error
     bounded by 2^-69, T[i][1]+T[i][2] approximates cosh(T[i][0]) with
     relative error bounded by 2^-69, and we have to compute:
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
     the absolute error on svh+svl is bounded by 2^-65.40*7.95 < 2^-62.40, and
     the absolute error on cvh+cvl is bounded by 2^-66.98*8.02 < 2^-63.97. */

  fast_sum2 (&cvh, &cvl, cvh, cvl, svh, svl);
  /* absolute error on cvh+cvl bounded by (neglecting the error in fast_sum2):
     2^-62.40+2^-63.97 < 2^-61.98.
     Since |v| > 0.00542, the cancellation factor
     (cosh(v)+sinh(v))/(cosh(v)-sinh(v)) is bounded by 1.0109,
     thus the relative error on cvh+cvl is < 1.0109*2^-61.98 < 2^-61.96. */

  s_mul (&h1, &l1, T[i][1], cvh, cvl); /* T[i][1]*(cvh+cvl+svh+svl) */
  /* |T[i][1] - sinh(T[i][0])| < 2^-17 ulp(T[i][1]) <= 2^-69 |T[i][1]|
     and |cvh + cvl - (cosh(v)+sinh(v))| < 2^-61.96*|cvh + cvl| thus
     |h1+l1-sinh(T[i][0])*(cosh(v)+sinh(v))| < 2^-61.94*|h1+l1| */
  s_mul (&h2, &l2, T[i][2], svh, svl); /* T[i][2]*(svh+svl) */
  /* |T[i][2] - exp(T[i][0])| < 2^-17 ulp(T[i][2]) <= 2^-69 |T[i][2]|
     and |svh + svl - sinh(v)| < 2^-65.40*|svh + svl| thus
     |h2+l2-exp(T[i][0])*sinh(v)| < 2^-65.28*|h2+l2| */
  fast_sum2 (h, l, h1, l1, h2, l2);

  /* Warning: h2 might be negative if j=0 and w<0, thus v=w */

  /* 2^-61.96 < 0x1.08p-62 and 2^-65.28 < 0x1.a6p-66 */
  return 0x1.08p-62 * h1 + 0x1.a6p-66 * (h2 > 0 ? h2 : -h2);
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
        return;
      }
  }

  double swh, swl, cwh, cwl;
  swh = *h;
  swl = *l;
  eval_C2 (&cwh, &cwl, w);
  double h1, l1, h2, l2;
  d_mul (&h1, &l1, U[j][1], Ul[j][0], cwh, cwl);
  d_mul (&h2, &l2, U[j][2], Ul[j][1], swh, swl);
  fast_sum2 (h, l, h1, l1, h2, l2);
  if (i == 0)
  {
    static double exceptions[][3] = {
      {0x1.19e03c96f0997p-6, 0x1.19e3cbe7ef607p-6, -0x1.fffffffffffffp-60},
      {0x1.8c154465149ep-8, 0x1.8c15e26bbaa2p-8, -0x1.fffffffffffffp-62},
      {0x1.9147ff03dfb3p-1, 0x1.bba4dc4067a68p-1, 0x1p-54},
      {0x1.9b88da8cd4e51p-3, 0x1.9e4f4a0396a4cp-3, 0x1.fffffffffffffp-57},
      {0x1.e434bb080b343p-5, 0x1.e47ceb971ab09p-5, 0x1.ffffffffffffep-59},
      {0x1.6022fb964b354p-7, 0x1.6024b7c5f26c9p-7, 0x1.d1bf08b11060ep-113},
      {0x1.fae401ab52294p-7, 0x1.fae92e8cd6c2ap-7, 0x1.7e854fd21daf9p-112},
      {0x1.135e31fdd05d3p-5, 0x1.136b78b25cc57p-5, 0x1.71b117cb2c81fp-113},
      {0x1.78f6aa58a8224p-5, 0x1.7918b9deef9c8p-5, 0x1.779e541d49396p-111},
      {0x1.dc00abdea0eaep-5, 0x1.dc4540dcebb2bp-5, 0x1.82c003c1e103dp-111},
      {0x1.249b20103ea68p-4, 0x1.24dada634b148p-4, 0x1.8bb30e252299ap-110},
      {0x1.616cc75d49226p-2, 0x1.687bd068c1c1ep-2, 0x1.9c0f093da51fdp-111},
      {0x1.e6be9678237a2p-2, 0x1.f94840422b64p-2, 0x1.5b943a99b8006p-106},
      {0x1.a3fc7e4dd47d1p-1, 0x1.d4b21ebf542fp-1, 0x1.ded6353f4d832p-107},
      {0x1.aa3b649a96091p-1, 0x1.dd32c5ed1e93p-1, 0x1.89a6dcd3d65ap-106},
      {0x1.dc5059d4e507dp+0, 0x1.9168c60ed5256p+1, 0x1.c584c5cd3ecfep-104},
      {0x1.dd2e3c85d2b62p-8, 0x1.dd2f50d8b91ddp-8, -0x1.728dfac1eb294p-112},
      {0x1.e2b6e387ef5bp-8, 0x1.e2b8019489d28p-8, -0x1.fa0df3fbc6b35p-114},
      {0x1.f887f6af43f73p-8, 0x1.f8893d4c49738p-8, -0x1.0e529966a8342p-111},
      {0x1.c47a6981ae88fp-6, 0x1.c4892321f02f5p-6, -0x1.1d067422c8b6bp-111},
      {0x1.e4720b8441721p-5, 0x1.e4ba57841e459p-5, -0x1.5bf569ebe7ebfp-110},
      {0x1.a00735384ad44p-3, 0x1.a2e533e5b2ea3p-3, -0x1.6fc3324320414p-109},
      {0x1.e90f16eb88c09p-2, 0x1.fbdd4a37760b7p-2, -0x1.f6968f01db399p-108},
      {0x1.2f5d3b178914ap+0, 0x1.7b8516ffd2406p+0, -0x1.27e918302273ep-104},
    };
    for (int i = 0; i < 24; i++)
      if (x == exceptions[i][0])
      {
        *h = exceptions[i][1];
        *l = exceptions[i][2];
        return;
      }
  }

  static double exceptions[][3] = {
    {0x1.1c11f1687d68fp+7, 0x1.e21f461cfa82bp+203, 0x1.ffffffffffffep+149},
    {0x1.7f0046225d651p+1, 0x1.3e11487da075dp+3, -0x1.fffffffffffffp-51},
    {0x1.8560fe96e572bp+1, 0x1.4e65b385b7c53p+3, 0x1.ad8aed0994c2p-101},
    {0x1.bc3c2d0c95f52p+1, 0x1.00fef7383a978p+4, 0x1.61182ae91b723p-100},
    {0x1.5ce42e9d0df03p+2, 0x1.d22c2e0cd6bf4p+6, 0x1.4b4a00c36a2adp-97},
    {0x1.c95ba20925c4bp+2, 0x1.3d52e798431b6p+9, 0x1.b34b670026b37p-95},
    {0x1.0a19aebb51e9p+3, 0x1.fee8f69c4cd25p+10, 0x1.48bbf8f59b3e9p-95},
    {0x1.1a3441bc1a6b7p+3, 0x1.a68ae777d2cc9p+11, 0x1.cf1f898fa67ebp-92},
    {0x1.20e29ea8b51e2p+4, 0x1.08b8abba28abcp+25, 0x1.9b157420749bdp-79},
    {0x1.c089fcf166171p+4, 0x1.5c452e0e37569p+39, 0x1.3bf07320cd829p-69},
    {0x1.e42a98b3a0be5p+4, 0x1.938768ca4f8aap+42, 0x1.6d39a3108f2b9p-62},
    {0x1.04db52248cbb8p+5, 0x1.0794072349523p+46, 0x1.0eaf76f95f673p-57},
    {0x1.21bc021eeb97ep+5, 0x1.3065064a170fbp+51, 0x1.084083e9e9dd2p-52},
    {0x1.39fc4d3bb711p+5, 0x1.8a4e90733b95ep+55, 0x1.6e767ed14e448p-50},
    {0x1.51d0f4f0a901cp+5, 0x1.e4a01c9ddbc87p+59, 0x1.676f36e53deedp-45},
    {0x1.6ce8416ec0a3fp+5, 0x1.bfa6fb8bb994fp+64, 0x1.2b1c5ece9e743p-39},
    {0x1.94925476814e9p+5, 0x1.f1b76b88f075p+71, 0x1.b710ac550e46ep-32},
    {0x1.96d81955b1fd7p+5, 0x1.4a9d153103f65p+72, 0x1.776d7bf831017p-32},
    {0x1.638099048d4f6p+6, 0x1.2a3f3c4fdc462p+127, 0x1.daed2d108ae4ap+23},
    {0x1.1f0da93354198p+7, 0x1.0bd73b73fc74cp+206, 0x1.588526e93304cp+102},
    {0x1.226b70c1a9d7bp+7, 0x1.686ab849bc518p+208, 0x1.658e06041ef5p+105},
    {0x1.0bc04af1b09f5p+9, 0x1.7b1d97c902985p+771, 0x1.551dfecc05bd4p+666},
    {0x1.8c0a26d055288p+1, 0x1.6056b06a21918p+3, -0x1.be3740fe7b06dp-102},
    {0x1.326e2606c8c86p+3, 0x1.c26eeb4c2c18ap+12, -0x1.5ee2f66583a45p-90},
    {0x1.16369cd53bb69p+4, 0x1.0fbc6c02b1c9p+24, -0x1.90094b292b292p-81},
    {0x1.4e15feac0d3a6p+4, 0x1.16fa4aa224fdap+29, -0x1.37f58912f1c6fp-74},
    {0x1.a1e4f11b513d7p+4, 0x1.9a65b6c2e2185p+36, -0x1.bbe5072d0f2f1p-70},
    {0x1.cb08ac3f4cdf7p+4, 0x1.4f8be4f84d663p+40, -0x1.de033ff070d6bp-63},
    {0x1.3c895d86e96c9p+5, 0x1.0f33837882a6p+56, -0x1.277dbae956fc4p-49},
    {0x1.5140272441b93p+5, 0x1.c38b159d44744p+59, -0x1.56d2d5f6a1678p-43},
    {0x1.6fc71838701e6p+5, 0x1.406f375086cc1p+65, -0x1.644d8f437d059p-39},
    {0x1.48a6374622156p+6, 0x1.72f930dc88767p+117, -0x1.a66033893e8a6p+13},
    {0x1.6474c604cc0d7p+6, 0x1.7a8f65ad009bdp+127, -0x1.0b611158ec877p+20},
    {0x1.07fa29aa7f36cp+8, 0x1.c9ce85bf1aaf8p+379, -0x1.044fc6419c125p+278},
    {0x1.1debbc1b6dd4cp+8, 0x1.692e3b58bcc99p+411, -0x1.bccab9bd92516p+308},
    {0x1.204684c1167e9p+8, 0x1.db9797d3d32e8p+414, -0x1.51e78c6bad663p+310},
    {0x1.01ee19aead26ap+9, 0x1.2c03f281caee1p+743, -0x1.e313f0eaa37cap+640},
  };
  for (int i = 0; i < 37; i++)
    if (x == exceptions[i][0])
    {
      *h = exceptions[i][1];
      *l = exceptions[i][2];
      return;
    }

  double svh, svl, cvh, cvl;
  svh = *h;
  svl = *l;
  /* svh+svl approximates sinh(v) */
  d_mul (&h1, &l1, U[j][1], Ul[j][0], swh, swl);
  d_mul (&h2, &l2, U[j][2], Ul[j][1], cwh, cwl);
  fast_sum2 (&cvh, &cvl, h2, l2, h1, l1); /* cvh+cvl approximates cosh(v) */
  fast_sum2 (&cvh, &cvl, cvh, cvl, svh, svl);
  /* now cvh+cvl approximates cosh(v)+sinh(v) */
  d_mul (&h1, &l1, T[i][1], Tl[i][0], cvh, cvl);
  d_mul (&h2, &l2, T[i][2], Tl[i][1], svh, svl);
  fast_sum2 (h, l, h1, l1, h2, l2);
  //*h = *l = 0;
}

#define MASK 0x7fffffffffffffff /* to mask the sign bit */

#if 0
static void
printT ()
{
  printf ("LOG T0=dict()\n");
  printf ("LOG T1=dict()\n");
  printf ("LOG T2=dict()\n");
  for (int i = 0; i < 256; i++)
  {
    printf ("LOG T0[%d]='%la'\n", i, T[i][0]);
    printf ("LOG T1[%d]='%la'\n", i, T[i][1]);
    printf ("LOG T2[%d]='%la'\n", i, T[i][2]);
  }
}

static void
printU ()
{
  printf ("LOG U0=dict()\n");
  printf ("LOG U1=dict()\n");
  printf ("LOG U2=dict()\n");
  for (int i = 0; i < 256; i++)
  {
    printf ("LOG U0[%d]='%la'\n", i, U[i][0]);
    printf ("LOG U1[%d]='%la'\n", i, U[i][1]);
    printf ("LOG U2[%d]='%la'\n", i, U[i][2]);
  }
}
#endif

double
cr_sinh (double x)
{
  d64u64 v = {.f = x};
  int e = (v.u >> 52) - 0x3ff;
  int s = v.u >> 63; /* sign bit */
  v.u &= (uint64_t) MASK; /* get absolute value */

#if 0
  static int count = 0;
  if (count++ == 0) printU();
  //if (count++ == 0) printT();
#endif

  if (e == 0x400 || e == 0xc00 || v.f >= 0x1.633ce8fb9f87ep+9)
    /* NaN or overflow */
  {
    /* this will return NaN for x=NaN, Inf with the correct sign for x=+Inf
       or -Inf, and Inf/DBL_MAX or -Inf/-DBL_MAX for other |x| >= 2^10 */
    return x * 0x1p1023;
  }
  
  double h, l;
  double err = cr_sinh_fast (&h, &l, v.f);
  double sign[] = { 1.0, -1.0 };
  h *= sign[s];
  l *= sign[s];
  
  double left  = h + (l - err);
  double right = h + (l + err);
  if (left == right)
    return left;

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
