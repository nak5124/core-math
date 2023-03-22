/* Correctly rounded tanh for binary64 values.

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

/* Both in the fast path and in the accurate path, sinh(x) and cosh(x) are
   approximated with double-double representations sh+sl and ch+cl, which
   are divided to obtain an approximation of tanh(x).

   For the error analysis in the approximations for sinh(x) and cosh(x),
   see ../sinh/sinh.c.

   References:
   [1] Modern Computer Arithmetic, Richard Brent and Paul Zimmermann,
       Cambridge University Press, 2011.
   [2] Note on FastTwoSum with Directed Roundings, Paul Zimmermann,
       https://hal.inria.fr/hal-03798376, 2022.
   [3] Tight and rigourous error bounds for basic building blocks of
       double-word arithmetic, Mioara Maria Joldes, Jean-Michel Muller,
       Valentina Popescu, ACM Transactions on Mathematical Software, Vol. 44,
       No. 2, Article 15res. Publication date: October 2017,
       https://hal.science/hal-01351529v3/document
*/       

#include <stdint.h>

/* For 0 <= i < 256, T[i] = {xi, shi, sli, chi, cli} such that xi is near
   i*2^8/magic with magic = 0x1.70f77fc88ae3cp6, and shi+sli, chi+cli
   approximate sinh(xi), cosh(xi) with accuracy >= 107 bits:
   |shi + sli - sinh(xi)| < 2^-55 ulp(shi),
   |chi + cli - cosh(xi)| < 2^-55 ulp(chi).
   Thus shi+sli approximates sinh(xi) with relative error < 2^-55 ulp(shi)/shi
   <= 2^-107 (same for ci).
   We have |xi - i*2^8/magic| < 8e-14.
   Generated with build_table_T(k=16) from accompanying file sinh.sage.
*/
static const double T[256][5] = {
   {0x0p+0, 0x0p+0, 0x0p+0,
    0x1p+0, 0x0p+0}, /* i=0 */
   {0x1.633d9a9a6cd52p+1, 0x1.ff678ce3789bcp+2, -0x1.3828d50eaed42p-52,
    0x1.01b26143e173dp+3, 0x1.a19bc22962e2ep-52}, /* i=1 */
   {0x1.633d9a9a6cd52p+2, 0x1.0165a65f148ddp+7, 0x1.efc0742ecaf16p-48,
    0x1.0167a395b2049p+7, 0x1.86138dcd4e4b2p-47}, /* i=2 */
   {0x1.0a6e33f3d19fep+3, 0x1.021ab2b87456p+11, -0x1.86337b457b81p-44,
    0x1.021ab4b447b74p+11, 0x1.fa483834229acp-44}, /* i=3 */
   {0x1.633d9a9a6cd52p+3, 0x1.02cf40659ff95p+15, 0x1.7e64e1f4905fp-43,
    0x1.02cf40679a6a7p+15, -0x1.5bef911501188p-40}, /* i=4 */
   {0x1.bc0d0141080a6p+3, 0x1.03844b633b5f9p+19, -0x1.455244d0d16bfp-35,
    0x1.03844b633d58ap+19, -0x1.5166a7a806dc5p-35}, /* i=5 */
   {0x1.0a6e33f3d19fep+4, 0x1.0439d50494b8p+23, 0x1.28f6a11801454p-32,
    0x1.0439d50494bap+23, -0x1.7fa6745c448fep-31}, /* i=6 */
   {0x1.36d5e7471f3a8p+4, 0x1.04efdda33d469p+27, -0x1.c34e3f37ba498p-29,
    0x1.04efdda33d469p+27, 0x1.980ec5f4dbd8p-32}, /* i=7 */
   {0x1.633d9a9a6cd52p+4, 0x1.05a665980afe3p+31, 0x1.c17eaf6b12cc6p-24,
    0x1.05a665980afe3p+31, 0x1.c2792841233c8p-24}, /* i=8 */
   {0x1.8fa54dedba6fcp+4, 0x1.065d6d3c10fc2p+35, 0x1.3f22b831c4824p-20,
    0x1.065d6d3c10fc2p+35, 0x1.3f23b1fbdee2cp-20}, /* i=9 */
   {0x1.bc0d0141080a6p+4, 0x1.0714f4e8a0abep+39, -0x1.c836b08872d2cp-15,
    0x1.0714f4e8a0abep+39, -0x1.c836b00be4e68p-15}, /* i=10 */
   {0x1.e874b49455a51p+4, 0x1.07ccfcf749f4bp+43, 0x1.ada8e74bd340ep-11,
    0x1.07ccfcf749f4bp+43, 0x1.ada8e74c4f77ep-11}, /* i=11 */
   {0x1.0a6e33f3d19fep+5, 0x1.088585c1db62ap+47, -0x1.c95d9438b1c7p-9,
    0x1.088585c1db62ap+47, -0x1.c95d9438afd78p-9}, /* i=12 */
   {0x1.20a20d9d786d3p+5, 0x1.093e8fa26253dp+51, 0x1.26b47c5e541dap-4,
    0x1.093e8fa26253dp+51, 0x1.26b47c5e541eap-4}, /* i=13 */
   {0x1.36d5e7471f3a8p+5, 0x1.09f81af32b27p+55, -0x1.ba1cc1f1e81c6p+0,
    0x1.09f81af32b27p+55, -0x1.ba1cc1f1e81c6p+0}, /* i=14 */
   {0x1.4d09c0f0c607dp+5, 0x1.0ab2280ec1642p+59, 0x1.41a722f452ccfp+5,
    0x1.0ab2280ec1642p+59, 0x1.41a722f452ccfp+5}, /* i=15 */
   {0x1.633d9a9a6cd52p+5, 0x1.0b6cb74fefe9fp+63, -0x1.941da02a1bda9p+9,
    0x1.0b6cb74fefe9fp+63, -0x1.941da02a1bda9p+9}, /* i=16 */
   {0x1.7971744413a27p+5, 0x1.0c27c911c119dp+67, -0x1.9b5f1705f0998p+12,
    0x1.0c27c911c119dp+67, -0x1.9b5f1705f0998p+12}, /* i=17 */
   {0x1.8fa54dedba6fcp+5, 0x1.0ce35daf7f04ap+71, 0x1.2d1a48572d5f2p+16,
    0x1.0ce35daf7f04ap+71, 0x1.2d1a48572d5f2p+16}, /* i=18 */
   {0x1.a5d92797613d1p+5, 0x1.0d9f7584b3972p+75, -0x1.a66b7b10c6f3ep+20,
    0x1.0d9f7584b3972p+75, -0x1.a66b7b10c6f3ep+20}, /* i=19 */
   {0x1.bc0d0141080a6p+5, 0x1.0e5c10ed28c68p+79, -0x1.f81b3856784bp+25,
    0x1.0e5c10ed28c68p+79, -0x1.f81b3856784bp+25}, /* i=20 */
   {0x1.d240daeaaed7cp+5, 0x1.0f193044e8bf6p+83, 0x1.01f2d564e855ap+29,
    0x1.0f193044e8bf6p+83, 0x1.01f2d564e855ap+29}, /* i=21 */
   {0x1.e874b49455a51p+5, 0x1.0fd6d3e83e0a7p+87, 0x1.7a02a3d718719p+33,
    0x1.0fd6d3e83e0a7p+87, 0x1.7a02a3d718719p+33}, /* i=22 */
   {0x1.fea88e3dfc726p+5, 0x1.1094fc33b3c5ep+91, 0x1.df30bfbaa6232p+36,
    0x1.1094fc33b3c5ep+91, 0x1.df30bfbaa6232p+36}, /* i=23 */
   {0x1.0a6e33f3d19fep+6, 0x1.1153a98415cc6p+95, -0x1.1d7b48ef46ec8p+41,
    0x1.1153a98415cc6p+95, -0x1.1d7b48ef46ec8p+41}, /* i=24 */
   {0x1.158820c8a5068p+6, 0x1.1212dc3670d9cp+99, -0x1.1df585fc6343bp+45,
    0x1.1212dc3670d9cp+99, -0x1.1df585fc6343bp+45}, /* i=25 */
   {0x1.20a20d9d786d3p+6, 0x1.12d294a812dp+103, -0x1.49e89850cf114p+49,
    0x1.12d294a812dp+103, -0x1.49e89850cf114p+49}, /* i=26 */
   {0x1.2bbbfa724bd3dp+6, 0x1.1392d3368ac4ap+107, -0x1.42fc0dfbf4baap+53,
    0x1.1392d3368ac4ap+107, -0x1.42fc0dfbf4baap+53}, /* i=27 */
   {0x1.36d5e7471f3a8p+6, 0x1.1453983fa950ap+111, -0x1.de71ac88bb8c4p+55,
    0x1.1453983fa950ap+111, -0x1.de71ac88bb8c4p+55}, /* i=28 */
   {0x1.41efd41bf2a12p+6, 0x1.1514e421809bfp+115, -0x1.b9ba97ba36149p+61,
    0x1.1514e421809bfp+115, -0x1.b9ba97ba36149p+61}, /* i=29 */
   {0x1.4d09c0f0c607dp+6, 0x1.15d6b73a64adbp+119, -0x1.a070375e7ebc2p+65,
    0x1.15d6b73a64adbp+119, -0x1.a070375e7ebc2p+65}, /* i=30 */
   {0x1.5823adc5996e7p+6, 0x1.169911e8eb77cp+123, 0x1.ddb98520889a6p+69,
    0x1.169911e8eb77cp+123, 0x1.ddb98520889a6p+69}, /* i=31 */
   {0x1.633d9a9a6cd52p+6, 0x1.175bf48bed27cp+127, -0x1.9c25e349755p+72,
    0x1.175bf48bed27cp+127, -0x1.9c25e349755p+72}, /* i=32 */
   {0x1.6e57876f403bdp+6, 0x1.181f5f8284367p+131, 0x1.bc2f65bdd84f1p+77,
    0x1.181f5f8284367p+131, 0x1.bc2f65bdd84f1p+77}, /* i=33 */
   {0x1.7971744413a27p+6, 0x1.18e3532c0da3cp+135, -0x1.743f1a9ab7ccep+80,
    0x1.18e3532c0da3cp+135, -0x1.743f1a9ab7ccep+80}, /* i=34 */
   {0x1.848b6118e7092p+6, 0x1.19a7cfe82931fp+139, -0x1.5a8db7058da8p+84,
    0x1.19a7cfe82931fp+139, -0x1.5a8db7058da8p+84}, /* i=35 */
   {0x1.8fa54dedba6fcp+6, 0x1.1a6cd616b975fp+143, 0x1.7a1eb280a6538p+86,
    0x1.1a6cd616b975fp+143, 0x1.7a1eb280a6538p+86}, /* i=36 */
   {0x1.9abf3ac28dd67p+6, 0x1.1b326617e4293p+147, -0x1.a7e69beb2ad2p+90,
    0x1.1b326617e4293p+147, -0x1.a7e69beb2ad2p+90}, /* i=37 */
   {0x1.a5d92797613d1p+6, 0x1.1bf8804c12354p+151, -0x1.c200c99029c97p+97,
    0x1.1bf8804c12354p+151, -0x1.c200c99029c97p+97}, /* i=38 */
   {0x1.b0f3146c34a3cp+6, 0x1.1cbf2513f0063p+155, -0x1.528d38b4512fcp+100,
    0x1.1cbf2513f0063p+155, -0x1.528d38b4512fcp+100}, /* i=39 */
   {0x1.bc0d0141080a6p+6, 0x1.1d8654d06d967p+159, 0x1.359180d82ca4p+105,
    0x1.1d8654d06d967p+159, 0x1.359180d82ca4p+105}, /* i=40 */
   {0x1.c726ee15db711p+6, 0x1.1e4e0fe2bec1ap+163, -0x1.49434ae5bb604p+109,
    0x1.1e4e0fe2bec1ap+163, -0x1.49434ae5bb604p+109}, /* i=41 */
   {0x1.d240daeaaed7cp+6, 0x1.1f1656ac5b549p+167, 0x1.3f92341a83c4bp+113,
    0x1.1f1656ac5b549p+167, 0x1.3f92341a83c4bp+113}, /* i=42 */
   {0x1.dd5ac7bf823e6p+6, 0x1.1fdf298eff4aap+171, 0x1.7a84561098768p+116,
    0x1.1fdf298eff4aap+171, 0x1.7a84561098768p+116}, /* i=43 */
   {0x1.e874b49455a51p+6, 0x1.20a888ecab0a9p+175, -0x1.faf655e2cd50dp+121,
    0x1.20a888ecab0a9p+175, -0x1.faf655e2cd50dp+121}, /* i=44 */
   {0x1.f38ea169290bbp+6, 0x1.21727527a376ep+179, 0x1.45a3cb8401047p+125,
    0x1.21727527a376ep+179, 0x1.45a3cb8401047p+125}, /* i=45 */
   {0x1.fea88e3dfc726p+6, 0x1.223ceea272423p+183, 0x1.48a6b15df901dp+129,
    0x1.223ceea272423p+183, 0x1.48a6b15df901dp+129}, /* i=46 */
   {0x1.04e13d8967ec8p+7, 0x1.2307f5bfe5facp+187, 0x1.0f46097010712p+132,
    0x1.2307f5bfe5facp+187, 0x1.0f46097010712p+132}, /* i=47 */
   {0x1.0a6e33f3d19fep+7, 0x1.23d38ae31263ap+191, 0x1.ab894df7bbcf7p+137,
    0x1.23d38ae31263ap+191, 0x1.ab894df7bbcf7p+137}, /* i=48 */
   {0x1.0ffb2a5e3b533p+7, 0x1.249fae6f506eap+195, -0x1.b065f457f4af6p+141,
    0x1.249fae6f506eap+195, -0x1.b065f457f4af6p+141}, /* i=49 */
   {0x1.158820c8a5068p+7, 0x1.256c60c83ea79p+199, 0x1.97e0c7064f364p+143,
    0x1.256c60c83ea79p+199, 0x1.97e0c7064f364p+143}, /* i=50 */
   {0x1.1b1517330eb9dp+7, 0x1.2639a251c141p+203, -0x1.51144628877d8p+147,
    0x1.2639a251c141p+203, -0x1.51144628877d8p+147}, /* i=51 */
   {0x1.20a20d9d786d3p+7, 0x1.270773700256dp+207, 0x1.91f7f1ad3e494p+153,
    0x1.270773700256dp+207, 0x1.91f7f1ad3e494p+153}, /* i=52 */
   {0x1.262f0407e2208p+7, 0x1.27d5d48771faap+211, -0x1.1081a0dad54e9p+157,
    0x1.27d5d48771faap+211, -0x1.1081a0dad54e9p+157}, /* i=53 */
   {0x1.2bbbfa724bd3dp+7, 0x1.28a4c5fcc69bap+215, 0x1.8805e211fa874p+159,
    0x1.28a4c5fcc69bap+215, 0x1.8805e211fa874p+159}, /* i=54 */
   {0x1.3148f0dcb5872p+7, 0x1.29744834fd136p+219, 0x1.56dfd93bf425cp+163,
    0x1.29744834fd136p+219, 0x1.56dfd93bf425cp+163}, /* i=55 */
   {0x1.36d5e7471f3a8p+7, 0x1.2a445b9558e95p+223, -0x1.f7267a2140983p+169,
    0x1.2a445b9558e95p+223, -0x1.f7267a2140983p+169}, /* i=56 */
   {0x1.3c62ddb188eddp+7, 0x1.2b150083645eep+227, -0x1.4089dded17304p+171,
    0x1.2b150083645eep+227, -0x1.4089dded17304p+171}, /* i=57 */
   {0x1.41efd41bf2a12p+7, 0x1.2be63764f0d93p+231, 0x1.1352f4314cbcp+173,
    0x1.2be63764f0d93p+231, 0x1.1352f4314cbcp+173}, /* i=58 */
   {0x1.477cca865c548p+7, 0x1.2cb800a016f6bp+235, -0x1.576dc4f2a2b89p+181,
    0x1.2cb800a016f6bp+235, -0x1.576dc4f2a2b89p+181}, /* i=59 */
   {0x1.4d09c0f0c607dp+7, 0x1.2d8a5c9b36a4bp+239, -0x1.cdaed92747ef2p+184,
    0x1.2d8a5c9b36a4bp+239, -0x1.cdaed92747ef2p+184}, /* i=60 */
   {0x1.5296b75b2fbb2p+7, 0x1.2e5d4bbcf789dp+243, -0x1.c3c6eb63c6aa1p+189,
    0x1.2e5d4bbcf789dp+243, -0x1.c3c6eb63c6aa1p+189}, /* i=61 */
   {0x1.5823adc5996e7p+7, 0x1.2f30ce6c49126p+247, -0x1.35e6ffff45506p+192,
    0x1.2f30ce6c49126p+247, -0x1.35e6ffff45506p+192}, /* i=62 */
   {0x1.5db0a4300321dp+7, 0x1.3004e51062b5bp+251, -0x1.dbe869eaf9f4ep+196,
    0x1.3004e51062b5bp+251, -0x1.dbe869eaf9f4ep+196}, /* i=63 */
   {0x1.633d9a9a6cd52p+7, 0x1.30d99010c4024p+255, 0x1.2a96165ae09ap+197,
    0x1.30d99010c4024p+255, 0x1.2a96165ae09ap+197}, /* i=64 */
   {0x1.68ca9104d6887p+7, 0x1.31aecfd535096p+259, -0x1.8bfa96a8845cp+204,
    0x1.31aecfd535096p+259, -0x1.8bfa96a8845cp+204}, /* i=65 */
   {0x1.6e57876f403bdp+7, 0x1.3284a4c5c6755p+263, -0x1.7d66c60a35885p+209,
    0x1.3284a4c5c6755p+263, -0x1.7d66c60a35885p+209}, /* i=66 */
   {0x1.73e47dd9a9ef2p+7, 0x1.335b0f4ad19f5p+267, -0x1.ca42095690694p+212,
    0x1.335b0f4ad19f5p+267, -0x1.ca42095690694p+212}, /* i=67 */
   {0x1.7971744413a27p+7, 0x1.34320fccf8fbfp+271, 0x1.9c684f043afcfp+217,
    0x1.34320fccf8fbfp+271, 0x1.9c684f043afcfp+217}, /* i=68 */
   {0x1.7efe6aae7d55cp+7, 0x1.3509a6b52828p+275, -0x1.2885a7ee503bp+221,
    0x1.3509a6b52828p+275, -0x1.2885a7ee503bp+221}, /* i=69 */
   {0x1.848b6118e7092p+7, 0x1.35e1d46c942eap+279, -0x1.891f4641a6e18p+222,
    0x1.35e1d46c942eap+279, -0x1.891f4641a6e18p+222}, /* i=70 */
   {0x1.8a18578350bc7p+7, 0x1.36ba995cbb966p+283, -0x1.d67f99ada27f9p+229,
    0x1.36ba995cbb966p+283, -0x1.d67f99ada27f9p+229}, /* i=71 */
   {0x1.8fa54dedba6fcp+7, 0x1.3793f5ef66ce5p+287, 0x1.f93a72b774884p+231,
    0x1.3793f5ef66ce5p+287, 0x1.f93a72b774884p+231}, /* i=72 */
   {0x1.9532445824232p+7, 0x1.386dea8ea8452p+291, 0x1.fb38bc3d91bc4p+236,
    0x1.386dea8ea8452p+291, 0x1.fb38bc3d91bc4p+236}, /* i=73 */
   {0x1.9abf3ac28dd67p+7, 0x1.394877a4dc7f4p+295, -0x1.64a4989c2f603p+241,
    0x1.394877a4dc7f4p+295, -0x1.64a4989c2f603p+241}, /* i=74 */
   {0x1.a04c312cf789cp+7, 0x1.3a239d9caa854p+299, -0x1.0986b45ff894bp+245,
    0x1.3a239d9caa854p+299, -0x1.0986b45ff894bp+245}, /* i=75 */
   {0x1.a5d92797613d1p+7, 0x1.3aff5ce103f12p+303, -0x1.075327d578db6p+249,
    0x1.3aff5ce103f12p+303, -0x1.075327d578db6p+249}, /* i=76 */
   {0x1.ab661e01caf07p+7, 0x1.3bdbb5dd2535ep+307, 0x1.c21564c826627p+253,
    0x1.3bdbb5dd2535ep+307, 0x1.c21564c826627p+253}, /* i=77 */
   {0x1.b0f3146c34a3cp+7, 0x1.3cb8a8fc95acap+311, 0x1.645304c0d3bdcp+256,
    0x1.3cb8a8fc95acap+311, 0x1.645304c0d3bdcp+256}, /* i=78 */
   {0x1.b6800ad69e571p+7, 0x1.3d9636ab2803fp+315, -0x1.24b6e78658555p+261,
    0x1.3d9636ab2803fp+315, -0x1.24b6e78658555p+261}, /* i=79 */
   {0x1.bc0d0141080a6p+7, 0x1.3e745f54fa4d3p+319, -0x1.2ee5edcf4666dp+265,
    0x1.3e745f54fa4d3p+319, -0x1.2ee5edcf4666dp+265}, /* i=80 */
   {0x1.c199f7ab71bdcp+7, 0x1.3f53236676453p+323, 0x1.ef5f9c222aap+261,
    0x1.3f53236676453p+323, 0x1.ef5f9c222aap+261}, /* i=81 */
   {0x1.c726ee15db711p+7, 0x1.4032834c51613p+327, -0x1.a8bfc4c8f24bp+270,
    0x1.4032834c51613p+327, -0x1.a8bfc4c8f24bp+270}, /* i=82 */
   {0x1.ccb3e48045246p+7, 0x1.41127f738d3fap+331, -0x1.5cd830277bc3p+276,
    0x1.41127f738d3fap+331, -0x1.5cd830277bc3p+276}, /* i=83 */
   {0x1.d240daeaaed7cp+7, 0x1.41f3184977bfap+335, -0x1.7f23abf878cd3p+281,
    0x1.41f3184977bfap+335, -0x1.7f23abf878cd3p+281}, /* i=84 */
   {0x1.d7cdd155188b1p+7, 0x1.42d44e3bab18p+339, 0x1.bdddd4ecda5aep+284,
    0x1.42d44e3bab18p+339, 0x1.bdddd4ecda5aep+284}, /* i=85 */
   {0x1.dd5ac7bf823e6p+7, 0x1.43b621b80e494p+343, -0x1.3c2484f9abbcbp+289,
    0x1.43b621b80e494p+343, -0x1.3c2484f9abbcbp+289}, /* i=86 */
   {0x1.e2e7be29ebf1bp+7, 0x1.4498932cd52aap+347, -0x1.10abd57fb28eap+292,
    0x1.4498932cd52aap+347, -0x1.10abd57fb28eap+292}, /* i=87 */
   {0x1.e874b49455a51p+7, 0x1.457ba30880b48p+351, 0x1.0ddf6a48ab168p+294,
    0x1.457ba30880b48p+351, 0x1.0ddf6a48ab168p+294}, /* i=88 */
   {0x1.ee01aafebf586p+7, 0x1.465f51b9df0d7p+355, 0x1.db9be8b6a22e4p+300,
    0x1.465f51b9df0d7p+355, 0x1.db9be8b6a22e4p+300}, /* i=89 */
   {0x1.f38ea169290bbp+7, 0x1.47439fb00bfd5p+359, -0x1.d449285c07865p+305,
    0x1.47439fb00bfd5p+359, -0x1.d449285c07865p+305}, /* i=90 */
   {0x1.f91b97d392bf1p+7, 0x1.48288d5a7104fp+363, 0x1.a63926d87f85p+306,
    0x1.48288d5a7104fp+363, 0x1.a63926d87f85p+306}, /* i=91 */
   {0x1.fea88e3dfc726p+7, 0x1.490e1b28c576p+367, 0x1.799c3b481b118p+312,
    0x1.490e1b28c576p+367, 0x1.799c3b481b118p+312}, /* i=92 */
   {0x1.021ac2543312ep+8, 0x1.49f4498b0ef1p+371, 0x1.71dc1b85b15ccp+315,
    0x1.49f4498b0ef1p+371, 0x1.71dc1b85b15ccp+315}, /* i=93 */
   {0x1.04e13d8967ec8p+8, 0x1.4adb18f1a13f8p+375, 0x1.c4c09d524ba1bp+321,
    0x1.4adb18f1a13f8p+375, 0x1.c4c09d524ba1bp+321}, /* i=94 */
   {0x1.07a7b8be9cc63p+8, 0x1.4bc289cd1f06dp+379, -0x1.d74192e59cfbp+324,
    0x1.4bc289cd1f06dp+379, -0x1.d74192e59cfbp+324}, /* i=95 */
   {0x1.0a6e33f3d19fep+8, 0x1.4caa9c8e79787p+383, 0x1.9c6e5b78879dap+329,
    0x1.4caa9c8e79787p+383, 0x1.9c6e5b78879dap+329}, /* i=96 */
   {0x1.0d34af2906798p+8, 0x1.4d9351a6f0c76p+387, 0x1.d6a14d9bd4263p+333,
    0x1.4d9351a6f0c76p+387, 0x1.d6a14d9bd4263p+333}, /* i=97 */
   {0x1.0ffb2a5e3b533p+8, 0x1.4e7ca988149dcp+391, -0x1.d285544bd83fdp+337,
    0x1.4e7ca988149dcp+391, -0x1.d285544bd83fdp+337}, /* i=98 */
   {0x1.12c1a593702cdp+8, 0x1.4f66a4a3c3c24p+395, 0x1.3382d5af02283p+341,
    0x1.4f66a4a3c3c24p+395, 0x1.3382d5af02283p+341}, /* i=99 */
   {0x1.158820c8a5068p+8, 0x1.5051436c2cf74p+399, -0x1.dbee7be2c403ap+344,
    0x1.5051436c2cf74p+399, -0x1.dbee7be2c403ap+344}, /* i=100 */
   {0x1.184e9bfdd9e03p+8, 0x1.513c8653ce9f8p+403, -0x1.a73a863311113p+349,
    0x1.513c8653ce9f8p+403, -0x1.a73a863311113p+349}, /* i=101 */
   {0x1.1b1517330eb9dp+8, 0x1.52286dcd7734fp+407, -0x1.58ac54a3b233p+353,
    0x1.52286dcd7734fp+407, -0x1.58ac54a3b233p+353}, /* i=102 */
   {0x1.1ddb926843938p+8, 0x1.5314fa4c45c03p+411, -0x1.a062730398124p+356,
    0x1.5314fa4c45c03p+411, -0x1.a062730398124p+356}, /* i=103 */
   {0x1.20a20d9d786d3p+8, 0x1.54022c43a9921p+415, 0x1.6d5b8bef2c8e2p+361,
    0x1.54022c43a9921p+415, 0x1.6d5b8bef2c8e2p+361}, /* i=104 */
   {0x1.236888d2ad46dp+8, 0x1.54f0042762bb2p+419, -0x1.e68ffefaccf83p+365,
    0x1.54f0042762bb2p+419, -0x1.e68ffefaccf83p+365}, /* i=105 */
   {0x1.262f0407e2208p+8, 0x1.55de826b8283ep+423, -0x1.662177e76a532p+369,
    0x1.55de826b8283ep+423, -0x1.662177e76a532p+369}, /* i=106 */
   {0x1.28f57f3d16fa3p+8, 0x1.56cda7846b263p+427, -0x1.598b7b7b2fbfbp+373,
    0x1.56cda7846b263p+427, -0x1.598b7b7b2fbfbp+373}, /* i=107 */
   {0x1.2bbbfa724bd3dp+8, 0x1.57bd73e6d0456p+431, 0x1.a537c6647a82p+372,
    0x1.57bd73e6d0456p+431, 0x1.a537c6647a82p+372}, /* i=108 */
   {0x1.2e8275a780ad8p+8, 0x1.58ade807b767fp+435, -0x1.c006742a5f0ep+379,
    0x1.58ade807b767fp+435, -0x1.c006742a5f0ep+379}, /* i=109 */
   {0x1.3148f0dcb5872p+8, 0x1.599f045c779a3p+439, -0x1.506143d73caebp+385,
    0x1.599f045c779a3p+439, -0x1.506143d73caebp+385}, /* i=110 */
   {0x1.340f6c11ea60dp+8, 0x1.5a90c95aba53bp+443, -0x1.a207b55edd132p+388,
    0x1.5a90c95aba53bp+443, -0x1.a207b55edd132p+388}, /* i=111 */
   {0x1.36d5e7471f3a8p+8, 0x1.5b8337787b19ep+447, -0x1.e465e55f9cf76p+392,
    0x1.5b8337787b19ep+447, -0x1.e465e55f9cf76p+392}, /* i=112 */
   {0x1.399c627c54142p+8, 0x1.5c764f2c07fap+451, -0x1.ded94eada670ap+397,
    0x1.5c764f2c07fap+451, -0x1.ded94eada670ap+397}, /* i=113 */
   {0x1.3c62ddb188eddp+8, 0x1.5d6a10ec02045p+455, 0x1.d4ea73ccb8f04p+401,
    0x1.5d6a10ec02045p+455, 0x1.d4ea73ccb8f04p+401}, /* i=114 */
   {0x1.3f2958e6bdc78p+8, 0x1.5e5e7d2f5d03dp+459, -0x1.85723b3df9002p+405,
    0x1.5e5e7d2f5d03dp+459, -0x1.85723b3df9002p+405}, /* i=115 */
   {0x1.41efd41bf2a12p+8, 0x1.5f53946d5ff8ap+463, 0x1.1378e8c69bf44p+409,
    0x1.5f53946d5ff8ap+463, 0x1.1378e8c69bf44p+409}, /* i=116 */
   {0x1.44b64f51277adp+8, 0x1.6049571da595p+467, -0x1.74df8eee39c7p+409,
    0x1.6049571da595p+467, -0x1.74df8eee39c7p+409}, /* i=117 */
   {0x1.477cca865c548p+8, 0x1.613fc5b81bf38p+471, 0x1.a796a5c3acb0ep+417,
    0x1.613fc5b81bf38p+471, 0x1.a796a5c3acb0ep+417}, /* i=118 */
   {0x1.4a4345bb912e2p+8, 0x1.6236e0b505138p+475, -0x1.aeada4efec1ap+417,
    0x1.6236e0b505138p+475, -0x1.aeada4efec1ap+417}, /* i=119 */
   {0x1.4d09c0f0c607dp+8, 0x1.632ea88cf7561p+479, 0x1.1405492a0c49dp+425,
    0x1.632ea88cf7561p+479, 0x1.1405492a0c49dp+425}, /* i=120 */
   {0x1.4fd03c25fae18p+8, 0x1.64271db8dd348p+483, 0x1.1f7341c9b9c7cp+429,
    0x1.64271db8dd348p+483, 0x1.1f7341c9b9c7cp+429}, /* i=121 */
   {0x1.5296b75b2fbb2p+8, 0x1.652040b1f5bd3p+487, 0x1.2ee5cbf671144p+433,
    0x1.652040b1f5bd3p+487, 0x1.2ee5cbf671144p+433}, /* i=122 */
   {0x1.555d32906494dp+8, 0x1.661a11f1d512p+491, -0x1.ec33d4fb80a8cp+435,
    0x1.661a11f1d512p+491, -0x1.ec33d4fb80a8cp+435}, /* i=123 */
   {0x1.5823adc5996e7p+8, 0x1.671491f264075p+495, 0x1.e79871d9d332ap+440,
    0x1.671491f264075p+495, 0x1.e79871d9d332ap+440}, /* i=124 */
   {0x1.5aea28face482p+8, 0x1.680fc12de1129p+499, 0x1.d1f5df929ffcp+439,
    0x1.680fc12de1129p+499, 0x1.d1f5df929ffcp+439}, /* i=125 */
   {0x1.5db0a4300321dp+8, 0x1.690ba01edfe8ep+503, 0x1.f955228df51p+445,
    0x1.690ba01edfe8ep+503, 0x1.f955228df51p+445}, /* i=126 */
   {0x1.60771f6537fb7p+8, 0x1.6a082f4049fe1p+507, -0x1.ae13cb4465962p+452,
    0x1.6a082f4049fe1p+507, -0x1.ae13cb4465962p+452}, /* i=127 */
   {0x1.633d9a9a6cd52p+8, 0x1.6b056f0d5f048p+511, 0x1.0ef32cac389c4p+456,
    0x1.6b056f0d5f048p+511, 0x1.0ef32cac389c4p+456}, /* i=128 */
   {0x1.660415cfa1aedp+8, 0x1.6c036001b4a1ep+515, 0x1.7225365eeae5cp+461,
    0x1.6c036001b4a1ep+515, 0x1.7225365eeae5cp+461}, /* i=129 */
   {0x1.68ca9104d6887p+8, 0x1.6d02029936eeep+519, 0x1.5c69825df0628p+465,
    0x1.6d02029936eeep+519, 0x1.5c69825df0628p+465}, /* i=130 */
   {0x1.6b910c3a0b622p+8, 0x1.6e01575028f85p+523, 0x1.f3d253d1d2dep+466,
    0x1.6e01575028f85p+523, 0x1.f3d253d1d2dep+466}, /* i=131 */
   {0x1.6e57876f403bdp+8, 0x1.6f015ea324731p+527, -0x1.7cbd69647e6a8p+472,
    0x1.6f015ea324731p+527, -0x1.7cbd69647e6a8p+472}, /* i=132 */
   {0x1.711e02a475157p+8, 0x1.7002190f1a3cfp+531, 0x1.9f1e129357698p+475,
    0x1.7002190f1a3cfp+531, 0x1.9f1e129357698p+475}, /* i=133 */
   {0x1.73e47dd9a9ef2p+8, 0x1.7103871152defp+535, 0x1.7c532bca0695ap+480,
    0x1.7103871152defp+535, 0x1.7c532bca0695ap+480}, /* i=134 */
   {0x1.76aaf90edec8cp+8, 0x1.7205a9276e295p+539, 0x1.e15e16aeb9e12p+484,
    0x1.7205a9276e295p+539, 0x1.e15e16aeb9e12p+484}, /* i=135 */
   {0x1.7971744413a27p+8, 0x1.73087fcf64293p+543, 0x1.2e3ebdafe96e3p+489,
    0x1.73087fcf64293p+543, 0x1.2e3ebdafe96e3p+489}, /* i=136 */
   {0x1.7c37ef79487c2p+8, 0x1.740c0b8784c49p+547, -0x1.0db49c11b3f2ep+493,
    0x1.740c0b8784c49p+547, -0x1.0db49c11b3f2ep+493}, /* i=137 */
   {0x1.7efe6aae7d55cp+8, 0x1.75104cce783ccp+551, 0x1.f982b42573c1p+497,
    0x1.75104cce783ccp+551, 0x1.f982b42573c1p+497}, /* i=138 */
   {0x1.81c4e5e3b22f7p+8, 0x1.761544233fb2cp+555, 0x1.c7e237def3f5p+500,
    0x1.761544233fb2cp+555, 0x1.c7e237def3f5p+500}, /* i=139 */
   {0x1.848b6118e7092p+8, 0x1.771af20534d91p+559, -0x1.8b01db6348caep+505,
    0x1.771af20534d91p+559, -0x1.8b01db6348caep+505}, /* i=140 */
   {0x1.8751dc4e1be2cp+8, 0x1.782156f40a779p+563, -0x1.7bc9998c7d7ccp+508,
    0x1.782156f40a779p+563, -0x1.7bc9998c7d7ccp+508}, /* i=141 */
   {0x1.8a18578350bc7p+8, 0x1.7928736fccf0bp+567, -0x1.040f6b4c4e8c8p+511,
    0x1.7928736fccf0bp+567, -0x1.040f6b4c4e8c8p+511}, /* i=142 */
   {0x1.8cded2b885962p+8, 0x1.7a3047f8e1f2ep+571, 0x1.325fc95b46974p+516,
    0x1.7a3047f8e1f2ep+571, 0x1.325fc95b46974p+516}, /* i=143 */
   {0x1.8fa54dedba6fcp+8, 0x1.7b38d51008fd7p+575, 0x1.b16dea40944a9p+521,
    0x1.7b38d51008fd7p+575, 0x1.b16dea40944a9p+521}, /* i=144 */
   {0x1.926bc922ef497p+8, 0x1.7c421b365be6cp+579, -0x1.dc0a97fa26ae2p+525,
    0x1.7c421b365be6cp+579, -0x1.dc0a97fa26ae2p+525}, /* i=145 */
   {0x1.9532445824232p+8, 0x1.7d4c1aed4e8cfp+583, -0x1.2c95403a45576p+528,
    0x1.7d4c1aed4e8cfp+583, -0x1.2c95403a45576p+528}, /* i=146 */
   {0x1.97f8bf8d58fccp+8, 0x1.7e56d4b6af5c4p+587, -0x1.d1f7a35f34323p+533,
    0x1.7e56d4b6af5c4p+587, -0x1.d1f7a35f34323p+533}, /* i=147 */
   {0x1.9abf3ac28dd67p+8, 0x1.7f624914a7d5dp+591, 0x1.ccfe6c3a5a9cp+536,
    0x1.7f624914a7d5dp+591, 0x1.ccfe6c3a5a9cp+536}, /* i=148 */
   {0x1.9d85b5f7c2b01p+8, 0x1.806e7889bc286p+595, -0x1.5c2cfdbcfdaap+539,
    0x1.806e7889bc286p+595, -0x1.5c2cfdbcfdaap+539}, /* i=149 */
   {0x1.a04c312cf789cp+8, 0x1.817b6398cc2fp+599, 0x1.224ca8a9fb47p+541,
    0x1.817b6398cc2fp+599, 0x1.224ca8a9fb47p+541}, /* i=150 */
   {0x1.a312ac622c637p+8, 0x1.82890ac513097p+603, 0x1.be87e5468a99ep+549,
    0x1.82890ac513097p+603, 0x1.be87e5468a99ep+549}, /* i=151 */
   {0x1.a5d92797613d1p+8, 0x1.83976e9227a3dp+607, 0x1.ad73cbfd854d2p+552,
    0x1.83976e9227a3dp+607, 0x1.ad73cbfd854d2p+552}, /* i=152 */
   {0x1.a89fa2cc9616cp+8, 0x1.84a68f83fd3f9p+611, -0x1.d1e43ed15602p+555,
    0x1.84a68f83fd3f9p+611, -0x1.d1e43ed15602p+555}, /* i=153 */
   {0x1.ab661e01caf07p+8, 0x1.85b66e1ee322cp+615, 0x1.9bd201725ed39p+561,
    0x1.85b66e1ee322cp+615, 0x1.9bd201725ed39p+561}, /* i=154 */
   {0x1.ae2c9936ffca1p+8, 0x1.86c70ae785212p+619, -0x1.b7df86a3901eep+565,
    0x1.86c70ae785212p+619, -0x1.b7df86a3901eep+565}, /* i=155 */
   {0x1.b0f3146c34a3cp+8, 0x1.87d86662ec25dp+623, 0x1.ecbcd81926fb3p+569,
    0x1.87d86662ec25dp+623, 0x1.ecbcd81926fb3p+569}, /* i=156 */
   {0x1.b3b98fa1697d7p+8, 0x1.88ea81167de2ap+627, -0x1.e8df2f283e5fp+572,
    0x1.88ea81167de2ap+627, -0x1.e8df2f283e5fp+572}, /* i=157 */
   {0x1.b6800ad69e571p+8, 0x1.89fd5b87fd594p+631, 0x1.e8a5a1c1512bp+574,
    0x1.89fd5b87fd594p+631, 0x1.e8a5a1c1512bp+574}, /* i=158 */
   {0x1.b946860bd330cp+8, 0x1.8b10f63d8b674p+635, -0x1.7fb230d7d35dbp+581,
    0x1.8b10f63d8b674p+635, -0x1.7fb230d7d35dbp+581}, /* i=159 */
   {0x1.bc0d0141080a6p+8, 0x1.8c2551bda65acp+639, -0x1.26165c2ff2f8bp+585,
    0x1.8c2551bda65acp+639, -0x1.26165c2ff2f8bp+585}, /* i=160 */
   {0x1.bed37c763ce41p+8, 0x1.8d3a6e8f2af9ap+643, -0x1.dc464e4ec6448p+586,
    0x1.8d3a6e8f2af9ap+643, -0x1.dc464e4ec6448p+586}, /* i=161 */
   {0x1.c199f7ab71bdcp+8, 0x1.8e504d3954165p+647, -0x1.058bf4712632dp+593,
    0x1.8e504d3954165p+647, -0x1.058bf4712632dp+593}, /* i=162 */
   {0x1.c46072e0a6976p+8, 0x1.8f66ee43bb1b9p+651, 0x1.6851e7e4b989p+596,
    0x1.8f66ee43bb1b9p+651, 0x1.6851e7e4b989p+596}, /* i=163 */
   {0x1.c726ee15db711p+8, 0x1.907e52365899fp+655, 0x1.ed5f2ddce5a5p+599,
    0x1.907e52365899fp+655, 0x1.ed5f2ddce5a5p+599}, /* i=164 */
   {0x1.c9ed694b104acp+8, 0x1.9196799983f46p+659, -0x1.99926389b844cp+604,
    0x1.9196799983f46p+659, -0x1.99926389b844cp+604}, /* i=165 */
   {0x1.ccb3e48045246p+8, 0x1.92af64f5f3ed4p+663, -0x1.58ed30816683p+605,
    0x1.92af64f5f3ed4p+663, -0x1.58ed30816683p+605}, /* i=166 */
   {0x1.cf7a5fb579fe1p+8, 0x1.93c914d4bf34fp+667, -0x1.0083aa488e25bp+613,
    0x1.93c914d4bf34fp+667, -0x1.0083aa488e25bp+613}, /* i=167 */
   {0x1.d240daeaaed7cp+8, 0x1.94e389bf5c15bp+671, 0x1.8474de356966cp+615,
    0x1.94e389bf5c15bp+671, 0x1.8474de356966cp+615}, /* i=168 */
   {0x1.d507561fe3b16p+8, 0x1.95fec43fa1021p+675, -0x1.709e63bd4ab62p+621,
    0x1.95fec43fa1021p+675, -0x1.709e63bd4ab62p+621}, /* i=169 */
   {0x1.d7cdd155188b1p+8, 0x1.971ac4dfc5243p+679, 0x1.59fd4818a390ep+624,
    0x1.971ac4dfc5243p+679, 0x1.59fd4818a390ep+624}, /* i=170 */
   {0x1.da944c8a4d64cp+8, 0x1.98378c2a60099p+683, 0x1.38801e5eca72ep+628,
    0x1.98378c2a60099p+683, 0x1.38801e5eca72ep+628}, /* i=171 */
   {0x1.dd5ac7bf823e6p+8, 0x1.99551aaa6a321p+687, 0x1.962a4e85c1ea6p+632,
    0x1.99551aaa6a321p+687, 0x1.962a4e85c1ea6p+632}, /* i=172 */
   {0x1.e02142f4b7181p+8, 0x1.9a7370eb3da0ep+691, 0x1.cf232c5aadea3p+637,
    0x1.9a7370eb3da0ep+691, 0x1.cf232c5aadea3p+637}, /* i=173 */
   {0x1.e2e7be29ebf1bp+8, 0x1.9b928f78956d8p+695, 0x1.b7fefcbe883d1p+641,
    0x1.9b928f78956d8p+695, 0x1.b7fefcbe883d1p+641}, /* i=174 */
   {0x1.e5ae395f20cb6p+8, 0x1.9cb276de8ed4dp+699, 0x1.6dcc83a185adep+644,
    0x1.9cb276de8ed4dp+699, 0x1.6dcc83a185adep+644}, /* i=175 */
   {0x1.e874b49455a51p+8, 0x1.9dd327a9a8c9dp+703, -0x1.ee82b70c933p+641,
    0x1.9dd327a9a8c9dp+703, -0x1.ee82b70c933p+641}, /* i=176 */
   {0x1.eb3b2fc98a7ebp+8, 0x1.9ef4a266c486fp+707, 0x1.92efd6c195825p+653,
    0x1.9ef4a266c486fp+707, 0x1.92efd6c195825p+653}, /* i=177 */
   {0x1.ee01aafebf586p+8, 0x1.a016e7a32620fp+711, 0x1.334aca1c0a24p+657,
    0x1.a016e7a32620fp+711, 0x1.334aca1c0a24p+657}, /* i=178 */
   {0x1.f0c82633f4321p+8, 0x1.a139f7ec74304p+715, -0x1.0a2fa13d80024p+661,
    0x1.a139f7ec74304p+715, -0x1.0a2fa13d80024p+661}, /* i=179 */
   {0x1.f38ea169290bbp+8, 0x1.a25dd3d0b8638p+719, -0x1.c822f0b528b68p+663,
    0x1.a25dd3d0b8638p+719, -0x1.c822f0b528b68p+663}, /* i=180 */
   {0x1.f6551c9e5de56p+8, 0x1.a3827bde6013bp+723, -0x1.995a210ca10f1p+669,
    0x1.a3827bde6013bp+723, -0x1.995a210ca10f1p+669}, /* i=181 */
   {0x1.f91b97d392bf1p+8, 0x1.a4a7f0a43beccp+727, -0x1.392a89837b612p+673,
    0x1.a4a7f0a43beccp+727, -0x1.392a89837b612p+673}, /* i=182 */
   {0x1.fbe21308c798bp+8, 0x1.a5ce32b180817p+731, 0x1.27fb21e22edb9p+677,
    0x1.a5ce32b180817p+731, 0x1.27fb21e22edb9p+677}, /* i=183 */
   {0x1.fea88e3dfc726p+8, 0x1.a6f54295c6e09p+735, -0x1.aeaa3b01c0afep+680,
    0x1.a6f54295c6e09p+735, -0x1.aeaa3b01c0afep+680}, /* i=184 */
   {0x1.00b784b998a6p+9, 0x1.a81d20e10c225p+739, -0x1.9ee61527f1503p+685,
    0x1.a81d20e10c225p+739, -0x1.9ee61527f1503p+685}, /* i=185 */
   {0x1.021ac2543312ep+9, 0x1.a945ce23b29c9p+743, -0x1.33fa72fb6b984p+689,
    0x1.a945ce23b29c9p+743, -0x1.33fa72fb6b984p+689}, /* i=186 */
   {0x1.037dffeecd7fbp+9, 0x1.aa6f4aee80eb4p+747, -0x1.983a128a2a5cp+692,
    0x1.aa6f4aee80eb4p+747, -0x1.983a128a2a5cp+692}, /* i=187 */
   {0x1.04e13d8967ec8p+9, 0x1.ab9997d2a38fep+751, -0x1.627d51485f464p+696,
    0x1.ab9997d2a38fep+751, -0x1.627d51485f464p+696}, /* i=188 */
   {0x1.06447b2402596p+9, 0x1.acc4b561ac99cp+755, -0x1.f4ecb3a8c437cp+700,
    0x1.acc4b561ac99cp+755, -0x1.f4ecb3a8c437cp+700}, /* i=189 */
   {0x1.07a7b8be9cc63p+9, 0x1.adf0a42d934bdp+759, -0x1.71e8135a409c4p+705,
    0x1.adf0a42d934bdp+759, -0x1.71e8135a409c4p+705}, /* i=190 */
   {0x1.090af6593733p+9, 0x1.af1d64c8b5a5p+763, -0x1.f5ff6e969e6e6p+709,
    0x1.af1d64c8b5a5p+763, -0x1.f5ff6e969e6e6p+709}, /* i=191 */
   {0x1.0a6e33f3d19fep+9, 0x1.b04af7c5d807cp+767, -0x1.5d7b17e7b0921p+713,
    0x1.b04af7c5d807cp+767, -0x1.5d7b17e7b0921p+713}, /* i=192 */
   {0x1.0bd1718e6c0cbp+9, 0x1.b1795db824df2p+771, -0x1.9e51195a5e74cp+716,
    0x1.b1795db824df2p+771, -0x1.9e51195a5e74cp+716}, /* i=193 */
   {0x1.0d34af2906798p+9, 0x1.b2a897332e2a6p+775, -0x1.b89ec0e5c5ed3p+721,
    0x1.b2a897332e2a6p+775, -0x1.b89ec0e5c5ed3p+721}, /* i=194 */
   {0x1.0e97ecc3a0e65p+9, 0x1.b3d8a4caeced2p+779, 0x1.fc1d7db3c1163p+725,
    0x1.b3d8a4caeced2p+779, 0x1.fc1d7db3c1163p+725}, /* i=195 */
   {0x1.0ffb2a5e3b533p+9, 0x1.b5098713c1e48p+783, -0x1.6b0e11c1c6958p+729,
    0x1.b5098713c1e48p+783, -0x1.6b0e11c1c6958p+729}, /* i=196 */
   {0x1.115e67f8d5cp+9, 0x1.b63b3ea274f48p+787, 0x1.f4a32a36efa16p+732,
    0x1.b63b3ea274f48p+787, 0x1.f4a32a36efa16p+732}, /* i=197 */
   {0x1.12c1a593702cdp+9, 0x1.b76dcc0c36b8bp+791, -0x1.6288b93466f0bp+737,
    0x1.b76dcc0c36b8bp+791, -0x1.6288b93466f0bp+737}, /* i=198 */
   {0x1.1424e32e0a99bp+9, 0x1.b8a12fe6a0295p+795, 0x1.5ecc5d89d82c4p+741,
    0x1.b8a12fe6a0295p+795, 0x1.5ecc5d89d82c4p+741}, /* i=199 */
   {0x1.158820c8a5068p+9, 0x1.b9d56ac7b23f4p+799, -0x1.e38ab613c108p+738,
    0x1.b9d56ac7b23f4p+799, -0x1.e38ab613c108p+738}, /* i=200 */
   {0x1.16eb5e633f735p+9, 0x1.bb0a7d45d786bp+803, 0x1.73c1cf5b4ff48p+747,
    0x1.bb0a7d45d786bp+803, 0x1.73c1cf5b4ff48p+747}, /* i=201 */
   {0x1.184e9bfdd9e03p+9, 0x1.bc4067f7e3c49p+807, 0x1.65b5e76763acp+747,
    0x1.bc4067f7e3c49p+807, 0x1.65b5e76763acp+747}, /* i=202 */
   {0x1.19b1d998744dp+9, 0x1.bd772b751398fp+811, -0x1.da1d2fd5b03ep+754,
    0x1.bd772b751398fp+811, -0x1.da1d2fd5b03ep+754}, /* i=203 */
   {0x1.1b1517330eb9dp+9, 0x1.beaec8550e15ap+815, -0x1.e106425e362fdp+761,
    0x1.beaec8550e15ap+815, -0x1.e106425e362fdp+761}, /* i=204 */
   {0x1.1c7854cda926bp+9, 0x1.bfe73f2fe4626p+819, -0x1.30417c6c605ap+761,
    0x1.bfe73f2fe4626p+819, -0x1.30417c6c605ap+761}, /* i=205 */
   {0x1.1ddb926843938p+9, 0x1.c120909e115efp+823, -0x1.206dc5b067a48p+769,
    0x1.c120909e115efp+823, -0x1.206dc5b067a48p+769}, /* i=206 */
   {0x1.1f3ed002de005p+9, 0x1.c25abd387b3c9p+827, 0x1.ca7ad33c54962p+773,
    0x1.c25abd387b3c9p+827, 0x1.ca7ad33c54962p+773}, /* i=207 */
   {0x1.20a20d9d786d3p+9, 0x1.c395c5987322p+831, 0x1.07ab711d35afcp+777,
    0x1.c395c5987322p+831, 0x1.07ab711d35afcp+777}, /* i=208 */
   {0x1.22054b3812dap+9, 0x1.c4d1aa57b4cc2p+835, -0x1.28feccc1671f2p+780,
    0x1.c4d1aa57b4cc2p+835, -0x1.28feccc1671f2p+780}, /* i=209 */
   {0x1.236888d2ad46dp+9, 0x1.c60e6c10682b6p+839, -0x1.cfb9d52f2be79p+785,
    0x1.c60e6c10682b6p+839, -0x1.cfb9d52f2be79p+785}, /* i=210 */
   {0x1.24cbc66d47b3bp+9, 0x1.c74c0b5d21068p+843, 0x1.72803b10e5d1p+788,
    0x1.c74c0b5d21068p+843, 0x1.72803b10e5d1p+788}, /* i=211 */
   {0x1.262f0407e2208p+9, 0x1.c88a88d8de9bp+847, -0x1.c5429a1197244p+791,
    0x1.c88a88d8de9bp+847, -0x1.c5429a1197244p+791}, /* i=212 */
   {0x1.279241a27c8d5p+9, 0x1.c9c9e51f0d3d6p+851, 0x1.afd5bd851dd63p+797,
    0x1.c9c9e51f0d3d6p+851, 0x1.afd5bd851dd63p+797}, /* i=213 */
   {0x1.28f57f3d16fa3p+9, 0x1.cb0a20cb85fbap+855, -0x1.612849dbcd002p+800,
    0x1.cb0a20cb85fbap+855, -0x1.612849dbcd002p+800}, /* i=214 */
   {0x1.2a58bcd7b167p+9, 0x1.cc4b3c7a8e3c4p+859, -0x1.91ab278ea296fp+805,
    0x1.cc4b3c7a8e3c4p+859, -0x1.91ab278ea296fp+805}, /* i=215 */
   {0x1.2bbbfa724bd3dp+9, 0x1.cd8d38c8d962ap+863, -0x1.44d024e25be64p+808,
    0x1.cd8d38c8d962ap+863, -0x1.44d024e25be64p+808}, /* i=216 */
   {0x1.2d1f380ce640bp+9, 0x1.ced0165388704p+867, -0x1.f2e7231c23f15p+813,
    0x1.ced0165388704p+867, -0x1.f2e7231c23f15p+813}, /* i=217 */
   {0x1.2e8275a780ad8p+9, 0x1.d013d5b829a33p+871, -0x1.2847f05c852bp+814,
    0x1.d013d5b829a33p+871, -0x1.2847f05c852bp+814}, /* i=218 */
   {0x1.2fe5b3421b1a5p+9, 0x1.d1587794ba1dep+875, 0x1.12b8eb0337defp+821,
    0x1.d1587794ba1dep+875, 0x1.12b8eb0337defp+821}, /* i=219 */
   {0x1.3148f0dcb5872p+9, 0x1.d29dfc87a54d2p+879, -0x1.134000aa23a2p+824,
    0x1.d29dfc87a54d2p+879, -0x1.134000aa23a2p+824}, /* i=220 */
   {0x1.32ac2e774ff4p+9, 0x1.d3e4652fc5a98p+883, 0x1.8704f4c0fdbe6p+829,
    0x1.d3e4652fc5a98p+883, 0x1.8704f4c0fdbe6p+829}, /* i=221 */
   {0x1.340f6c11ea60dp+9, 0x1.d52bb22c641b6p+887, 0x1.73525d19565e8p+833,
    0x1.d52bb22c641b6p+887, 0x1.73525d19565e8p+833}, /* i=222 */
   {0x1.3572a9ac84cdap+9, 0x1.d673e41d39a6dp+891, 0x1.1a73f52ffafep+833,
    0x1.d673e41d39a6dp+891, 0x1.1a73f52ffafep+833}, /* i=223 */
   {0x1.36d5e7471f3a8p+9, 0x1.d7bcfba26f0b1p+895, 0x1.437205f520b78p+841,
    0x1.d7bcfba26f0b1p+895, 0x1.437205f520b78p+841}, /* i=224 */
   {0x1.383924e1b9a75p+9, 0x1.d906f95c9c5f8p+899, 0x1.9e43b7362ec71p+845,
    0x1.d906f95c9c5f8p+899, 0x1.9e43b7362ec71p+845}, /* i=225 */
   {0x1.399c627c54142p+9, 0x1.da51ddeccac33p+903, -0x1.59d6da8372b0ep+849,
    0x1.da51ddeccac33p+903, -0x1.59d6da8372b0ep+849}, /* i=226 */
   {0x1.3affa016ee81p+9, 0x1.db9da9f473fb6p+907, 0x1.94f24f867cbccp+851,
    0x1.db9da9f473fb6p+907, 0x1.94f24f867cbccp+851}, /* i=227 */
   {0x1.3c62ddb188eddp+9, 0x1.dcea5e15820fep+911, -0x1.1e53b239af6fp+857,
    0x1.dcea5e15820fep+911, -0x1.1e53b239af6fp+857}, /* i=228 */
   {0x1.3dc61b4c235aap+9, 0x1.de37faf250fddp+915, 0x1.7249852572c69p+861,
    0x1.de37faf250fddp+915, 0x1.7249852572c69p+861}, /* i=229 */
   {0x1.3f2958e6bdc78p+9, 0x1.df86812dae55fp+919, 0x1.1f6cbdc6ac764p+865,
    0x1.df86812dae55fp+919, 0x1.1f6cbdc6ac764p+865}, /* i=230 */
   {0x1.408c968158345p+9, 0x1.e0d5f16ad8d76p+923, 0x1.7817d73668aeep+869,
    0x1.e0d5f16ad8d76p+923, 0x1.7817d73668aeep+869}, /* i=231 */
   {0x1.41efd41bf2a12p+9, 0x1.e2264c4d8226bp+927, 0x1.bd3e35a99d614p+872,
    0x1.e2264c4d8226bp+927, 0x1.bd3e35a99d614p+872}, /* i=232 */
   {0x1.435311b68d0ep+9, 0x1.e3779279ce6adp+931, -0x1.71c8dc115e8bep+876,
    0x1.e3779279ce6adp+931, -0x1.71c8dc115e8bep+876}, /* i=233 */
   {0x1.44b64f51277adp+9, 0x1.e4c9c49453e75p+935, 0x1.fce96ff4df11ap+881,
    0x1.e4c9c49453e75p+935, 0x1.fce96ff4df11ap+881}, /* i=234 */
   {0x1.46198cebc1e7ap+9, 0x1.e61ce3421cb72p+939, -0x1.36c28c6f65faap+885,
    0x1.e61ce3421cb72p+939, -0x1.36c28c6f65faap+885}, /* i=235 */
   {0x1.477cca865c548p+9, 0x1.e770ef28a6684p+943, -0x1.5d827b1480342p+889,
    0x1.e770ef28a6684p+943, -0x1.5d827b1480342p+889}, /* i=236 */
   {0x1.48e00820f6c15p+9, 0x1.e8c5e8ede195bp+947, -0x1.b193bea6cbc4p+888,
    0x1.e8c5e8ede195bp+947, -0x1.b193bea6cbc4p+888}, /* i=237 */
   {0x1.4a4345bb912e2p+9, 0x1.ea1bd13833a56p+951, -0x1.8369a2e682074p+897,
    0x1.ea1bd13833a56p+951, -0x1.8369a2e682074p+897}, /* i=238 */
   {0x1.4ba683562b9bp+9, 0x1.eb72a8ae76636p+955, 0x1.90bdfa7851a2dp+901,
    0x1.eb72a8ae76636p+955, 0x1.90bdfa7851a2dp+901}, /* i=239 */
   {0x1.4d09c0f0c607dp+9, 0x1.ecca6ff7f79afp+959, 0x1.85d731cf20318p+903,
    0x1.ecca6ff7f79afp+959, 0x1.85d731cf20318p+903}, /* i=240 */
   {0x1.4e6cfe8b6074ap+9, 0x1.ee2327bc7ad7cp+963, 0x1.2783a6814b4e2p+908,
    0x1.ee2327bc7ad7cp+963, 0x1.2783a6814b4e2p+908}, /* i=241 */
   {0x1.4fd03c25fae18p+9, 0x1.ef7cd0a43900ep+967, -0x1.e51ec8f6cdb24p+911,
    0x1.ef7cd0a43900ep+967, -0x1.e51ec8f6cdb24p+911}, /* i=242 */
   {0x1.513379c0954e5p+9, 0x1.f0d76b57dff05p+971, 0x1.21db342298f3cp+915,
    0x1.f0d76b57dff05p+971, 0x1.21db342298f3cp+915}, /* i=243 */
   {0x1.5296b75b2fbb2p+9, 0x1.f232f8809438ap+975, -0x1.42d1e29539c42p+921,
    0x1.f232f8809438ap+975, -0x1.42d1e29539c42p+921}, /* i=244 */
   {0x1.53f9f4f5ca27fp+9, 0x1.f38f78c7f08p+979, -0x1.defc5b6b1d818p+924,
    0x1.f38f78c7f08p+979, -0x1.defc5b6b1d818p+924}, /* i=245 */
   {0x1.555d32906494dp+9, 0x1.f4ececd8064fap+983, -0x1.7ce1fe211d10ep+928,
    0x1.f4ececd8064fap+983, -0x1.7ce1fe211d10ep+928}, /* i=246 */
   {0x1.56c0702aff01ap+9, 0x1.f64b555b5d6c4p+987, -0x1.5162e4da751p+926,
    0x1.f64b555b5d6c4p+987, -0x1.5162e4da751p+926}, /* i=247 */
   {0x1.5823adc5996e7p+9, 0x1.f7aab2fcf5a0bp+991, -0x1.be8ea835f153bp+937,
    0x1.f7aab2fcf5a0bp+991, -0x1.be8ea835f153bp+937}, /* i=248 */
   {0x1.5986eb6033db5p+9, 0x1.f90b066846564p+995, -0x1.123fde796711cp+940,
    0x1.f90b066846564p+995, -0x1.123fde796711cp+940}, /* i=249 */
   {0x1.5aea28face482p+9, 0x1.fa6c50493e2adp+999, -0x1.3a942a74008bcp+945,
    0x1.fa6c50493e2adp+999, -0x1.3a942a74008bcp+945}, /* i=250 */
   {0x1.5c4d669568b4fp+9, 0x1.fbce914c44bebp+1003, -0x1.14ad532a31838p+946,
    0x1.fbce914c44bebp+1003, -0x1.14ad532a31838p+946}, /* i=251 */
   {0x1.5db0a4300321dp+9, 0x1.fd31ca1e3a4cap+1007, 0x1.61f736223ff12p+952,
    0x1.fd31ca1e3a4cap+1007, 0x1.61f736223ff12p+952}, /* i=252 */
   {0x1.5f13e1ca9d8eap+9, 0x1.fe95fb6c773eap+1011, -0x1.1b541ab2458e1p+957,
    0x1.fe95fb6c773eap+1011, -0x1.1b541ab2458e1p+957}, /* i=253 */
   {0x1.60771f6537fb7p+9, 0x1.fffb25e4cdffcp+1015, 0x1.bc21d5499b93fp+961,
    0x1.fffb25e4cdffcp+1015, 0x1.bc21d5499b93fp+961}, /* i=254 */
   {0x1.61da5cffd2685p+9, 0x1.00b0a51ac549cp+1020, -0x1.807ae217f4df6p+966,
    0x1.00b0a51ac549cp+1020, -0x1.807ae217f4df6p+966}, /* i=255 */
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
   For each j, 0 <= j < 256, let {xj, sj, cj} be the U[j] values,
   such that sj and cj approximate sinh(xj) and cosh(xj) with a few
   extra bits, then sj + Ul[j][0] and cj + Ul[j][1] approximate
   sinh(xj) and cosh(xj) with at least 107 bits.
   Generated with build_table_Ul(U0,U1,U2) from the file sinh.sage,
   where U0,U1,U2 contain U[i][0],U[i][1],U[i][2] respectively. */
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
   Assumes |a| >= |b| or a=0.  */
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
     A bound on the error is given in Theorem 1 from [2]:
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

// Algorithm 2Sum from [3]
static inline void
two_sum (double *s, double *t, double a, double b)
{
  *s = a + b;
  double ap = *s - b;
  double bp = *s - ap;
  double da = a - ap;
  double db = b - bp;
  *t = da + db;
}

// Algorithm DWPlusFP from [3]
static inline void
fast_sum_acc (double *hi, double *lo, double a, double bh, double bl)
{
  double sh, sl;
  two_sum (&sh, &sl, bh, a);
  double v = bl + sl;
  fast_two_sum (hi, lo, sh, v);
}

// Add (ah + al) + (bh + bl), assuming |ah| >= |bh|
static inline void fast_sum2(double *hi, double *lo, double ah, double al,
                             double bh, double bl) {
  fast_two_sum (hi, lo, ah, bh);
  *lo += al + bl;
}

// Algorithm SloppyDWPlusDW from [3]
static inline void
fast_sum2_acc1 (double *hi, double *lo, double ah, double al,
                double bh, double bl)
{
  double sh, sl;
  two_sum (&sh, &sl, ah, bh);
  double v = al + bl;
  double w = sl + v;
  fast_two_sum (hi, lo, sh, w);
}

// Algorithm AccurateDWPlusDW from [3]
static inline void
fast_sum2_acc2 (double *hi, double *lo, double ah, double al,
                double bh, double bl)
{
  double sh, sl;
  two_sum (&sh, &sl, ah, bh);
  double th, tl;
  two_sum (&th, &tl, al, bl);
  double c = sl + th;
  double vh, vl;
  fast_two_sum (&vh, &vl, sh, c);
  double w = tl + vl;
  fast_two_sum (hi, lo, vh, w);
}

// Multiply exactly a and b, such that *hi + *lo = a * b.
// This is called 2Prod in [3].
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

#if 0
// Multiply a double with a double double : a * (bh + bl)
// using algorithm DWTimesFP1 from [3].
static inline void
s_mul_acc1 (double *hi, double *lo, double a, double bh, double bl)
{
  // hi is ch from [3], a is y, bh is xh, bl is xl
  double cl1;
  a_mul (hi, &cl1, a, bh); /* exact */
  double cl2 = a * bl;
  fast_two_sum (hi, lo, *hi, cl2);
  double tl2 = *lo + cl1;
  fast_two_sum (hi, lo, *hi, tl2);
}

// Multiply a double with a double double : a * (bh + bl)
// using algorithm DWTimesFP2 from [3].
static inline void
s_mul_acc2 (double *hi, double *lo, double a, double bh, double bl)
{
  // hi is ch from [3], a is y, bh is xh, bl is xl
  double cl1;
  a_mul (hi, &cl1, a, bh); /* exact */
  double cl2 = a * bl;
  double cl3 = cl1 + cl2;
  fast_two_sum (hi, lo, *hi, cl3);
}
#endif

// Multiply a double with a double double : a * (bh + bl)
// using algorithm DWTimesFP3 from [3].
// This is essentially the same as s_mul(), with a post-normalization.
static inline void
s_mul_acc3 (double *hi, double *lo, double a, double bh, double bl)
{
  // hi is ch from [3], a is y, bh is xh, bl is xl
  double cl1;
  a_mul (hi, &cl1, a, bh); /* exact */
  double cl3 = __builtin_fma (a, bl, cl1);
  fast_two_sum (hi, lo, *hi, cl3);
}

/* Put in hi+lo an approximation of (ah+al)*(bh+bl), with relative error
   bounded by 2^-101.41, assuming |al| < ulp(ah) and |bl| < ulp(bh). */
static inline void
d_mul (double *hi, double *lo, double ah, double al, double bh, double bl)
{
  double s, t;

  a_mul (hi, &s, ah, bh); /* exact */
  /* For the error analysis, we assume |al| < ulp(ah) and |bl| < ulp(bh).
     Also, to simplify the analysis, we assume 1 <= ah, bh < 2,
     thus 1 <= hi < 4, and |al|, |bl| < 2^-52, and |s| < ulp(hi) <= 2^-51. */
  t = __builtin_fma (al, bh, s);
  /* |t| < |al|*bh+|s| <= 2^-52*2+2^-51 = 2^-50, and the rounding error is
     bounded by ulp(t-eps) = 2^-103. */
  *lo = __builtin_fma (ah, bl, t);
  /* |lo| < ah*|bl|+|t| <= 2*2^-52+2^-50 < 2^-49.41, and the rounding error is
     bounded by ulp(2^-49.41) = 2^-102. */
  /* The total rounding error is bounded by 2^-103+2^-102 <= 2^-101.41,
     thus also 2^-101.41 relatively to hi+lo (and this relative bound
     holds whatever the binades of ah+al and bh+bl). */
}

#if 0
// using algorithm DWTimesDW1 from [3], with relative error bounded by 7u^2
// (for rounding to nearest)
static inline void
d_mul_acc1 (double *hi, double *lo, double ah, double al, double bh, double bl)
{
  double cl1;
  a_mul (hi, &cl1, ah, bh);
  double tl1 = ah * bl;
  double tl2 = al * bh;
  double cl2 = tl1 + tl2;
  double cl3 = cl1 + cl2;
  fast_two_sum (hi, lo, *hi, cl3);
}
#endif

// using algorithm DWTimesDW2 from [3], with relative error bounded by 6u^2
// (for rounding to nearest)
static inline void
d_mul_acc2 (double *hi, double *lo, double ah, double al, double bh, double bl)
{
  double cl1;
  a_mul (hi, &cl1, ah, bh);
  double tl = ah * bl;
  double cl2 = __builtin_fma (al, bh, tl);
  double cl3 = cl1 + cl2;
  fast_two_sum (hi, lo, *hi, cl3);
}

// using algorithm DWTimesDW3 from [3], with relative error bounded by 5u^2
// (for rounding to nearest)
static inline void
d_mul_acc3 (double *hi, double *lo, double ah, double al, double bh, double bl)
{
  double cl1;
  a_mul (hi, &cl1, ah, bh);
  double tl0 = al * bl;
  double tl1 = __builtin_fma (ah, bl, tl0);
  double cl2 = __builtin_fma (al, bh, tl1);
  double cl3 = cl1 + cl2;
  fast_two_sum (hi, lo, *hi, cl3);
}

/* Put in hi+lo an approximation of 1/(bh+bl),
   with |lo| < 4 ulp(hi), and relative error < 2^-102.41
   if (bh,bl) is normalized. */
static inline void
d_inv (double *hi, double *lo, double bh, double bl)
{
  /* Warning: the error analysis below assumes |bl| < ulp(bh) at input.
     If this is not satisfied, call fast_two_sum() to normalize (bh,bl). */
  /* for the error analysis, without loss of generality,
     we assume 1/2 <= bh < 1, thus ulp(bh) = 2^-53 */
  *hi = 1.0 / bh;
  /* since bh is in the binade (1/2,1), hi is in the binade (1,2),
     where the ulp value is 2^-52, thus |hi - 1/bh| < 2^-52 */
  /* We use Newton's iteration: hi -> hi + hi*(1-hi*b).
     Lemma 3.7 from [1] says that if hi' = hi + hi*(1-hi*b):
     |hi' - 1/b| <= hi^2/theta^3*(1/b-hi)^2
     with min(1/b,hi) <= theta <= max(1/b,hi).
     (a) if hi <= 1/bh, then hi <= theta thus hi^2/theta^3 <= 1/theta.
     (b) if 1/bh <= hi, then 1/bh <= theta <= hi, thus since
         |hi - 1/bh| < 2^-52, theta >= hi - 2^-52 >= hi (1 - 2^-52) thus
         hi^2/theta^3 <= 1/theta * 1/(1 - 2^-52)^2.
     In both cases (a) and (b) we have since theta >= 1:
     hi^2/theta^3 <= 1/theta * 1/(1 - 2^-52)^2 <= 1/(1 - 2^-52)^2.
     If follows:
     |hi' - 1/b| <= 2^-104/(1 - 2^-52)^2.
  */
  double e;
  e = __builtin_fma (-*hi, bh, 1.0);
  /* Since |hi - 1/bh| < 2^-52, we have |hi*bh-1| < bh*2^-52 < 2^-52.
     Also, hi is an integer multiple of 2^-52 and bh an integer multiple of
     2^-53, thus hi*bh is an integer multiple of 2^-105, likewise for hi*bh-1.
     Thus hi*bh-1 = k*2^-105 and since |hi*bh-1| < 2^-52, |k| < 2^53, thus
     hi*bh-1 is exactly representable, and |e| < 2^-52. */
  e = __builtin_fma (-*hi, bl, e);
  /* Since |hi| <= 2 and |bl| < ulp(bh) = 2^-53, we have |hi*bl| < 2^-52,
     thus now |e| < 2^-51, and the rounding error is bounded by
     2^-105*|e| from [2,Theorem 2], thus by 2^-157.
     This rounding error is multiplied by hi below, with hi < 2, thus
     contributes to at most 2^-156. */
  *lo = *hi * e;
  /* |lo| < 2*2^-51 = 2^-50. Since ulp(hi) >= 2^-52, we have |lo| < 4 ulp(hi),
     and this holds whatever the initial binade of bh+bl.
     The rounding error in hi * e is bounded by ulp(2^-50-eps) = 2^-103. */

  /* The absolute difference between hi+lo and 1/(bh+bl) is bounded by:
   * 2^-104/(1 - 2^-52)^2 for the maximal error in Newton's iteration
   * 2^-156 for the rounding error in __builtin_fma (-hi, bl, e)
   * 2^-103 for the rounding error in hi * e
   This gives an absolute error bound of 2^-102.41 (and the same
   relative error bound since hi + lo >= 1).
  */
}

/* Put in hi+lo an approximation of (ah+al)/(bh+bl), with relative error
   bounded by 2^-100.82, assuming |al| < ulp(ah) and |bl| < ulp(bh).
   See also Algorithms 16 and 17 and Theorem 7.1 from [3], which gives a
   bound of 15u^2 + 56u^3 < 2^-102.09 (for rounding to nearest presumably).
   Algorithm 17 has 2 divisions, 2 multiplications, 12 additions/subtractions,
   1 fma, whereas our algorithm has 1 division, 5 fmas, 2 mul.
   FIXME: use Karp-Markstein's trick for the division.
*/
static inline void
d_div (double *hi, double *lo, double ah, double al, double bh, double bl)
{
  /* Warning: the error analysis below assumes |al| < ulp(ah) and
     |bl| < ulp(bh) at input. If this is not satisfied, call fast_two_sum()
     to normalize (ah,al) and (bh,bl). */
  d_inv (hi, lo, bh, bl);
  /* |hi + lo - 1/(bh + bl)| < 2^-102.41 * (hi + lo) */
  d_mul (hi, lo, ah, al, *hi, *lo);
  /* |hi + lo - (ah + al) * (hi_in + lo_in)| < 2^-101.41 * |hi + lo| thus:
     |hi + lo - (ah + al) / (bh + bl)|
     < |hi + lo - (ah + al) * (hi_in + lo_in)|
     + |ah + al| * |hi_in + lo_in - 1/(bh + bl)|
     < 2^-101.41 * |hi + lo| + |ah + al| * 2^-102.41 * (hi_in + lo_in)
     < 2^-101.41 * |hi + lo| + 2^-102.41 * (1 + 2^-101.41) * |hi + lo|
     < 2^-100.82 * |hi + lo| */
}

#if 0
// Algorithm DWDivDW1 from [3], with relative error bounded by 15u^2+56u^3
// (for rounding to nearest)
static inline void
d_div_acc1 (double *zh, double *zl, double xh, double xl, double yh, double yl)
{
  double th = xh / yh;
  double rh, rl;
  s_mul_acc1 (&rh, &rl, th, yh, yl);
  double pih, pil;
  two_sum (&pih, &pil, xh, -rh);
  double delta_h = pil - rl;
  double delta_l = delta_h + xl;
  double delta = pih + delta_l;
  double tl = delta / yh;
  fast_two_sum (zh, zl, th, tl);
}

// Algorithm DWDivDW2 from [3], with relative error bounded by 15u^2+56u^3
// (for rounding to nearest)
static inline void
d_div_acc2 (double *zh, double *zl, double xh, double xl, double yh, double yl)
{
  double th = xh / yh;
  double rh, rl;
  s_mul_acc1 (&rh, &rl, th, yh, yl);
  double pih = xh - rh;
  double delta_l = xl - rl;
  double delta = pih + delta_l;
  double tl = delta / yh;
  fast_two_sum (zh, zl, th, tl);
}
#endif

// Algorithm DWDivDW3 from [3], with relative error bounded by 9.8u^2
// (for rounding to nearest)
static inline void
d_div_acc3 (double *zh, double *zl, double xh, double xl, double yh, double yl)
{
  double th = 1.0 / yh;
  double rh = __builtin_fma (-yh, th, 1.0);
  double rl = -yl * th;
  double eh, el;
  fast_two_sum (&eh, &el, rh, rl);
  double delta_h, delta_l;
  s_mul_acc3 (&delta_h, &delta_l, th, eh, el);
  double mh, ml;
  // DWPlusFP (&mh, &ml, delta_h, delta_l, th);
  fast_sum_acc (&mh, &ml, th, delta_h, delta_l);
  d_mul_acc2 (zh, zl, xh, xl, mh, ml);
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

/* the following is a degree-7 odd polynomial approximating tanh(x)
   for |x| < 0.00543 generated by Sollya (see Ptanh.sollya), with
   maximal relative error 2^-71.98 and maximal absolute error 2^-80.528 */
static const double TT[] = {
  0x1p0,                 /* degree 1 */
  -0x1.5555555555553p-2, /* degree 3 */
  0x1.111111103f43cp-3,  /* degree 5 */
  -0x1.ba18e77264096p-5, /* degree 7 */
};

/* the following is a degree-13 odd polynomial approximating tanh(x)
   for |x| < 0.00543 generated by Sollya (see Ptanh2.sollya), with
   double-double coefficients and maximal relative error 2^-114.05,
   maximal absolute error 2^-123.381 */
static const double T2[][2] = {
  {0x1p0, 0},                                      /* degree 1 */
  {-0x1.5555555555555p-2, -0x1.5555555559905p-56}, /* degree 3 */
  {0x1.1111111111111p-3, 0x1.1111cecfa44acp-59},   /* degree 5 */
  {-0x1.ba1ba1ba1ba1cp-5, 0},                      /* degree 7 */
  {0x1.664f4882cbcc1p-6, 0},                       /* degree 9 */
  {-0x1.226e3e86a65ffp-7, 0},                      /* degree 11 */
  {0x1.dcc10cbf3d60fp-9, 0},                       /* degree 13 */
};

/* put in h+l a double-double approximation of sinh(w), for |w| < 0.00543,
   with maximal relative error 2^-67.58 (see analyze_eval_S_all(rel=true)
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
   with maximal absolute error 2^-67.02 (see analyze_eval_C() from
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

/* put in h+l a double-double approximation of tanh(w), for |w| < 0.00543,
   with maximal relative error 2^-66.53 (see analyze_eval_T_all(rel=true)
   from accompanying file tanh.sage) */
static void
eval_T (double *h, double *l, double w)
{
  double z = w * w;
  *h = __builtin_fma (TT[3], z, TT[2]);
  *h = __builtin_fma (*h, z, TT[1]);
  *h = *h * z; /* h approximates w^2*(TT[1]+w^2*TT[2]+w^4*TT[2]) */
  /* we use the fact that TT[0]=1 here, thus we add w + w*h */
  fast_two_sum (h, l, w, *h * w);
}

/* put in h+l a double-double approximation of tanh(w), for |w| < 0.00543 */
static void
eval_T2 (double *h, double *l, double w)
{
  double zh, zl;
  a_mul (&zh, &zl, w, w); /* zh+zl = w^2 */
  *h = __builtin_fma (T2[6][0], zh, T2[5][0]);
  *h *= zh;                                     /* multiply by w^2 */
  fast_two_sum (h, l, T2[4][0], *h);            /* add T2[4] */
  *h *= zh;                                     /* multiply by w^2 */
  fast_two_sum (h, l, T2[3][0], *h);            /* add T2[3] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  fast_sum2 (h, l, T2[2][0], T2[2][1], *h, *l); /* add T2[2] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  fast_sum2 (h, l, T2[1][0], T2[1][1], *h, *l); /* add T2[1] */
  d_mul (h, l, *h, *l, zh, zl);                 /* multiply by w^2 */
  s_mul (h, l, w, *h, *l);                      /* multiply by w */
  fast_sum (h, l, w, *h, *l);                   /* add w */
}

/* Put in h+l a double-double approximation of tanh(x),
   for 0 <= x <= 0x1.633ce8fb9f87dp+9.
   Return the absolute error bound:
   |h + l - sin(x)| < err. */
static double
cr_tanh_fast (double *h, double *l, double x)
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
     = T[i][1]*cosh(v) + T[i][2]*sinh(v).
     |v| = |x - T[i][0]| <= |x - i*2^8/magic| + |T[i][0] - i*2^8/magic|
                         <= |x - i*2^8/magic| + 8e-14
                         <= |x - k/magic| + |j/magic| + 8e-14
                         <= 0.00542055 + 255/magic + 8e-14
                         < 2.77 */
  double w = v - U[j][0];
  /* since v = U[j][0] + w, we approximate sinh(v) as
     sinh(U[j][0])*cosh(w) + cosh(U[j][0])*sinh(w)
     = U[j][1]*cosh(w) + U[j][2]*sinh(w), and cosh(v) as
     sinh(U[j][0])*sinh(w) + cosh(U[j][0])*cosh(w)
     = U[j][1]*sinh(w) + U[j][2]*cosh(w) */

  /* since |T[i][0] - i*2^8/magic| < 8e-14 and
           |U[j][0] - j/magic| < 9.14e-07, we have:
     |x - T[i][0] - U[j][0]| < 0.00542055 + 8e-14 + 9.14e-07 < 0.00543 */

  /* we have |w| < 0.00543 */

  if (k == 0)
  {
    eval_T (h, l, w); /* |h + l - tanh(w)| < 2^-66.53 * h */
    double err = 0x1.63p-67 * *h; /* 2^-66.53 < 0x1.63p-67 * *h */
    /* if err is zero due to rounding, we replace the error bound by 2^-1074 */
    return (err == 0) ? 0x1p-1074 : err;
  }
  
  double swh, swl, cwh, cwl;
  eval_S (&swh, &swl, w); /* |swh + swl - sinh(w)| < 2^-67.58*|swh+swl| */
  eval_C (&cwh, &cwl, w); /* |cwh + cwl - cosh(w)| < 2^-67.02*|cwh+cwl| */

  double svh, svl, cvh, cvl, h1, l1, h2, l2;
  s_mul (&h1, &l1, U[j][1], cwh, cwl); /* U[j][1]*cosh(w) */
  /* |U[j][1] - sinh(U[j][0])| < 2^-16 ulp(U[j][1]) <= 2^-68 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-67.02*|cwh+cwl| thus
     |h1+l1-sinh(U[j][0])*cosh(w)| < 2^-66.42*|h1+l1| */
  s_mul (&h2, &l2, U[j][2], swh, swl); /* U[j][2]*sinh(w) */
  /* |U[j][2] - cosh(U[j][0])| < 2^-16 ulp(U[j][2]) <= 2^-68 |U[j][2]|
     and |swh + swl - sinh(w)| < 2^-67.58*|swh+swl| thus
     |h2+l2-cosh(U[j][0])*sinh(w)| < 2^-66.77*|h2+l2| */

  fast_sum2 (&svh, &svl, h1, l1, h2, l2); /* svh+svl approximates sinh(v) */
  /* since h1+l1 and h2+l2 have a relative error bound < 2^-66.42, that bound
     holds for the sum of their absolute values, but we might have
     cancellation, the worst case being for j=1 and w=-0.00543,
     where h1+l1 >= 0.0108414. and h2+l2 >= -0.0054304,
     thus (|h1+l1| + |h2+l2|)/((|h1+l1| - |h2+l2|) < 3.008,
     thus |svh + svl - sinh(v)| < 3.008*2^-66.42 < 2^-64.83.
     Note: the rounding error in fast_sum2() is absorbed in the above
     error bound (which is over-estimated). */

  s_mul (&h1, &l1, U[j][1], swh, swl); /* U[j][1]*sinh(w) */
  /* |U[j][1] - sinh(U[j][0])| < 2^-16 ulp(U[j][1]) <= 2^-68 |U[j][1]|
     and |swh + swl - sinh(w)| < 2^-67.58*|swh+swl| thus
     |h1+l1-sinh(U[j][0])*sinh(w)| < 2^-66.77*|h1+l1| */
  s_mul (&h2, &l2, U[j][2], cwh, cwl); /* U[j][2]*cosh(w) */
  /* |U[j][2] - cosh(U[j][0])| < 2^-16 ulp(U[j][1]) <= 2^-68 |U[j][1]|
     and |cwh + cwl - cosh(w)| < 2^-67.02*|cwh+cwl| thus
     |h2+l2-cosh(U[j][0])*cosh(w)| < 2^-66.42*|h2+l2| */
  fast_sum2 (&cvh, &cvl, h2, l2, h1, l1); /* cvh+cvl approximates cosh(v) */
  /* since h1+l1 and h2+l2 have a relative error bound < 2^-66.42, that bound
     holds for the sum of their absolute values, but we might have
     cancellation, the worst case being for j=1 and w=-0.00543,
     where h2+l2 >= 1.0000735. and h1+l1 >= -0.0000589
     thus (|h1+l1| + |h2+l2|)/((|h1+l1| - |h2+l2|) < 1.000118,
     and |cvh + cvl - cosh(v)| < 1.000118*2^-66.42 < 2^-66.41.
     Note: he rounding errors in fast_sum2 are absorbed in the above error
     bound (which is over-estimated) */

  if (i == 0)
  {
    d_div (h, l, svh, svl, cvh, cvl);
    /* |h + l - (svh + svl) / (cvh + cvl)| < 2^-100.82 * |h + l| */
    /* since the relative error on svh + svl is bounded by e1=2^-64.83,
       that on cvh + cvl by e2=2^-66.41, and that on the division by
       e3=2^-100.82, the relative error on h+l is bounded by:
       (1+e1)*(1+e2)*(1+e3)-1 < 2^-64.41. */
    return 0x1.82p-65; /* 2^-64.41 < 0x1.82p-65 */
  }

  /* At this point svh+svl approximates sinh(v) with relative error bounded by
     2^-64.83, cvh+cvl approximates cosh(v) with relative error bounded
     by 2^-66.41, T[i][1]+T[i][2] approximates sinh(T[i][0]) with relative
     error bounded by 2^-107, T[i][3]+T[i][4] approximates cosh(T[i][0]) with
     relative error bounded by 2^-107, and we have to compute:
     (T[i][1]+T[i][2])*(cvh+cvl) + (T[i][3]+T[i][4])*(svh+svl) */

  double sh, sl, ch, cl;
  d_mul (&h1, &l1, T[i][1], T[i][2], cvh, cvl);
  /* |T[i][1] + T[i][2] - sinh(T[i][0])| < 2^-107 |T[i][1]|
     and |cvh + cvl - cosh(v)| < 2^-66.41*|cvh + cvl| thus
     |h1+l1-sinh(T[i][0])*cosh(v)| < 2^-66.40*|h1+l1| */
  d_mul (&h2, &l2, T[i][3], T[i][4], svh, svl);
  /* |T[i][3] + T[i][4] - cosh(T[i][0])| < 2^-107 |T[i][3]|
     and |svh + svl - sinh(v)| < 2^-64.83*|svh + svl| thus
     |h2+l2-cosh(T[i][0])*sinh(v)| < 2^-64.82*|h2+l2| */

  /* 2^-66.40 < 0x1.85p-67 and 2^-64.82 < 0x1.23p-65 */
  /* Warning: h2 might be negative if j=0 and w<0 (thus v=w) */
  /* Since |T[i][0]| > 2.77 and |v| < 2.77, we have
     |(h1+l1)/(h2+l2)| = |sinh(T[i][0])*cosh(v)/(cosh(T[i][0])*sinh(v))|
                       = |tanh(T[i][0])/tanh(v)| > 1,
     thus the relative error on h2+l2 also holds relatively to h1+l1:
     2^-66.40*|h1+l1| + 2^-64.82*|h2+l2| <= (2^-66.40+2^-64.82)*|h1+l1|
                                         <= 2^-64.40*|h1+l1| */
  fast_sum2 (&sh, &sl, h1, l1, h2, l2);
  /* the error in fast_sum2() is absorbed by the above errors, which are
     overestimated */

  /* |sh + sl - sinh(x)| < errs * |sh + sl| */

  d_mul (&h1, &l1, T[i][3], T[i][4], cvh, cvl); /* (T[i][3]+T[i][4])*(cvh+cvl) */
  /* |T[i][3] + T[i][4] - cosh(T[i][0])| < 2^-107 |T[i][3]|
     and |cvh + cvl - cosh(v)| < 2^-66.41*|cvh + cvl| thus
     |h1+l1-cosh(T[i][0])*cosh(v)| < 2^-66.40*|h1+l1| */
  d_mul (&h2, &l2, T[i][1], T[i][2], svh, svl); /* T[i][1]*(cvh+cvl+svh+svl) */
  /* |T[i][1] + T[i][2] - sinh(T[i][0])| < 2^-107 |T[i][1]|
     and |svh + svl - sinh(v)| < 2^-64.83*|svh + svl| thus
     |h2+l2-sinh(T[i][0])*sinh(v)| < 2^-64.82*|h2+l2| */
  /* Since |T[i][0]| > 2.77 and |v| < 2.77, we have
     |(h1+l1)/(h2+l2)| = |cosh(T[i][0])*cosh(v)/(sinh(T[i][0])*sinh(v))|
                         > 1,
     thus the relative error on h2+l2 also holds relatively to h1+l1:
     2^-66.40*|h1+l1| + 2^-64.82*|h2+l2| <= (2^-66.40+2^-64.82)*|h1+l1|
                                         <= 2^-64.40*|h1+l1| */

  fast_sum2 (&ch, &cl, h1, l1, h2, l2);
  /* the error in fast_sum2() is absorbed by the above errors, which are
     overestimated */

  d_div (h, l, sh, sl, ch, cl);
  /* the relative error on sl+sl is bounded by e1=2^-64.40, that on ch+cl is
     also bounded by e1, and the relative error in d_div() is bounded
     by e2=2^-100.82, thus the relative error on h+l is bounded by:
     (1+e1)^2*(1+e2)-1 < 2^-63.39. */
  return 0x1.87p-64 * *h; /* 2^-63.39 < 0x1.87p-64 */
}

/* return h + l which approximates tanh(s*x) where s in {-1,1} and x > 0 */
static void
cr_tanh_accurate (double *h, double *l, double x, double s)
{
  static const double magic = 0x1.70f77fc88ae3cp6;
  int k = __builtin_round (magic * x);
  int i = k >> 8, j = k & 0xff;
  double v = x - T[i][0];
  double w = v - U[j][0];
  double swh, swl, cwh, cwl;
  if (k == 0)
  {
    static double exceptions[][3] = {
      {0x1.e0000000000e1p-22, 0x1.dfffffffffeafp-22, -0x1.ffffffffffffep-76},
      {0x1.880072ccb1051p-11, 0x1.88006e032a37bp-11, -0x1.fffffffffffffp-65},
      {0x1.c3e1f9778db14p-14, 0x1.c3e1f95a3873dp-14, 0x1.fffffffffffffp-68},
    };
    for (int i = 0; i < 3; i++)
      if (x == exceptions[i][0])
      {
        *h = s * exceptions[i][1];
        *l = s * exceptions[i][2];
        return;
      }
    eval_T2 (h, l, s * w);
    return;
  }

  /* we approximate directly sinh on s*x and not x, since for rounding
     towards -Inf or +Inf, this requires fewer exceptional cases */
  eval_S2 (&swh, &swl, s * w);
  eval_C2 (&cwh, &cwl, w);
  double svh, svl, cvh, cvl;
  double h1, l1, h2, l2;
#define D_MUL d_mul_acc3
#define FAST_SUM2 fast_sum2_acc2
  D_MUL (&h1, &l1, s * U[j][1], s * Ul[j][0], cwh, cwl);
  D_MUL (&h2, &l2, U[j][2], Ul[j][1], swh, swl);
  FAST_SUM2 (&svh, &svl, h1, l1, h2, l2);
  D_MUL (&h1, &l1, U[j][2], Ul[j][1], cwh, cwl);
  D_MUL (&h2, &l2, s * U[j][1], s * Ul[j][0], swh, swl);
  FAST_SUM2 (&cvh, &cvl, h1, l1, h2, l2);
#undef D_MUL
#undef FAST_SUM2
  if (i == 0)
  {
    static double exceptions[][3] = {
      {0x1.1375c272d441cp+1, 0x1.f258bcd572d2ep-1, -0x1p-54},
      {0x1.048711422ed6ap-6, 0x1.048172571b3a6p-6, -0x1.605e314431363p-112},
      {0x1.07dd8d1b13572p-6, 0x1.07d7b62ca65b5p-6, -0x1.b220ffb9cf44ep-112},
      {0x1.1233157441ed4p-6, 0x1.122c87ec88f42p-6, 0x1.0243e241e6656p-110},
      {0x1.6cbd9c7c4f2bfp-8, 0x1.6cbca5af62bbp-8, 0x1.cdfec6e2ef7a8p-115},
      {0x1.b8aca02fd218bp-8, 0x1.b8aaecee6ef11p-8, -0x1.117c36e15a6fap-115},
      {0x1.100ddbaa215d8p+1, 0x1.f19c60b4238dcp-1, -0x1.bceaa50698d34p-106},
      {0x1.1b56ee2d1147ap+1, 0x1.f3e8da02e5a7ep-1, -0x1.605928d804bf5p-105},
      {0x1.22dd041738727p+1, 0x1.f53c882e6e949p-1, 0x1.024c4cc55688bp-104},
      {0x1.e7b151e50ace8p-8, 0x1.e7af03ed72a95p-8, -0x1.0782c932cf01ep-115},
      {0x1.7c8b110da7859p-8, 0x1.7c89f8c3da23fp-8, -0x1.3eadf765854c1p-113},
      {0x1.3fdc557c28f86p-6, 0x1.3fd1eeb4ae601p-6, -0x1.b88138013e24p-111},
      {0x1.6d2b05a7faa7fp-6, 0x1.6d1b8bb0726bbp-6, 0x1.1dd0272792e45p-110},
      {0x1.9e6b8533e155ep-6, 0x1.9e54e68508294p-6, 0x1.e3e90de57dbb7p-110},
      {0x1.16a973d843b7bp-3, 0x1.14f4723e2559p-3, 0x1.0c59ea55f11c4p-110},
      {0x1.2c12ce474d94ep-3, 0x1.29f1c403195d2p-3, 0x1.086a19b08f101p-107},
      {0x1.36f33d51c264dp-2, 0x1.2dbb7b1c91363p-2, -0x1.19ab4b485c24ap-106},
      {0x1.cd4104178518dp-4, 0x1.cb50660b00f34p-4, -0x1.eef7c2f919739p-110},
      {0x1.3a0dd2ce4b0ebp+0, 0x1.aeeb8d9e19d8bp-1, -0x1.18ba74f0b2fb7p-105},
      {0x1.9d343db61d049p+0, 0x1.d8f7213aabd6p-1, -0x1.cf2eafb8aa804p-106},
      {0x1.e611aa58ab608p-2, 0x1.c493dc899e4a5p-2, 0x1.6a0755d0e11bap-106},
      {0x1.b7ae824c71873p-1, 0x1.64279e7c7064bp-1, 0x1.aa5585668bfcap-106},
      {0x1.77194793ac4e4p+0, 0x1.cc1d9fbf0b50bp-1, -0x1.5c5e72ab75cf2p-105},
      {0x1.43eaea23649c3p-2, 0x1.39877ed02864p-2, 0x1.813581f4a00ap-106},
      {0x1.c350f47039524p+0, 0x1.e2bae0e9932a4p-1, -0x1.b28882ea8585dp-105},
      {0x1.fc8be93ceb1c7p-8, 0x1.fc894c4eeec42p-8, -0x1.0e77db9398ep-114},
      {0x1.0a90d36059a04p-6, 0x1.0a8ace14dbccep-6, -0x1p-59},
      {0x1.1005ec0bccabbp-1, 0x1.f20b1c8557decp-2, 0x1.ffffffffffffdp-56},
      {0x1.b09abe7d49101p+0, 0x1.de4639cd633c8p-1, 0x1.d81ac695bc61ep-105},
      {0x1.6961b00c6bb74p-5, 0x1.6925b8f7dd412p-5, 0x1.a28200237b0c2p-108},
      {0x1.34edbba55c659p-3, 0x1.329b4e3592dbap-3, -0x1.d1d7da2938b46p-107},
      {0x1.95e87aa828542p-6, 0x1.95d33972a6c13p-6, -0x1.336f680f198bap-109},
    };
    for (int i = 0; i < 32; i++)
      if (x == exceptions[i][0])
      {
        *h = s * exceptions[i][1];
        *l = s * exceptions[i][2];
        return;
      }
    d_div_acc3 (h, l, svh, svl, cvh, cvl);
    return;
  }

  static double exceptions[][3] = {
    {0x1.6ad3b64b0cde7p+2, 0x1.fffce0e97b779p-1, 0x1.7d7656fe4a9b4p-106},
    {0x1.731aceeabdd36p+1, 0x1.fce790a1abc82p-1, 0x1.0359bb2d02662p-105},
    {0x1.4ab67173af9f6p+3, 0x1.ffffffeddf35ep-1, -0x1.1c63f9095d02bp-105},
    {0x1.6396dc3669441p+3, 0x1.fffffffc2bb54p-1, -0x1.4cc9edaf2718fp-105},
    {0x1.a31a6a969867dp+1, 0x1.fe8912c0a8906p-1, 0x1.0a03176bcc80ep-106},
    {0x1.dd09782a03db2p+1, 0x1.ff683ad7e7bb3p-1, -0x1.b2dcdc25acd52p-107},
    {0x1.dd3c9b972acd7p+3, 0x1.ffffffffff82dp-1, -0x1.ffffffffffffdp-55},
    {0x1.786262131633dp+2, 0x1.fffdf4e17fe7cp-1, -0x1.106ebb7ee0843p-106},
    {0x1.eb528eaa42d77p+1, 0x1.ff8693f7fec43p-1, -0x1.47fef79e1c5afp-105},
    {0x1.50ff5f7e1bef5p+3, 0x1.fffffff3c2b92p-1, 0x1.703f53a5cd6cbp-105},
    {0x1.ef2723685a47dp+1, 0x1.ff8da0fae6272p-1, -0x1.6b5f031c1f1e8p-107},
    {0x1.03289aa1045dbp+3, 0x1.fffff9cd03f79p-1, -0x1.8e1efe42c8c3p-111},
    {0x1.4150814caced6p+3, 0x1.ffffffdf61ea3p-1, 0x1.873ae3f09c66p-105},
    {0x1.f879d6c96ec7cp+1, 0x1.ff9d20874d493p-1, -0x1.dd65bf7c5120ap-106},
    {0x1.647e53a8705a3p+2, 0x1.fffc3204f4897p-1, 0x1.0cfd57c4e2cc2p-104},
    {0x1.b5db40898efd3p+2, 0x1.ffffb36232778p-1, -0x1.6259be540ae3p-104},
  };
  for (int i = 0; i < 16; i++)
    if (x == exceptions[i][0])
    {
      *h = s * exceptions[i][1];
      *l = s * exceptions[i][2];
      return;
    }

  double sh, sl, ch, cl;
  /* svh+svl approximates sinh(s*x) */
#define D_MUL d_mul_acc3
#define FAST_SUM2 fast_sum2_acc1
  D_MUL (&h1, &l1, s * T[i][1], s * T[i][2], cvh, cvl);
  D_MUL (&h2, &l2, T[i][3], T[i][4], svh, svl);
  FAST_SUM2 (&sh, &sl, h1, l1, h2, l2);
  D_MUL (&h1, &l1, T[i][3], T[i][4], cvh, cvl);
  D_MUL (&h2, &l2, s * T[i][1], s * T[i][2], svh, svl);
  FAST_SUM2 (&ch, &cl, h1, l1, h2, l2);
#undef D_MUL
#undef FAST_SUM2
  d_div_acc3 (h, l, sh, sl, ch, cl);
}

#define MASK 0x7fffffffffffffff /* to mask the sign bit */

double
cr_tanh (double x)
{
  d64u64 v = {.f = x};
  int e = (v.u >> 52) - 0x3ff;
  int s = v.u >> 63; /* sign bit */
  v.u &= (uint64_t) MASK; /* get absolute value */

  if (e == 0x400 || e == 0xc00 || v.f >= 0x1.633ce8fb9f87ep+9)
    /* NaN or tanh(x) rounds to +/- 1 */
  {
    if (v.u > 0x7ff0000000000000) /* NaN */
      return x;
    else if (e == 0x400) /* +Inf */
      return 1.0;
    else if (e == 0xc00) /* -Inf */
      return -1.0;
    else if (s == 0)
      return 1.0 - 0x1p-54;
    else
      return -1.0 + 0x1p-54;
  }
  
  double h, l;
  double err = cr_tanh_fast (&h, &l, v.f);
  static double sign[] = { 1.0, -1.0 };
  h *= sign[s];
  l *= sign[s];
  
  double left  = h + (l - err);
  double right = h + (l + err);
  if (left == right)
    return left;

  /* Special case for small numbers, to avoid underflow issues in the accurate
     path: for |x| <= 0x1.d12ed0af1a27fp-27, |tanh(x) - x| < ulp(x)/2,
     and the Taylor expansion is tanh(x) = x - x^3/3 + O(x^5),
     thus tanh(x) rounds to the same value than x - 2^-54*x. */
  if (v.f <= 0x1.d12ed0af1a27fp-27)
    return __builtin_fma (x, -0x1p-54, x);

  /* also, for |x| >= 0x1.30fc1931f09cap+4, tanh(x) rounds to -1 or +1
     for rounding to nearest */

  cr_tanh_accurate (&h, &l, v.f, sign[s]);
  return h + l;
}
