/* Correctly-rounded error function for binary64 value.

Copyright (c) 2023 Paul Zimmermann

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

typedef union {double f; uint64_t u;} b64u64_u;

/* for 1 <= i < 95, p = C[i-1] is a degree-10 polynomial approximation of
   erf(i/16+1/32+x) for -1/32 <= x <= 1/32, where the coefficients of
   degree 0 and 1 are double-double: the coefficient of degree 0 is p[0]+p[1],
   that of degree 1 is p[2]+p[3], that of degree 2 is p[4], ..., that of
   degree 10 is p[12] */
static const double C[94][13] = {
   {0x1.207d480e90658p-4, 0x1.4c37de3fbd9b4p-58, 0x1.1fbd27cdc72d3p+0, -0x1.212bce2e695f8p-54, -0x1.1fbd27cdc72d3p-4, -0x1.7ca791fd8f7e3p-2, 0x1.1efd545de992p-5, 0x1.c532b7bbd364cp-4, -0x1.7da7bd04d6058p-7, -0x1.ac35a6ffdaa2bp-6, 0x1.7cab0e27fdefbp-9, 0x1.4a4282240d152p-8, -0x1.3d349039dbaaep-11}, /* j=1 68.936 */
   {0x1.1f5e1a35c3b89p-3, 0x1.d0b53b5ac6b55p-57, 0x1.1c62fa1e869b6p+0, 0x1.ce5fa4eebb5e8p-55, -0x1.1c62fa1e869b6p-3, -0x1.6f552dbcc3334p-2, 0x1.196c9cd8df1d6p-4, 0x1.aaba623e058c6p-4, -0x1.734ea67e0c512p-6, -0x1.89258c8fcff08p-6, 0x1.6f661b3744008p-8, 0x1.2773b9d46b4edp-8, -0x1.16a7e2a85be09p-10}, /* j=2 69.481 */
   {0x1.ac45e37fe2526p-3, 0x1.48d76c0005fdbp-57, 0x1.16e2d7093cd8cp+0, 0x1.978d07b54d8a2p-54, -0x1.a254428ddb453p-3, -0x1.59b3da8e1e173p-2, 0x1.988648fe89781p-4, 0x1.80342730feea1p-4, -0x1.09e7bcffeafbbp-5, -0x1.516b1e8701d59p-6, 0x1.038efd1a46cf9p-7, 0x1.e168905c4e108p-9, -0x1.aa1779bb095e3p-10}, /* j=3 69.188 */
   {0x1.1af54e232d609p-2, -0x1.bee802eb7031bp-56, 0x1.0f5d1602f7e41p+0, -0x1.3e4b7912a20c6p-55, -0x1.0f5d1602f7e41p-2, -0x1.3c974458cbdf5p-2, 0x1.040e8a6d83692p-3, 0x1.47e5cfee39e14p-4, -0x1.4c0b2553174bap-5, -0x1.08d945d3c34e6p-6, 0x1.3db856084d345p-7, 0x1.554342e09d916p-9, -0x1.f5f108c7333b6p-10}, /* j=4 69.974 */
   {0x1.5da9f415ff23fp-2, -0x1.a744e51e032a8p-59, 0x1.05fd3ecbec298p+0, -0x1.f178df9b57eep-54, -0x1.477c8e7ee733dp-2, -0x1.1917b60acab73p-2, 0x1.322a728d4d829p-3, 0x1.04c50a9cd3017p-4, -0x1.7ce764bad5915p-5, -0x1.68aac53a7eb66p-7, 0x1.62a71cfa217c8p-7, 0x1.6bf09fb398fcep-10, -0x1.e5d686641405fp-10}, /* j=5 68.952 */
   {0x1.9dd0d2b721f39p-2, -0x1.1673a032b6936p-56, 0x1.f5f0cdaf15313p-1, 0x1.e15128140592p-60, -0x1.78749a434fe4ep-2, -0x1.e106c51d1ef9fp-3, 0x1.5529abccff891p-3, 0x1.7488b8a7fa43ap-5, -0x1.9a7945aaefd5ep-5, -0x1.65c10d6ad286p-8, 0x1.70986f2e53cc3p-7, 0x1.04e9fd70c91b4p-13, -0x1.f326511f01137p-10}, /* j=6 69.783 */
   {0x1.db081ce6e2a48p-2, -0x1.7ff0d60f33e5dp-56, 0x1.dd167c4cf9d2ap-1, 0x1.44fb295629733p-55, -0x1.a173acc35a985p-2, -0x1.889a80f4ad958p-3, 0x1.6c2eea0d179bfp-3, 0x1.b064543900df7p-6, -0x1.a3fd9fc806c1dp-5, 0x1.060afffc88f26p-13, 0x1.678b12e27752dp-7, -0x1.1dc84e67cd219p-10, -0x1.e34391d1af066p-10}, /* j=7 73.148 */
   {0x1.0a7ef5c18edd2p-1, 0x1.5e7e49ed35672p-56, 0x1.c1efca49a5011p-1, 0x1.4c1996c347cap-55, -0x1.c1efca49a5011p-2, -0x1.2bf531866e01p-3, 0x1.76f27de8086cep-3, 0x1.dfeeb5a496f12p-8, -0x1.99f13afbc09c8p-5, 0x1.623c5b784ef76p-8, 0x1.493a75e3a10d2p-7, -0x1.1be811ffcb31dp-9, -0x1.796af6b8030bp-10}, /* j=8 69.833 */
   {0x1.25b8a88b6dd7fp-1, 0x1.9533a35946784p-55, 0x1.a5074e215762p-1, 0x1.fafdcd682cfbep-56, -0x1.d9a837e5824e4p-2, -0x1.9c41d1d5fae5fp-4, 0x1.75bebc1b17e4cp-3, -0x1.6410ad92c8a13p-7, -0x1.7df888e63cab2p-5, 0x1.4a547dcb549d8p-7, 0x1.18f102d648762p-7, -0x1.8d2d7a0516941p-9, -0x1.18fb96848e37bp-10}, /* j=9 70.192 */
   {0x1.3f196dcd0f135p-1, -0x1.f260226cfe6eep-56, 0x1.86e9694134b9ep-1, -0x1.3bc8c63379853p-55, -0x1.e8a3c39181e85p-2, -0x1.c8105021682f7p-5, 0x1.6963c8a39d07dp-3, -0x1.c1242dffc4ea6p-6, -0x1.52b2667f67d31p-5, 0x1.c7cd9838961b4p-7, 0x1.b62d4de056806p-8, -0x1.dbfdb57fa9141p-9, -0x1.6928387cac79ap-11}, /* j=10 71.578 */
   {0x1.569243d2b3a9bp-1, 0x1.8eec6a5ba0b5ep-56, 0x1.681ff24b4ab04p-1, -0x1.db063d51d4bf8p-58, -0x1.ef2bed2786b25p-2, -0x1.a4254557d7287p-7, 0x1.532415c266308p-3, -0x1.558b4c558b57ap-5, -0x1.1b7ad57fe6d5p-5, 0x1.1201d1bfd61acp-6, 0x1.298e96a40b074p-8, -0x1.0275ad9bbed08p-8, -0x1.980dbf2fd3a47p-14}, /* j=11 69.822 */
   {0x1.6c1c9759d0e5fp-1, 0x1.b1437873b033bp-55, 0x1.492e42d78d2c5p-1, -0x1.8bc60c4b6e49cp-55, -0x1.edc5644353c27p-2, 0x1.b6e8591f66e39p-6, 0x1.349b5eaa149dap-3, -0x1.b42a1890b5e54p-5, -0x1.b847797be92a5p-6, 0x1.2e0af94c2eb6ap-6, 0x1.2db5fe8653516p-9, -0x1.03f44f015b8d7p-8, 0x1.42777b1a2d7bcp-13}, /* j=12 72.292 */
   {0x1.7fb9bfaed8078p-1, 0x1.66d21928af29ap-56, 0x1.2a8dcede3673bp-1, -0x1.7367ce246f87ep-56, -0x1.e5267029187cp-2, 0x1.fe0796bb9d04dp-5, 0x1.0fa23021ae704p-3, -0x1.fa21ebca6342ap-5, -0x1.31546dcb782e3p-6, 0x1.37e545333bc77p-6, 0x1.0a634c702bf4dp-13, -0x1.e7fc385c708c9p-9, 0x1.b2dc1d23ddc02p-12}, /* j=13 69.99 */
   {0x1.91724951b8fc6p-1, -0x1.2792e4f680df4p-55, 0x1.0cab61f084b93p-1, 0x1.09a062926ac56p-56, -0x1.d62beb64e8441p-2, 0x1.7c9d756a115b5p-4, 0x1.cc60567d7596ep-4, -0x1.1350f4b21a8e7p-4, -0x1.53bb495ebefe7p-7, 0x1.30ac207288b4dp-6, -0x1.e415f13cb09a7p-10, -0x1.aabea5c7828b5p-9, 0x1.1d2ac801f4996p-10}, /* j=14 69.87 */
   {0x1.a1551a16aaeafp-1, 0x1.a554984153968p-57, 0x1.dfca26f5bbf88p-2, -0x1.ddb365f66c39cp-57, -0x1.c1cd84866038fp-2, 0x1.e4c9975da0983p-4, 0x1.747e31bf45d07p-4, -0x1.1d1f00109975fp-4, -0x1.47653f2307809p-9, 0x1.1a817bab992edp-6, -0x1.cb948fe45c5cep-9, -0x1.57a1921f1253cp-9, 0x1.418b627f83584p-10}, /* j=15 70.687 */
   {0x1.af767a741088bp-1, -0x1.c97d9c45b7eb2p-56, 0x1.a911f096fbc26p-2, -0x1.086b9634ce92bp-56, -0x1.a911f096fbc26p-2, 0x1.1b614b0f52819p-3, 0x1.1b614b0f54397p-4, -0x1.1b614b0f5182ep-4, 0x1.2e45a45416256p-8, 0x1.f096fd0fc74f2p-7, -0x1.390cc6d6dded5p-8, -0x1.ee22d5cc817fcp-10, 0x1.264f2dd5038f6p-10}, /* j=16 70.86 */
   {0x1.bbef0fbde6221p-1, -0x1.322bbc4f2422p-55, 0x1.75a91a7f4d2edp-2, 0x1.6eb6d97d284a8p-58, -0x1.8d03ac274201cp-2, 0x1.3954778d6a0dfp-3, 0x1.88e0f7b185374p-5, -0x1.0f7c15f75fa67p-4, 0x1.5e22cf70ab43bp-7, 0x1.9ad28c99210c1p-7, -0x1.704d28940a68ep-8, -0x1.233f7ed7764p-10, 0x1.3a0071f9d96aap-10}, /* j=17 72.385 */
   {0x1.c6dad2829ec62p-1, -0x1.ab77707fde6e4p-57, 0x1.45e99bcbb7915p-2, 0x1.7bc98a65f62cep-56, -0x1.6ea6cf452e838p-2, 0x1.4cb3cf0aa0b9cp-3, 0x1.ca5083167904ep-6, -0x1.f65d15f1d53f2p-5, 0x1.fd1c6c2837e1fp-7, 0x1.3acc7927ab9a5p-7, -0x1.8b447edbdc042p-8, -0x1.7a2a754f0aabbp-12, 0x1.32a4fa16169c4p-10}, /* j=18 73.545 */
   {0x1.d0580b2cfd249p-1, 0x1.4fc9c51b942b8p-55, 0x1.1a0dc51a9934dp-2, -0x1.ca98db3cd091p-57, -0x1.4ef05a0f95eebp-2, 0x1.5648b5dc47418p-3, 0x1.40fbaba43be34p-7, -0x1.c0db89d0ad671p-5, 0x1.388c3eedba984p-6, 0x1.aecb773596c4cp-8, -0x1.8bcd4e2ac3ebp-8, 0x1.4a21c8eafb7abp-12, 0x1.17b08478792d8p-10}, /* j=19 71.555 */
   {0x1.d8865d98abe01p-1, -0x1.fced1feaf1054p-55, 0x1.e4652fadcb6b2p-3, -0x1.eaa19413902cp-61, -0x1.2ebf3dcc9f22fp-2, 0x1.571d01c5c56cp-3, -0x1.93a9a7bb9777bp-8, -0x1.8281ce0b40a7ep-5, 0x1.5d0003eca5316p-6, 0x1.db43d0686d4dcp-9, -0x1.756b323d07612p-8, 0x1.cbafd6d2a81f6p-11, 0x1.d05a5db8099e1p-11}, /* j=20 71.149 */
   {0x1.df85ea8db188ep-1, -0x1.f71ed6189d2f4p-55, 0x1.9cb5bd549b111p-3, -0x1.97e38ef637808p-59, -0x1.0ed7443f85c33p-2, 0x1.5066cda84bbacp-3, -0x1.419fa10b7143ap-6, -0x1.3f41761d69b8ap-5, 0x1.6d1d7263bf2e3p-6, 0x1.e378168650dc3p-11, -0x1.4ccaabb453e45p-8, 0x1.54fee1c5b662cp-10, 0x1.392a8eede7093p-11}, /* j=21 72.507 */
   {0x1.e5768c3b4a3fcp-1, 0x1.8b631ce44fc8cp-57, 0x1.5ce595c455b0ap-3, 0x1.c2a536e9ef344p-59, -0x1.dfbbadedf5d2ep-3, 0x1.4374d82e04c69p-3, -0x1.f3b8d52d3416dp-6, -0x1.f572c4c8e1cefp-6, 0x1.6b16f51f51e91p-6, -0x1.73ff9aa31f7bep-10, -0x1.173f6c5d28505p-8, 0x1.9d670dc627a21p-10, 0x1.3eca19b3a6dap-12}, /* j=22 73.402 */
   {0x1.ea7730ed0bbb9p-1, 0x1.2c5bc77196fccp-55, 0x1.24a7b84d38971p-3, 0x1.a9f6a25c529aep-57, -0x1.a4b118ef01593p-3, 0x1.319c7a75f9189p-3, -0x1.3db5bed47fec1p-5, -0x1.7019bda6d9706p-6, 0x1.59d3aa44e6dbcp-6, -0x1.b324e4180e7eap-9, -0x1.b477ed1835d0ap-9, 0x1.bfdd351d007d7p-10, 0x1.c453a432f4a38p-14}, /* j=23 74.839 */
   {0x1.eea5557137aep-1, -0x1.385da06788ad3p-55, 0x1.e723726b824a9p-4, -0x1.22814164756fep-59, -0x1.6d5a95d0a1b7fp-3, 0x1.1c2a02beb6abap-3, -0x1.6d5a95d09f598p-5, -0x1.e723726baf415p-7, 0x1.3ca3d6fd2d50dp-6, -0x1.36d7375438959p-8, -0x1.35a81e797731bp-9, 0x1.c008ee392a3e8p-10, -0x1.576cd5827631bp-13}, /* j=24 71.591 */
   {0x1.f21c9f12f0677p-1, -0x1.7f01e661a7ffp-58, 0x1.92470a61b6965p-4, 0x1.c68a5ee00ee26p-58, -0x1.3a47801c56a57p-3, 0x1.0453f90d3bd36p-3, -0x1.8a7c6a2393c6ep-5, -0x1.075c0880503a8p-7, 0x1.16f9c9e2ad7f2p-6, -0x1.74c2fa71feb9p-8, -0x1.7615241d69fbcp-10, 0x1.a3ab1b99c4c6ap-10, -0x1.9b1c5d6416685p-13}, /* j=25 72.095 */
   {0x1.f4f693b67bd77p-1, -0x1.3a1f528bb29fp-56, 0x1.499d478bca735p-4, 0x1.3179dfdc2f59cp-60, -0x1.0bcfca21947dbp-3, 0x1.d6631e1a28eap-4, -0x1.974c036867fbbp-5, -0x1.17d430bda28a6p-9, 0x1.d857f2e562c11p-7, -0x1.954aaf96800d6p-8, -0x1.2e4c886d9eba3p-11, 0x1.71b71f2b61366p-10, -0x1.5b4028e33911bp-12}, /* j=26 73.143 */
   {0x1.f74a6d9a38383p-1, 0x1.c33a010bd3ef5p-55, 0x1.0bf97e95f2a64p-4, -0x1.44667883b7324p-58, -0x1.c435059d09788p-4, 0x1.a3687c1eaf1aep-4, -0x1.9647a30b1739ep-5, 0x1.6981061dda567p-9, 0x1.7e8755f6d887ap-7, -0x1.9be7315952735p-8, 0x1.3a7793562a09fp-13, 0x1.3194b00daf15ap-10, -0x1.a1c5e7ae94298p-12}, /* j=27 73.344 */
   {0x1.f92d077f8d56dp-1, 0x1.8b565a3623918p-56, 0x1.b055303221015p-5, 0x1.cbb9d60148e5ep-59, -0x1.7a4a8a2bdce13p-4, 0x1.7148c3d57c312p-4, -0x1.8a0da543054fap-5, 0x1.b22257dd049bdp-8, 0x1.25b378aa7c763p-7, -0x1.8d10fb7595e38p-8, 0x1.7ec9dea566a42p-11, 0x1.d4ce701d9bc6cp-11, -0x1.e8585a3339a1cp-12}, /* j=28 73.239 */
   {0x1.fab0dd89d1309p-1, -0x1.ae62205ede897p-55, 0x1.5a08e85af27ep-5, 0x1.e51eef216f609p-59, -0x1.399812926bc23p-4, 0x1.4140efb719cafp-4, -0x1.7535a61a43015p-5, 0x1.374c88c7f040fp-7, 0x1.a4070a51b5ef7p-8, -0x1.6dc0791cd1264p-8, 0x1.2edf58ea12d9cp-10, 0x1.4501acec55759p-11, -0x1.95ef9debfdc3bp-12}, /* j=29 72.358 */
   {0x1.fbe61eef4cf6ap-1, 0x1.15deab46393cp-55, 0x1.12ceb37ff9bc3p-5, 0x1.a3cddc85018c9p-59, -0x1.01a1c847fa207p-4, 0x1.143d1c6f4f092p-4, -0x1.5a316520b8c33p-5, 0x1.779b1e5710395p-7, 0x1.0d099caf484cdp-8, -0x1.42fcbb807cb2bp-8, 0x1.76fbea2261eadp-10, 0x1.7c0ffc00cf5e8p-12, -0x1.846bee5b5239bp-12}, /* j=30 73.487 */
   {0x1.fcdacca0bfb73p-1, -0x1.2c33df20498aap-55, 0x1.b1160991ff737p-6, -0x1.d8935d72775eep-61, -0x1.a38d59456f77dp-5, 0x1.d5bd91b6b012p-5, -0x1.3b35dcbc802cdp-5, 0x1.9d76b0a063994p-7, 0x1.14c887b77663dp-9, -0x1.117f43ea41e13p-8, 0x1.9b46fb3498777p-10, 0x1.1dacb8f3a280ep-13, -0x1.58a6cf31dea7ep-12}, /* j=31 76.169 */
   {0x1.fd9ae142795e3p-1, 0x1.9727d0864a694p-56, 0x1.529b9e8cf9a1ep-6, 0x1.b52f269adb3d8p-61, -0x1.529b9e8cf9a1ep-5, 0x1.8b0ae3a47892p-5, -0x1.1a2c59757b0c5p-5, 0x1.ace7404c3c8ebp-7, 0x1.e1936069aca65p-12, -0x1.bae0ae0900699p-9, 0x1.a112599fc37d4p-10, -0x1.a20d1adcec486p-15, -0x1.08682454c00b2p-12}, /* j=32 74.368 */
   {0x1.fe307f2b503dp-1, -0x1.8a55c80920bp-57, 0x1.06ae13b0d3255p-6, -0x1.88536e932e76ap-60, -0x1.0ee3844e59be7p-5, 0x1.48b127f8ed8a2p-5, -0x1.f155b4e7d9a1bp-6, 0x1.aa2c0753e6ae4p-7, -0x1.bbf7e1655a93dp-11, -0x1.5478d9eb75845p-9, 0x1.8eabc4a0dd4dfp-10, -0x1.91fda35d0719p-13, -0x1.75b4d79b340cap-13}, /* j=33 74.082 */
   {0x1.fea4218d6594ap-1, -0x1.e332874bae348p-58, 0x1.94624e78e0fafp-7, -0x1.40e7bc548360cp-61, -0x1.ada873606f0aap-6, 0x1.0ea475da3be7cp-5, -0x1.afe553fa444b5p-6, 0x1.9973b48a09e01p-7, -0x1.dd78ee7e27813p-10, -0x1.ea03bf89cab6ap-10, 0x1.6aa47a44f2b07p-10, -0x1.2f95fd4a016a5p-12, -0x1.16460fd0e4f1bp-13}, /* j=34 74.478 */
   {0x1.fefcce6813974p-1, -0x1.b27c99b610c8p-58, 0x1.34d7dbc76d7e5p-7, 0x1.37e4dfb157ea7p-61, -0x1.51cc18621fc23p-6, 0x1.b925a99886bb3p-6, -0x1.71e7d408c8726p-6, 0x1.7ea58080b455cp-7, -0x1.46eb9d3a915cep-9, -0x1.403336fa63325p-10, 0x1.3b38b3b2ab455p-10, -0x1.6ba72730d8eb8p-12, -0x1.18d663efafa5cp-14}, /* j=35 75.472 */
   {0x1.ff404760319b4p-1, 0x1.f142283eb5aa4p-56, 0x1.d4143a9dfe965p-8, -0x1.6d75137201f1ep-63, -0x1.074b60f8df349p-6, 0x1.63ef61e824421p-6, -0x1.38a983278893ap-6, 0x1.5d3b17bbeb369p-7, -0x1.7cae0d6be244p-9, -0x1.5f83199535244p-11, 0x1.06042501e336ap-10, -0x1.8323d9cdfa7e1p-12, -0x1.39fc3b29f67a4p-16}, /* j=36 74.941 */
   {0x1.ff733814af88cp-1, 0x1.0a873a5aba582p-56, 0x1.5ff2750fe782p-8, -0x1.5e7b2e4b7af26p-62, -0x1.96f0575a63ae5p-7, 0x1.1c5a643f0436p-6, -0x1.04f5caaf21428p-6, 0x1.382a146b03854p-7, -0x1.95cab954ecafbp-9, -0x1.d2fda2b719ddap-13, 0x1.9f52b6ad888eep-11, -0x1.7d5082c3acc74p-12, 0x1.781f01f90ed59p-16}, /* j=37 75.479 */
   {0x1.ff9960f3eb327p-1, -0x1.08b1cf628d264p-56, 0x1.06918b6355624p-8, 0x1.221c18b6f9fd3p-62, -0x1.37ccd585f564bp-7, 0x1.c1ec102e364edp-7, -0x1.ae59615f8ed62p-7, 0x1.11dae47356141p-7, -0x1.982b2745de004p-9, 0x1.0283d50ac535p-13, 0x1.377be760426c8p-11, -0x1.61d5696e00404p-12, 0x1.d36d81df7d985p-15}, /* j=38 77.625 */
   {0x1.ffb5bdf67fe6fp-1, 0x1.4e812085f55p-62, 0x1.84ba3004a50dp-9, -0x1.900c6c588d32p-64, -0x1.d9c2ea85a927dp-8, 0x1.60898536e1048p-7, -0x1.5eb1c899f0ee8p-7, 0x1.d854f73e79f52p-8, -0x1.897719a139c0cp-9, 0x1.88cdc5b467ceep-12, 0x1.b32481cff494fp-12, -0x1.38091246aa3ecp-12, 0x1.31bf1fdd937c1p-14}, /* j=39 77.078 */
   {0x1.ffcaa8f4c9beap-1, 0x1.b0cedc5c46b72p-55, 0x1.1d83170fbf6fbp-9, 0x1.ea7a26415dc6ep-63, -0x1.64e3dcd3af4bap-8, 0x1.119da0c46ccbp-7, -0x1.1a89b97ceb10ap-7, 0x1.90e81283fcfbbp-8, -0x1.6ecdbf5c4c494p-9, 0x1.1c610b39a5391p-11, 0x1.11539cd9dabeap-12, -0x1.066d673f433b4p-12, 0x1.5376fa55a886ep-14}, /* j=40 76.671 */
   {0x1.ffd9f78c7524ap-1, 0x1.04ed6bfff564fp-55, 0x1.a024365f771bdp-10, 0x1.3c6d88f601343p-64, -0x1.0a9732d5284ddp-8, 0x1.a4bf47a43042bp-8, -0x1.c23802d8a630ep-8, 0x1.4f400706670d3p-8, -0x1.4c9a2c94b0edbp-9, 0x1.4f7a510d34576p-11, 0x1.18adebe879de2p-13, -0x1.a4c92895c44acp-13, 0x1.5295ca3aa8179p-14}, /* j=41 77.008 */
   {0x1.ffe514bbdc197p-1, -0x1.cd9622e5f0e8p-58, 0x1.2ce898809244ep-10, 0x1.08b54892f8bbdp-64, -0x1.8af14828bffa7p-9, 0x1.407fbd18f1203p-8, -0x1.62d4c6d49c3e3p-8, 0x1.146c4b3e2585cp-8, -0x1.267f3bc8a7c1bp-9, 0x1.64f8929e6ba18p-11, 0x1.22a3f12abab0cp-15, -0x1.401599dae3c58p-13, 0x1.3262c520dba3ap-14}, /* j=42 77.958 */
   {0x1.ffed167b12ac2p-1, -0x1.ddc0cfaec02d1p-55, 0x1.afc85e0f82e12p-11, 0x1.414793dd2d2c2p-66, -0x1.221a9f326bef4p-9, 0x1.e3c9aab90bcfap-9, -0x1.14b1b981421c9p-8, 0x1.c1c19b9e5b3bep-9, -0x1.feac3db81b855p-10, 0x1.63e882b8946dap-11, -0x1.44453ddf7eb35p-15, -0x1.c8f75b09cfbecp-14, 0x1.14eecc82399b9p-14}, /* j=43 78.451 */
   {0x1.fff2cfb0453d9p-1, 0x1.9a913753bb66dp-55, 0x1.3360ccd23db3ap-11, 0x1.21b433a1a627p-69, -0x1.a6a519a114d7p-10, 0x1.69cf466ccdf69p-9, -0x1.ab0c273ac2344p-9, 0x1.693596065aa85p-9, -0x1.b2755c0127fe6p-10, 0x1.52b627b041bd1p-11, -0x1.75567893cc61p-14, -0x1.2ae91d4820426p-14, 0x1.ca65df142df07p-15}, /* j=44 79.188 */
   {0x1.fff6dee89352ep-1, 0x1.b96c0a45ff8e1p-55, 0x1.b23a5a23e421p-12, 0x1.6cbd07beb869ap-67, -0x1.315107613c673p-10, 0x1.0c243329a9ca8p-9, -0x1.4630116262589p-9, 0x1.1e84d10224829p-9, -0x1.6b418720d8781p-10, 0x1.36eddfc737bb2p-11, -0x1.f788b2fd2a97ep-14, -0x1.5248bd0aa4efdp-15, 0x1.7610402e71234p-15}, /* j=45 78.516 */
   {0x1.fff9ba420e834p-1, 0x1.1379ed3aba56ep-56, 0x1.30538fbb77ecdp-12, 0x1.77342a958ba08p-69, -0x1.b5781e9d7c647p-11, 0x1.89e17c074d38p-10, -0x1.ed4ac7daea105p-10, 0x1.c11f27064a03p-10, -0x1.2add1ce973472p-10, 0x1.151f79a9bc54ap-11, -0x1.1c63bc7b26648p-13, -0x1.0f6f6d09ea759p-16, 0x1.168a1bdac929fp-15}, /* j=46 79.895 */
   {0x1.fffbb8f1049c6p-1, 0x1.d2c62526d166ep-56, 0x1.a740684026555p-13, -0x1.9004337b41fb4p-69, -0x1.36d34c8f1c26ap-11, 0x1.1eb6e14974a3p-10, -0x1.714eb8cc09935p-10, 0x1.5bec08c00afa9p-10, -0x1.e4621d7cfa292p-11, 0x1.e1b7b7b6e2bedp-12, -0x1.2456aedf774dcp-13, 0x1.b63a2dcc2de01p-20, 0x1.941679c46252p-16}, /* j=47 79.545 */
   {0x1.fffd1ac4135f9p-1, 0x1.eeafa1faf5581p-55, 0x1.2408e9ba3327fp-13, -0x1.822e0b4049783p-67, -0x1.b60d5e974cbbep-12, 0x1.9db74b1d1dcf2p-11, -0x1.11c85b1e8fef4p-10, 0x1.0a7b5546a7434p-10, -0x1.82f235b0c5ffbp-11, 0x1.998b49bc4f6d2p-12, -0x1.1aa5db242158dp-13, 0x1.d1a967621f96p-17, 0x1.05c06da0d736p-16}, /* j=48 81.222 */
   {0x1.fffe0e0140857p-1, -0x1.6aa36e5a4dff4p-57, 0x1.8fdc1b2dcf7b9p-14, 0x1.3e5ce22c9f818p-71, -0x1.322484cf12daap-12, 0x1.27dc1bc6cff05p-11, -0x1.9202f465eafd4p-11, 0x1.93b4c97452985p-11, -0x1.30e9e616d72a3p-11, 0x1.555b9ee9faf49p-12, -0x1.055957c332f5bp-13, 0x1.6884616670e08p-16, 0x1.295370cefbe29p-17}, /* j=49 80.539 */
   {0x1.fffeb3ebb267bp-1, 0x1.e47f68c40f19cp-57, 0x1.0f9e1b4dd36dfp-14, -0x1.63e9331d998bp-71, -0x1.a8670aa99a5bcp-13, 0x1.a3737e2a2f2d3p-12, -0x1.24544f02d2d04p-11, 0x1.2e7e763d237c8p-11, -0x1.da496e50284eep-12, 0x1.176cf788cada9p-12, -0x1.d2aca30d06e03p-14, 0x1.aab46b9073c9p-16, 0x1.096232db6a97ep-18}, /* j=50 81.209 */
   {0x1.ffff2436a21dcp-1, -0x1.3607957e07cfdp-55, 0x1.6e2367dc27f95p-15, 0x1.6e5abc89dadfp-73, -0x1.23c436c36fdabp-13, 0x1.26bf00867a845p-12, -0x1.a51fb50b15bfp-12, 0x1.c0825378e6d1fp-12, -0x1.6c3dbfe2c092bp-12, 0x1.c1dd15d6bac58p-13, -0x1.94c34ab048a09p-14, 0x1.bed7130ef46b1p-16, 0x1.146b5564fe5fp-23}, /* j=51 81.766 */
   {0x1.ffff6f9f67e55p-1, 0x1.e1e448293b293p-55, 0x1.e9b5e8d00ce77p-16, -0x1.d3b907526e62ep-70, -0x1.8de3cd290a7cp-14, 0x1.9aa489e3cad35p-13, -0x1.2c7d5ef0540dcp-12, 0x1.490a4d23005f8p-12, -0x1.145464e8c87f6p-12, 0x1.647f7322dc914p-13, -0x1.567491b40fea2p-14, 0x1.b2a229b6ffc4fp-16, -0x1.33a4de0b1c1a9p-19}, /* j=52 82.405 */
   {0x1.ffffa1de8c582p-1, 0x1.832540293fac5p-55, 0x1.44f21e49054f2p-16, 0x1.e9ee218e4f664p-71, -0x1.0d18811478659p-14, 0x1.1b964d438f62ep-13, -0x1.a8d7851f266c3p-13, 0x1.ddd6df9b581eap-13, -0x1.9e52b7adef1dfp-13, 0x1.165b20c9027a2p-13, -0x1.1b75a9125f3fcp-14, 0x1.918fd666ca55ap-16, -0x1.0454fe44d87cdp-18}, /* j=53 82.397 */
   {0x1.ffffc316d9edp-1, -0x1.8b32f4429e0a9p-55, 0x1.abe09e9144b5ep-17, 0x1.2cfda42d58866p-71, -0x1.690585ca91f98p-15, 0x1.84522fe8815c5p-14, -0x1.298f8d45f622ap-13, 0x1.57757788590d5p-13, -0x1.330aab7996d18p-13, 0x1.ac9998548e49ep-14, -0x1.cc155c0462b26p-15, 0x1.64bdaaf6a6563p-16, -0x1.362a507697cdep-18}, /* j=54 83.289 */
   {0x1.ffffd8e1a2f22p-1, -0x1.c10adf68fa3e4p-55, 0x1.1783ceac2891p-17, -0x1.7eccec27a8863p-71, -0x1.e06a8b37e5b93p-16, 0x1.07978c7b8496bp-14, -0x1.9d039884f8afbp-14, 0x1.e8d1145e95209p-14, -0x1.c1f72511fa379p-14, 0x1.458b9e0508b87p-14, -0x1.6eb051981283ep-15, 0x1.330486317fa6bp-16, -0x1.42e1d73962dc9p-18}, /* j=55 85.9 */
   {0x1.ffffe710d565ep-1, 0x1.c9ea52d257bf2p-55, 0x1.6a597219a93dap-18, -0x1.c6b2bd58b33c2p-72, -0x1.3d0e43d67415ep-16, 0x1.62ccea63cb0b8p-15, -0x1.1c07721ac824bp-14, 0x1.586bafc9c2954p-14, -0x1.46153fb80580bp-14, 0x1.e827fa55ddbe4p-15, -0x1.1f6311fd1dbb6p-15, 0x1.013a48887f99ep-16, -0x1.36c930de7a7c1p-18}, /* j=56 84.504 */
   {0x1.fffff039f9e8fp-1, -0x1.9d1bcd5f68bfep-55, 0x1.d21397ead99cbp-19, -0x1.2a0d2ed05c128p-75, -0x1.9f19734d29cf9p-17, 0x1.d982bd41d892cp-16, -0x1.8320fc4836a0bp-15, 0x1.e0a1cb1d23356p-15, -0x1.d38422316850bp-15, 0x1.696dae63c387ap-15, -0x1.bb6e241768923p-16, 0x1.a50cb2bcb0836p-17, -0x1.1cc4094f8f8fep-18}, /* j=57 84.916 */
   {0x1.fffff618c3da6p-1, -0x1.19309cf55a478p-58, 0x1.296a70f414053p-19, 0x1.293c1f6062c78p-74, -0x1.0d88765d3224bp-17, 0x1.394b1fa67113ep-16, -0x1.05760ad1bd1d4p-15, 0x1.4c1fe48a7b57cp-15, -0x1.4b9820321f1a9p-15, 0x1.085c0e5800969p-15, -0x1.510a47ddec286p-16, 0x1.5177ee6971c44p-17, -0x1.f2872584c5847p-19}, /* j=58 84.703 */
   {0x1.fffff9d446cccp-1, -0x1.bb06babeac674p-57, 0x1.789fb715aae95p-20, 0x1.0956a09d29bep-77, -0x1.5b333cc7f98f1p-18, 0x1.9b12fdbf90f04p-17, -0x1.5e06923144fcbp-16, 0x1.c6a071929864p-16, -0x1.d178cb02136bp-16, 0x1.7e29d0dd1cddfp-16, -0x1.f92040553e854p-17, 0x1.09609e09cef46p-17, -0x1.a72a310ec7611p-19}, /* j=59 84.914 */
   {0x1.fffffc2f171e3p-1, 0x1.85edd0397c65p-55, 0x1.d9371e2ff7c35p-21, 0x1.5202c8f3298fbp-75, -0x1.bba3ac4cf8472p-19, 0x1.0b6a7b0f1b52dp-17, -0x1.d06f586093e29p-17, 0x1.3436bca002e59p-16, -0x1.4357b5549aae9p-16, 0x1.110de2606b03bp-16, -0x1.7566bd2eb15a3p-17, 0x1.9a198592bfeebp-18, -0x1.5d78cf21a4ba5p-19}, /* j=60 85.274 */
   {0x1.fffffda86faa9p-1, -0x1.d230252f2108ap-56, 0x1.26f9df8519bd7p-21, -0x1.a2e22f04c8321p-75, -0x1.1926290adc888p-19, 0x1.5900c02d97264p-18, -0x1.3166de6a8c956p-17, 0x1.9dfcc328e2d3ap-17, -0x1.bcab1ed40dd03p-17, 0x1.81cd70a45b74ap-17, -0x1.106e9d2a6220fp-17, 0x1.37b62a7c86c6ep-18, -0x1.194b07a12c0f3p-19}, /* j=61 85.21 */
   {0x1.fffffe92ced93p-1, -0x1.d2db2ed0002aep-55, 0x1.6ce1aa3fd7bddp-22, 0x1.13ed125d83be2p-77, -0x1.617a9cedd8ffep-20, 0x1.b95fa39b39e6cp-19, -0x1.8e1fc41538c87p-18, 0x1.13717218617b7p-17, -0x1.2eb290b15a9c8p-17, 0x1.0d8c325af8e02p-17, -0x1.88856aa17d9e3p-18, 0x1.d2934641015a7p-19, -0x1.bbc8c51978a49p-20}, /* j=62 85.63 */
   {0x1.ffffff233ee1dp-1, 0x1.db123ed19a6b7p-55, 0x1.bfd7555a3bd68p-23, 0x1.b8cc44d611cc5p-77, -0x1.b8d7f804d2e73p-21, 0x1.17f93e51491a4p-19, -0x1.013b0457d079bp-18, 0x1.6b245d7ebdbabp-18, -0x1.98077549bcd27p-18, 0x1.7491fed2828eap-18, -0x1.1750678fa9c9dp-18, 0x1.5817080d0e657p-19, -0x1.5736752eec0f6p-20}, /* j=63 85.884 */
   {0x1.ffffff7b91176p-1, 0x1.0b2865615e8p-56, 0x1.10b1488aeb235p-23, 0x1.dc1e199a7e4ap-79, -0x1.10b1488aeb235p-21, 0x1.603a5308c4f68p-20, -0x1.4980e25286c45p-19, 0x1.da5f10dd5691p-19, -0x1.105053770992ep-18, 0x1.fd7c5c2e8dd46p-19, -0x1.88c7abe64b422p-19, 0x1.f46e5cf7c229bp-20, -0x1.048ded88153cfp-20}, /* j=64 86.243 */
   {0x1.ffffffb127525p-1, 0x1.504f382dcd18ap-55, 0x1.4980cb3c8094ap-24, -0x1.96b1ba28bde55p-78, -0x1.4ea6ce697296fp-22, 0x1.b771d9b6f057p-21, -0x1.a26c653fad22bp-20, 0x1.3302bb8a04p-19, -0x1.67f42e53ab378p-19, 0x1.58b4a663135b5p-19, -0x1.10f56fb67b135p-19, 0x1.6704a8f8cd5afp-20, -0x1.852f31b66594fp-21}, /* j=65 86.536 */
   {0x1.ffffffd169d0cp-1, 0x1.70a2bfb074fbfp-55, 0x1.8b0cfce0579ep-25, -0x1.0023f8020d62p-80, -0x1.976564c75a5afp-23, 0x1.0fdac559b6d9bp-21, -0x1.07600ca6e0cabp-20, 0x1.89ca77472035bp-20, -0x1.d73aa4e62a3fep-20, 0x1.cd9e024805ec9p-20, -0x1.7710ad753d5f9p-20, 0x1.fc74497f0572p-21, -0x1.1e02399a89503p-21}, /* j=66 86.933 */
   {0x1.ffffffe4aed5ep-1, 0x1.389c0f32ef38p-59, 0x1.d5f3a8dea7358p-26, 0x1.c596eee9b353p-84, -0x1.ebfb14c9170cp-24, 0x1.4d9228525f19ep-22, -0x1.48b536addaa88p-21, 0x1.f48ccf2584ef8p-21, -0x1.3183b6141b8ep-20, 0x1.31efd596e5414p-20, -0x1.fd9ee08120945p-21, 0x1.63859edff2482p-21, -0x1.9e1b73c607077p-22}, /* j=67 87.356 */
   {0x1.fffffff01a8b6p-1, 0x1.23370ec9c328p-60, 0x1.155a09065d4f7p-26, 0x1.b1d47e175a544p-81, -0x1.26afa996c3246p-24, 0x1.95ea6fdffafb5p-23, -0x1.96ba7366c00dep-22, 0x1.3b46801b2143dp-21, -0x1.8868e1d29db17p-21, 0x1.916e8578d3f2fp-21, -0x1.566efe5102ed3p-21, 0x1.eb1aa01b58ffp-22, -0x1.277d5c8e286cdp-22}, /* j=68 87.751 */
   {0x1.fffffff6d1e56p-1, -0x1.64d969b4b984fp-55, 0x1.44d26de513198p-27, -0x1.ad8796c5e0e48p-84, -0x1.5e32de7af8978p-25, 0x1.e9e05b3c8edbap-24, -0x1.f2f6fa7db5486p-23, 0x1.899dcad0514ccp-22, -0x1.f34b7ef1db60cp-22, 0x1.04bdf9b51bc2dp-21, -0x1.c73bc28a82144p-22, 0x1.4f28119de4e21p-22, -0x1.a01fb8ae4e202p-23}, /* j=69 88.214 */
   {0x1.fffffffabd229p-1, -0x1.4dbe49bebc458p-57, 0x1.7974e743dea3dp-28, 0x1.a3cf40636ddbap-82, -0x1.9cd7dcf23b833p-26, 0x1.252af6f48bd3ep-24, -0x1.2f7354e6b687cp-23, 0x1.e7102f8b9917ep-23, -0x1.3ab0b0f21d5d6p-22, 0x1.4f60fae1aaf26p-22, -0x1.2b6315fb8aba9p-22, 0x1.c4229b3190406p-23, -0x1.210962848ed0ap-23}, /* j=70 88.705 */
   {0x1.fffffffd01f89p-1, -0x1.35e8e39882c1p-56, 0x1.b334fac4b9f9ap-29, 0x1.f1e97646ba55p-83, -0x1.e2cec6323e50fp-27, 0x1.5c027d5bb9d7ap-25, -0x1.6df4d024ffb49p-24, 0x1.2aaf7c226f95ep-23, -0x1.8902edfe055cdp-23, 0x1.ab2a9eb3a266fp-23, -0x1.85abd35b6e223p-23, 0x1.2d7ef483471c5p-23, -0x1.8c6a4a0e8ee46p-24}, /* j=71 89.206 */
   {0x1.fffffffe4fa3p-1, 0x1.d166bcb68683p-57, 0x1.f1e3523b41d7ep-30, 0x1.31524e7486ad3p-84, -0x1.180fde4155097p-27, 0x1.99b866561854cp-26, -0x1.b598cb4614605p-25, 0x1.6b1baf48524cap-24, -0x1.e650e348908d9p-24, 0x1.0d67823ddc90cp-23, -0x1.f5f305ae078bbp-24, 0x1.8d989f091de07p-24, -0x1.0c7bb77556ff5p-24}, /* j=72 89.714 */
   {0x1.ffffffff0dd2bp-1, 0x1.0df73e7d2fc7ap-55, 0x1.1a94ff5716551p-30, -0x1.6ff991f6f55b6p-84, -0x1.4251f33f5578fp-28, 0x1.de6bc1f75b029p-27, -0x1.036b5fd1c3edcp-25, 0x1.b58f1389e06bfp-25, -0x1.2a2347f12222cp-24, 0x1.508da617933c1p-24, -0x1.3ffe9dd458648p-24, 0x1.0348c2f59174cp-24, -0x1.67367b756fda5p-25}, /* j=73 90.25 */
   {0x1.ffffffff79626p-1, 0x1.5fbc52d650d66p-55, 0x1.3e44e45301b94p-31, -0x1.55d19cb6914f6p-86, -0x1.6fffa7fff9fe2p-29, 0x1.1508f768ead8cp-27, -0x1.30fd0c66a5622p-26, 0x1.055632869acp-25, -0x1.6a3a9d4ce10fep-25, 0x1.a06fb101b8373p-25, -0x1.93e25717d6d67p-25, 0x1.4e84ce9f0393fp-25, -0x1.db01034af2d16p-26}, /* j=74 90.814 */
   {0x1.ffffffffb5be5p-1, -0x1.729d6819c777ep-56, 0x1.63ac6b4edc89p-32, -0x1.373c8505dd9dp-90, -0x1.a0ce0dc06a708p-30, 0x1.3e380dd75893p-28, -0x1.638bc4fb0244p-27, 0x1.35753ad86e3c7p-26, -0x1.b41f33ceeb9e5p-26, 0x1.fe692cc3ebe24p-26, -0x1.f8aee6307d10ap-26, 0x1.aafcee0f4818ap-26, -0x1.366d57fdfe3cdp-26}, /* j=75 91.392 */
   {0x1.ffffffffd759dp-1, 0x1.f7bee7eb236b3p-55, 0x1.8a61745ec7d2p-33, -0x1.afe847deb72cp-91, -0x1.d453ba308d495p-31, 0x1.6a8aeba4757d2p-29, -0x1.9b017abbf09bep-28, 0x1.6b43c957dca81p-27, -0x1.042f2a68ece85p-26, 0x1.35dc6cec5fefcp-26, -0x1.3834d2895d0bep-26, 0x1.0da5a535bd341p-26, -0x1.911fe34462b89p-27}, /* j=76 91.985 */
   {0x1.ffffffffe9ebp-1, -0x1.ea527e0bee1dp-58, 0x1.b1e5acf351d8bp-34, -0x1.02b715d324efep-88, -0x1.05042a0a5f3c5p-31, 0x1.99ac8fd63b444p-30, -0x1.d72344378d096p-29, 0x1.a6be9a188ee78p-28, -0x1.33aacb4fb9192p-27, 0x1.74b715dbf2c49p-27, -0x1.7e7e920190a02p-27, 0x1.5109ca2c4657cp-27, -0x1.00470c3149729p-27}, /* j=77 92.593 */
   {0x1.fffffffff4188p-1, 0x1.7a2cb3d057019p-55, 0x1.d9a880f306bddp-35, -0x1.42ddf1d0c3027p-89, -0x1.20a2ae94181bap-32, 0x1.cb2a2e56403b6p-31, -0x1.0bc6ecf663a8cp-29, 0x1.e7ba576f04ffep-29, -0x1.689347704de3bp-28, 0x1.bc2a9e250f2e4p-28, -0x1.d00fdac39a30dp-28, 0x1.a0f4af0ff7dbbp-28, -0x1.43e129fdd063p-28}, /* j=78 93.219 */
   {0x1.fffffffff9a1bp-1, -0x1.6a87270d2428cp-57, 0x1.0084ff12563ap-35, -0x1.7f6624101d286p-90, -0x1.3ca42adaa26f8p-33, 0x1.fe73513c65deap-32, -0x1.2dd9aa5a2b25bp-30, 0x1.16ef6b98d6be2p-29, -0x1.a2d58ea11e479p-29, 0x1.0638838c8ffbap-28, -0x1.16cdc4f97d577p-28, 0x1.fe98aee84af83p-29, -0x1.94e03729bfe15p-29}, /* j=79 93.87 */
   {0x1.fffffffffc9e8p-1, -0x1.a759f7738927ep-56, 0x1.13af4f04f999cp-36, -0x1.d18475e09366fp-90, -0x1.589b22c638001p-34, 0x1.196da0aa68475p-32, -0x1.516d3cb76b1e3p-31, 0x1.3c51d0b14aca2p-30, -0x1.e23586e9bbe31p-30, 0x1.32c70aa088adep-29, -0x1.4bce8954c5bdfp-29, 0x1.357f910b2d259p-29, -0x1.f4c469ee27ebep-30}, /* j=80 94.532 */
   {0x1.fffffffffe38p-1, 0x1.7ce07114e501ep-55, 0x1.25f9ee0b923ep-37, 0x1.ea30ac3beee08p-91, -0x1.74105146a5166p-35, 0x1.33cde4f35c17fp-33, -0x1.760fe7b664f3cp-32, 0x1.63a70fdebdd97p-31, -0x1.1324f70052c8fp-30, 0x1.63a2f434fa172p-30, -0x1.8724a81206a9dp-30, 0x1.737a8ffb42b1bp-30, -0x1.326d03d8dc12cp-30}, /* j=81 95.208 */
   {0x1.ffffffffff11ap-1, -0x1.3eafccbc6e86ap-56, 0x1.370ab8327af64p-38, -0x1.9cfeae61fc22ap-92, -0x1.8e85bc00ad8b5p-36, 0x1.4decacbf852bap-34, -0x1.9b3c558006cf2p-33, 0x1.8c78e45259811p-32, -0x1.373cd6ea2344ap-31, 0x1.988a82c73b1a2p-31, -0x1.c8bfec71a33aep-31, 0x1.b97bbc2fe3fb6p-31, -0x1.732047d271b68p-31}, /* j=82 95.902 */
   {0x1.ffffffffff845p-1, 0x1.b0edc5a89abc8p-56, 0x1.46897d4b69fcep-39, -0x1.9ab5c2183b349p-93, -0x1.a77a4e7dcd73cp-37, 0x1.67543695da825p-35, -0x1.c05c1e2fc4f3dp-34, 0x1.b63941ac80ab5p-33, -0x1.5cfd7ec1be14fp-32, 0x1.d1153f1d46c19p-32, -0x1.082f82d6965bap-31, 0x1.03c7b8aadb259p-31, -0x1.bce2427697b96p-32}, /* j=83 96.613 */
   {0x1.ffffffffffc05p-1, 0x1.07ba96a6b2e2ap-55, 0x1.5422ef5d894a6p-40, -0x1.7424d9ca2e79ep-95, -0x1.be6dda2ac4316p-38, 0x1.7f8a0f3e20c1dp-36, -0x1.e4cb4aea6e8ap-35, 0x1.e044b3fa50853p-34, -0x1.83ea4babfcc2ep-33, 0x1.065935eb2c55fp-32, -0x1.2ec52cbf65a0ep-32, 0x1.2ec365d7b096ep-32, -0x1.07f00574bf6a1p-32}, /* j=84 97.337 */
   {0x1.ffffffffffdf8p-1, -0x1.dcf8b10ff973p-55, 0x1.5f8b87a31bd8fp-41, 0x1.cc132ce2e535fp-95, -0x1.d2e55024a0fbfp-39, 0x1.9612cc225ab44p-37, -0x1.03ee5f38b816p-35, 0x1.04f2f727492aap-34, -0x1.ab709a05d04cbp-34, 0x1.25548f59435c1p-33, -0x1.57c84aebf692fp-33, 0x1.5d77eb7f1a36p-33, -0x1.3609391e967dep-33}, /* j=85 98.078 */
   {0x1.ffffffffffef8p-1, 0x1.14be6226402d4p-56, 0x1.68823e52970cap-42, 0x1.95ff0771e009ap-96, -0x1.e46f03befaf8bp-40, 0x1.aa76120eaf294p-38, -0x1.146faeb88e85bp-36, 0x1.192d3b30a367fp-35, -0x1.d2eaae774c204p-35, 0x1.450d19c860284p-34, -0x1.82c2cf245b684p-34, 0x1.8f87e2cd02dbfp-34, -0x1.688c54c5042e3p-34}, /* j=86 98.833 */
   {0x1.fffffffffff7bp-1, 0x1.00fa07f7fb616p-55, 0x1.6ed2f2515e942p-43, -0x1.9fa86a9d62d32p-97, -0x1.f2a6c1669c91p-41, 0x1.bc42ba389cbbcp-39, -0x1.2391e135ad606p-37, 0x1.2c6c2461b350bp-36, -0x1.f9a3c1c24301ap-36, 0x1.65021aa3d0017p-35, -0x1.af21f4efb136ap-35, 0x1.c46fe92122461p-35, -0x1.9f2a55ab17d94p-35}, /* j=87 99.605 */
   {0x1.fffffffffffbep-1, -0x1.182b326b228d9p-55, 0x1.7258610b3b244p-44, -0x1.1f51e682a35efp-98, -0x1.fd39856f71518p-42, 0x1.cb12e2f5e6bcp-40, -0x1.31011e96bd173p-38, 0x1.3e4a1f981164p-37, -0x1.0f6e89d7edafp-36, 0x1.84a49df3142e4p-36, -0x1.dc3874482d5cp-36, 0x1.fb86adff48dcfp-36, -0x1.d9603855a25e9p-36}, /* j=88 100.387 */
   {0x1.fffffffffffdfp-1, 0x1.5669e670f914ep-56, 0x1.72fd93e036cefp-45, 0x1.900335ce573a8p-99, -0x1.01f450d1e61bbp-42, 0x1.d68fb81b28d48p-41, -0x1.3c706aa4cf311p-39, 0x1.4e647967300ddp-38, -0x1.20e9eb8f55b02p-37, 0x1.a35b503354ddp-37, -0x1.04a10c00c45c8p-36, 0x1.19ff67a4222d2p-36, -0x1.0b3f6dda2ba9fp-36}, /* j=89 101.187 */
   {0x1.ffffffffffffp-1, -0x1.20ef3618f2d52p-56, 0x1.70beaf9c7ffccp-46, -0x1.4fe56c6797e82p-101, -0x1.0346137a09fd8p-43, 0x1.de74c0dc33e3p-42, -0x1.459c8175bfc1fp-40, 0x1.5c5ee415c685ep-39, -0x1.30e3dc1e26e3cp-38, 0x1.c08722ceeb646p-38, -0x1.1ab45df5c3b08p-37, 0x1.36754f4f66492p-37, -0x1.2adb0c00b7e32p-37}, /* j=90 102.004 */
   {0x1.ffffffffffff8p-1, 0x1.0160ef15c497ep-56, 0x1.6ba91ac73479ep-47, 0x1.2c7efca8f4866p-102, -0x1.028a39099f4e5p-44, 0x1.e292863e0fe2ep-43, -0x1.4c4e690fba1e8p-41, 0x1.67e6e5c1ef392p-40, -0x1.3f00d819db8b5p-39, 0x1.db888bdda3829p-39, -0x1.2fe554dcabcb8p-38, 0x1.52a0e7d1ba0ecp-38, -0x1.4afcbe1054ad8p-38}, /* j=91 102.829 */
   {0x1.ffffffffffffcp-1, 0x1.8115fd1b12786p-56, 0x1.63daf8b4b1e28p-48, -0x1.eab97d93de2c9p-102, -0x1.ff8ac583bfb45p-46, 0x1.e2d06d6fcb6cfp-44, -0x1.505d95359ddadp-42, 0x1.70b7012e4e8d3p-41, -0x1.4aed68394188ap-40, 0x1.f3c524ce0ae47p-40, -0x1.43c0e92ad7198p-39, 0x1.6df475870ed29p-39, -0x1.6b1462bc0942p-39}, /* j=92 103.655 */
   {0x1.ffffffffffffep-1, 0x1.59ab24e589a3p-56, 0x1.5982008db1322p-49, 0x1.c1e4160369b8p-107, -0x1.f610e8cde57c6p-47, 0x1.df2dac2f23cfap-45, -0x1.51b17f95f733ep-43, 0x1.76996df71ac71p-42, -0x1.546155bcbf635p-41, 0x1.0456b0f2374e8p-40, -0x1.55d5e5098124dp-40, 0x1.87dd23632a8d6p-40, -0x1.8a9213148f2adp-40}, /* j=93 104.529 */
   {0x1.fffffffffffffp-1, 0x1.0439397b5f70ap-56, 0x1.4cd9c04158cf8p-50, 0x1.cf1236b38982p-104, -0x1.e8dfd25ffa70ep-48, 0x1.d7c149fc93adbp-46, -0x1.50429df37d26fp-44, 0x1.796a3a21a279fp-43, -0x1.5b228193b2a7bp-42, 0x1.0ce10b9168a33p-41, -0x1.65b92db24da5ep-41, 0x1.9fc80e506523ap-41, -0x1.a8c0fb19b0709p-41}, /* j=94 105.385 */
};

/* Add a + b, such that *hi + *lo approximates a + b.
   Assumes |a| >= |b|.
   For rounding to nearest we have hi + lo = a + b exactly.
   For directed rounding, we have
   (a) hi + lo = a + b exactly when the exponent difference between a and b
       is at most 53 (the binary64 precision)
   (b) otherwise |(a+b)-(hi+lo)| <= 2^-105 min(|a+b|,|hi|)
       (see https://hal.inria.fr/hal-03798376)
   We also have |lo| < ulp(hi). */
static inline void fast_two_sum(double *hi, double *lo, double a, double b) {
  double e;

  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
}

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma (a, b, -*hi);
}

/* Assuming 0.0625 0 <= z <= 0x1.7afb48dc96626p+2, put in h+l an approximation
   of erf(z) with relative error bounded by 2^-66.8 (cf analyze_p() in
   erf.sage). */
static double
cr_erf_fast (double *h, double *l, double z)
{
  double th, tl;
  if (z < 0.0625)
  {
    /* the following is a degree-11 minimax polynomial for erf(x) on [0,1/16]
       generated by Sollya, with double-double coefficients for degree 1 and 3,
       and double coefficients for degrees 5 to 11 */
    static const double c0[] = {
      0x1.20dd750429b6dp+0, 0x1.1ae3a7862d9c4p-56,  /* degree 1 */
      -0x1.812746b0379e7p-2, 0x1.f1a64d72722a2p-57. /* degree 3 */
      0x1.ce2f21a042b7fp-4,                         /* degree 5 */
      -0x1.b82ce31189904p-6,                        /* degree 7 */
      0x1.565bbf8a0fe0bp-8,                         /* degree 9 */
      -0x1.bf9f8d2c202e4p-11};                      /* degree 11 */
    double z2 = z*z, z4 = z2*z2;
    double c9 = __builtin_fma (c0[7], z2, c0[6]);
    double c5 = __builtin_fma (c0[5], z2, c0[4]);
    c5 = __builtin_fma (c9, z4, c5);
    /* compute c0[2] + c0[3] + z2*c5 */
    a_mul (&th, &tl, z2, c5);
    fast_two_sum (h, l, c0[2], th);
    l += tl + c0[3];
    /* compute c0[0] + c0[1] + z2*(h + l) */
    a_mul (&th, &tl, z2, h);
    tl = __builtin_fma (z2, l, tl); /* tl += z2*l */
    fast_two_sum (h, l, c0[0], th);
    l += tl + c0[1];
    /* multiply (h,l) by z */
    a_mul (h, &tl, h, z);
    *l = __builtin_fma (*l, z, tl);
    return;
  }
  double v = __builtin_floor (16.0 * z);
  uint32_t i = 16.0 * z;
  /* i/16 <= z < (i+1)/16 */
  /* For 0.0625 0 <= z <= 0x1.7afb48dc96626p+2, z - 0.03125 is exact:
     (1) either z - 0.03125 is in the same binade as z, then 0.03125 is
         an integer multiple of ulp(z), so is z - 0.03125
     (2) if z - 0.03125 is in a smaller binade, both z and 0.03125 are
         integer multiple of the ulp() of that smaller binade.
     Also, subtracting 0.0625 * v is exact. */
  z = (z - 0.03125) - 0.0625 * v;
  /* now |z| <= 1/32 */
  const double *c = C[i-1];
  double z2 = z*z, z4 = z2*z2;
  /* the degree-10 coefficient is c[12] */
  double c9 = __builtin_fma (c[12], z, c[11]);
  double c7 = __builtin_fma (c[10], z, c[9]);
  double c5 = __builtin_fma (c[8], z, c[7]);
  double c3 = __builtin_fma (c[6], z, c[5]);
  c7 = __builtin_fma (c9, z2, c7);
  c3 = __builtin_fma (c5, z2, c3);
  c3 = __builtin_fma (c7, z4, c3);
  /* the degree-2 coefficient is c[4] */
  double c2 = __builtin_fma (c3, z, c[4]);
  /* compute c[2] + c[3] + z*c2 */
  a_mul (&th, &tl, z, c2);
  fast_two_sum (h, l, c[2], th);
  l += tl + c[3];
  /* compute c[0] + c[1] + z*(h + l) */
  a_mul (&th, &tl, z, h);
  tl = __builtin_fma (z, l, tl); /* tl += z*l */
  fast_two_sum (h, l, c[0], th);
  l += tl + c[1];
}

double
cr_erf (double x)
{
  double z = __builtin_fabs (x);
  b64u64_u t = {.f = z};
  uint64_t ux = t.u;
  /* erf(x) rounds to +/-1 for RNDN for |x| > 0x1.7afb48dc96626p+2 */
  if (ux > 0x4017afb48dc96626) /* 0x4017afb48dc96626 == 0x1.7afb48dc96626p+2 */
  {
    double os = __builtin_copysign (1.0, x);
#define MASK (uint64_t) 0x7f80000000000000
    if (ux > MASK)
      return x; /* NaN */
    if (ux == MASK)
      return os; /* +/-Inf */
    return os - 0x1p-54 * os;
  }
  /* now |x| <= 0x1.7afb48dc96626p+2 */
  double h, l;
  cr_erf_fast (&h, &l, z);
}
