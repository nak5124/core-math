/* Correctly-rounded biased argument base-10 logarithm function for binary32 value.

Copyright (c) 2022 Alexei Sibidanov.

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

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;

float cr_log10p1f(float x) {
  static const double ix[] = {
    0x1p+0, 0x1.fc07f01fcp-1, 0x1.f81f81f82p-1, 0x1.f44659e4ap-1, 0x1.f07c1f07cp-1,
    0x1.ecc07b302p-1, 0x1.e9131abfp-1, 0x1.e573ac902p-1, 0x1.e1e1e1e1ep-1,
    0x1.de5d6e3f8p-1, 0x1.dae6076bap-1, 0x1.d77b654b8p-1, 0x1.d41d41d42p-1,
    0x1.d0cb58f6ep-1, 0x1.cd8568904p-1, 0x1.ca4b3055ep-1, 0x1.c71c71c72p-1,
    0x1.c3f8f01c4p-1, 0x1.c0e070382p-1, 0x1.bdd2b8994p-1, 0x1.bacf914c2p-1,
    0x1.b7d6c3ddap-1, 0x1.b4e81b4e8p-1, 0x1.b2036406cp-1, 0x1.af286bca2p-1,
    0x1.ac5701ac6p-1, 0x1.a98ef606ap-1, 0x1.a6d01a6dp-1, 0x1.a41a41a42p-1,
    0x1.a16d3f97ap-1, 0x1.9ec8e951p-1, 0x1.9c2d14ee4p-1, 0x1.99999999ap-1,
    0x1.970e4f80cp-1, 0x1.948b0fcd6p-1, 0x1.920fb49dp-1, 0x1.8f9c18f9cp-1,
    0x1.8d3018d3p-1, 0x1.8acb90f6cp-1, 0x1.886e5f0acp-1, 0x1.861861862p-1,
    0x1.83c977ab2p-1, 0x1.818181818p-1, 0x1.7f405fd02p-1, 0x1.7d05f417ep-1,
    0x1.7ad2208ep-1, 0x1.78a4c8178p-1, 0x1.767dce434p-1, 0x1.745d1745ep-1,
    0x1.724287f46p-1, 0x1.702e05c0cp-1, 0x1.6e1f76b44p-1, 0x1.6c16c16c2p-1,
    0x1.6a13cd154p-1, 0x1.681681682p-1, 0x1.661ec6a52p-1, 0x1.642c8590cp-1,
    0x1.623fa7702p-1, 0x1.605816058p-1, 0x1.5e75bb8dp-1, 0x1.5c9882b94p-1,
    0x1.5ac056b02p-1, 0x1.58ed23082p-1, 0x1.571ed3c5p-1, 0x1.555555556p-1,
    0x1.5390948f4p-1, 0x1.51d07eae2p-1, 0x1.501501502p-1, 0x1.4e5e0a73p-1,
    0x1.4cab88726p-1, 0x1.4afd6a052p-1, 0x1.49539e3b2p-1, 0x1.47ae147aep-1,
    0x1.460cbc7f6p-1, 0x1.446f86562p-1, 0x1.42d6625d6p-1, 0x1.414141414p-1,
    0x1.3fb013fbp-1, 0x1.3e22cbce4p-1, 0x1.3c995a47cp-1, 0x1.3b13b13b2p-1,
    0x1.3991c2c18p-1, 0x1.381381382p-1, 0x1.3698df3dep-1, 0x1.3521cfb2cp-1,
    0x1.33ae45b58p-1, 0x1.323e34a2cp-1, 0x1.30d19013p-1, 0x1.2f684bda2p-1,
    0x1.2e025c04cp-1, 0x1.2c9fb4d82p-1, 0x1.2b404ad02p-1, 0x1.29e4129e4p-1,
    0x1.288b01288p-1, 0x1.27350b882p-1, 0x1.25e22708p-1, 0x1.24924924ap-1,
    0x1.23456789ap-1, 0x1.21fb78122p-1, 0x1.20b470c68p-1, 0x1.1f7047dc2p-1,
    0x1.1e2ef3b4p-1, 0x1.1cf06ada2p-1, 0x1.1bb4a4046p-1, 0x1.1a7b9611ap-1,
    0x1.19453808cp-1, 0x1.181181182p-1, 0x1.16e068942p-1, 0x1.15b1e5f76p-1,
    0x1.1485f0e0ap-1, 0x1.135c81136p-1, 0x1.12358e75ep-1, 0x1.111111112p-1,
    0x1.0fef010fep-1, 0x1.0ecf56be6p-1, 0x1.0db20a89p-1, 0x1.0c9714fbcp-1,
    0x1.0b7e6ec26p-1, 0x1.0a6810a68p-1, 0x1.0953f3902p-1, 0x1.084210842p-1,
    0x1.073260a48p-1, 0x1.0624dd2f2p-1, 0x1.05197f7d8p-1, 0x1.041041042p-1,
    0x1.03091b52p-1, 0x1.020408102p-1, 0x1.01010101p-1, 0x1p-1};

  static const double lix[] = {
    0x0p+0, -0x1.bafd472256696p-9, -0x1.b9476a4fb1457p-8, -0x1.49b08514b989p-7,
    -0x1.b5e908eb4b0ffp-7, -0x1.10a83a8435685p-6, -0x1.45f4f5ad5f82fp-6, -0x1.7adc3df399ad7p-6,
    -0x1.af5f92b02a2c7p-6, -0x1.e3806acc4d6cap-6, -0x1.0ba01a8134efap-5, -0x1.25502c0fd7ed5p-5,
    -0x1.3ed1199a4969bp-5, -0x1.58238eeb914fbp-5, -0x1.71483427a2057p-5, -0x1.8a3fadebf1f07p-5,
    -0x1.a30a9d6083332p-5, -0x1.bba9a058dc5edp-5, -0x1.d41d5164db865p-5, -0x1.ec6647eb5bf9fp-5,
    -0x1.02428c1ef6a23p-4, -0x1.0e3d29d81e6d4p-4, -0x1.1a23445508744p-4, -0x1.25f5215ed6964p-4,
    -0x1.31b3055c3238ep-4, -0x1.3d5d335c2dbe2p-4, -0x1.48f3ed1e0e9e8p-4, -0x1.547773197adb2p-4,
    -0x1.5fe8048896fdcp-4, -0x1.6b45df6f5226dp-4, -0x1.769140a260558p-4, -0x1.81ca63d085b29p-4,
    -0x1.8cf1838848b55p-4, -0x1.9806d9417c816p-4, -0x1.a30a9d60df453p-4, -0x1.adfd0741aa6a4p-4,
    -0x1.b8de4d3abacc6p-4, -0x1.c3aea4a5cde2dp-4, -0x1.ce6e41e46060ep-4, -0x1.d91d586694046p-4,
    -0x1.e3bc1ab0bee19p-4, -0x1.ee4aba6145d8ep-4, -0x1.f8c968346f0bep-4, -0x1.019c2a06389bbp-3,
    -0x1.06cbd67a47c05p-3, -0x1.0bf3d0939f002p-3, -0x1.11142f0829878p-3, -0x1.162d082ae3016p-3,
    -0x1.1b3e71ec72395p-3, -0x1.204881df09e84p-3, -0x1.254b4d35d4b7ep-3, -0x1.2a46e8ca5d3c1p-3,
    -0x1.2f3b691c436abp-3, -0x1.3428e253eadfep-3, 0x1.2f7301cf65fb6p-3, 0x1.2a935ba61396cp-3,
    0x1.25ba8215d23e2p-3, 0x1.20e8624051cp-3, 0x1.1c1ce9955892fp-3, 0x1.175805d15502ap-3,
    0x1.1299a4fb64682p-3, 0x1.0de1b563712e3p-3, 0x1.093025a1b463ep-3, 0x1.0484e4941945p-3,
    0x1.ffbfc2bbff173p-4, 0x1.f68216c9c73c4p-4, 0x1.ed50a4a21d0ap-4, 0x1.e42b4c1704e13p-4,
    0x1.db11ed76be21bp-4, 0x1.d204698ccf18fp-4, 0x1.c902a19e24ca9p-4, 0x1.c00c7766dc01ap-4,
    0x1.b721cd170e8b5p-4, 0x1.ae428550a9d18p-4, 0x1.a56e8325ab159p-4, 0x1.9ca5aa1777424p-4,
    0x1.93e7de0fbcf52p-4, 0x1.8b350364bf329p-4, 0x1.828cfed25f10ep-4, 0x1.79efb57b2d087p-4,
    0x1.715d0ce3ad2c7p-4, 0x1.68d4eaf240544p-4, 0x1.605735eedf988p-4, 0x1.57e3d47c385cap-4,
    0x1.4f7aad9bed5fp-4, 0x1.471ba8a7f67d8p-4, 0x1.3ec6ad545e626p-4, 0x1.367ba3aa55189p-4,
    0x1.2e3a740bc4708p-4, 0x1.2603072a54df8p-4, 0x1.1dd5460cd9433p-4, 0x1.15b11a0999254p-4,
    0x1.0d966cc6491ccp-4, 0x1.0585283722143p-4, 0x1.fafa6d3a22492p-5, 0x1.eafd0502eca5cp-5,
    0x1.db11ed7711842p-5, 0x1.cb38fcccfbce3p-5, 0x1.bb7209d1e5c7cp-5, 0x1.abbcebd8806e1p-5,
    0x1.9c197abfae953p-5, 0x1.8c878eeb3c9e4p-5, 0x1.7d070145903bdp-5, 0x1.6d97ab3aec004p-5,
    0x1.5e3966b787e12p-5, 0x1.4eec0e23d865ep-5, 0x1.3faf7c66e884fp-5, 0x1.30838cdbd3adcp-5,
    0x1.21681b5d39d8fp-5, 0x1.125d04324518bp-5, 0x1.0362241e90b97p-5, 0x1.e8eeb0a0806fep-6,
    0x1.cb38fccf111e7p-6, 0x1.ada2e8e25cacdp-6, 0x1.902c31d5296ap-6, 0x1.72d4956dd3c4bp-6,
    0x1.559bd23f02e66p-6, 0x1.3881a7b8bc456p-6, 0x1.1b85d60432cf6p-6, 0x1.fd503c3c3fd76p-7,
    0x1.c3d083774d2d1p-7, 0x1.8a8c06bb38facp-7, 0x1.51824c76abc39p-7, 0x1.18b2dc8fe2981p-7,
    0x1.c03a80b5194dp-8, 0x1.4f820527a8cep-8, 0x1.be76bd777d654p-9, 0x1.bd96a1d61d141p-10,
    0x0p+0,
  };

  double z = x;
  b32u32_u t = {.f = x};
  unsigned ux = t.u, ax = ux&(~0u>>1);
  if (__builtin_expect(ux >= 0x17fu<<23, 0)) { // x <= -1
    if (ux==(0x17fu<<23)) return -1.0/0.0f; // -1.0
    if (ux>(0x1ffu<<23)) return x; // nan
    return __builtin_nanf("-"); // x < -1
  } else if(__builtin_expect(ax >= (0xff<<23), 0)){ // +inf, nan
    if(ax > (0xff<<23)) return x; // nan
    return __builtin_inff();
  } else if(__builtin_expect(ax<0x3cb7aa26u, 1)){ // 0x1.6f544cp-6
    double z2 = z*z, z4 = z2*z2;
    if(__builtin_expect(ax<0x3b9d9d34u, 1)){ // 0x1.3b3a68p-8
      if(__builtin_expect(ax<0x39638a7eu, 1)){ // 0x1.c714fcp-13
	if(__builtin_expect(ax<0x329c5639u, 1)){ // 0x1.38ac72p-26
	  if(__builtin_expect(ux == 0xa6aba8afu, 0)) return -0x1.2a33bcp-51f + 0x1p-104f;
	  if(__builtin_expect(ux == 0xaf39b9a7u, 0)) return -0x1.42a342p-34f + 0x1p-86f;
	  static const double c[] =
	    {0x1.bcb7b1526e50fp-2, -0x1.bcb7b1526e50fp-3};
	  return z*(c[0] + z*c[1]);
	} else {
	  static const double c[] =
	    {0x1.bcb7b1526e50ep-2, -0x1.bcb7b1526e50dp-3, 0x1.287a76c37d1edp-3, -0x1.bcb7b26b76b1p-4};
	  return z*((c[0] + z*c[1]) + z2*(c[2] + z*c[3]));
	}
      } else {
	if(__builtin_expect(ux == 0x399a7c00u, 0)) return 0x1.0c53cap-13f + 0x1.fp-38f;
	static const double c[] =
	  {0x1.bcb7b1526e50fp-2, -0x1.bcb7b1526e511p-3, 0x1.287a76364fe3ep-3, -0x1.bcb7b150c5959p-4,
	   0x1.63c86971690bp-4, -0x1.287cb84cbc498p-4};
	return z*((c[0] + z*c[1]) + z2*(c[2] + z*c[3]) + z4*((c[4] + z*c[5])));
      }
    } else {
      static const double c[] =
	{0x1.bcb7b1526e50ep-2, -0x1.bcb7b1526e50bp-3, 0x1.287a763701b2fp-3, -0x1.bcb7b152997b4p-4,
	 0x1.63c6235d07b0ep-4, -0x1.287a71622092bp-4, 0x1.fca571759497p-5, -0x1.bd1d7daf1ddadp-5};
      return z*((c[0] + z*c[1]) + z2*(c[2] + z*c[3]) + z4*((c[4] + z*c[5]) + z2*(c[6] + z*c[7])));
    }
  } else {
    if(__builtin_expect(ux == 0x7956ba5eu, 0)) return 0x1.16bebap+5f + 0x1.cp-53f;
    int32_t xe = ax>>23, s = xe - 119;
    if(__builtin_expect(s>=11 && s<32 && (ux<<s) == 0x80000000, 0)){
      int32_t xm = (ux&~0u>>9)|1<<23, k = xm>>(31-s);
      if(k == 9) return 1.0f;
      if(k == 99) return 2.0f;
      if(k == 999) return 3.0f;
      if(k == 9999) return 4.0f;
      if(k == 99999) return 5.0f;
      if(k == 999999) return 6.0f;
      if(k == 9999999) return 7.0f;
    }
    b64u64_u t = {.f = z + 1.0};
    uint64_t m = t.u&(~0ul>>12);
    int32_t e = (t.u>>52) - 0x3ff;
    int32_t j = (m + (1l<<(52-8)))>>(52-7), k = j>53;
    e += k;
    b64u64_u xd = {.u = m | 0x3fful<<52};
    z = __builtin_fma(xd.f, ix[j], -1.0);
    static const double c[] =
      {0x1.bcb7b1526e50ep-2, -0x1.bcb7b1526e50fp-3, 0x1.287a7636acbcap-3, -0x1.bcb7b151b5928p-4,
       0x1.63c7a4a6483a9p-4, -0x1.287bf378d551ep-4};
    double z2 = z*z;
    double c0 = c[0] + z*c[1];
    double c2 = c[2] + z*c[3];
    double c4 = c[4] + z*c[5];
    c0 += z2*(c2 + z2*c4);
    b64u64_u res = {.f = (z*c0 - lix[j]) + e*0x1.34413509f79ffp-2};
    float r = res.f;
    return r;
  }
}

/* just to compile since glibc does not contain this function*/
float log10p1f(float x){
  return cr_log10p1f(x);
}
