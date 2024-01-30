/* Correctly rounded log(1+x) for binary64 values.

Copyright (c) 2024 Alexei Sibidanov.

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

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef uint64_t u64;
typedef int64_t i64;
typedef unsigned short ushort;
typedef union {double f; uint64_t u;} b64u64_u;

static inline double fasttwosum(double x, double y, double *e){
  double s = x + y, z = s - x;
  *e = y - z;
  return s;
}

static inline double twosum(double xh, double ch, double *l){
  double s = xh + ch, d = s - xh;
  *l = (ch - d) + (xh + (d - s));
  return s;
}

static inline double fastsum(double xh, double xl, double yh, double yl, double *e){
  double sl, sh = fasttwosum(xh, yh, &sl);
  *e = (xl + yl) + sl;
  return sh;
}

static inline double sum(double xh, double xl, double ch, double cl, double *l){
  double s = xh + ch, d = s - xh;
  *l = ((ch - d) + (xh + (d - s))) + (xl + cl);
  return s;
}

static inline double muldd(double xh, double xl, double ch, double cl, double *l){
  double ahhh = ch*xh;
  *l = (cl*xh + ch*xl) + __builtin_fma(ch, xh, -ahhh);
  return ahhh;
}

static inline double mulddd(double x, double ch, double cl, double *l){
  double ahhh = ch*x;
  *l = cl*x + __builtin_fma(ch, x, -ahhh);
  return ahhh;
}

static inline double polydd(double xh, double xl, int n, const double c[][2], double *l){
  int i = n-1;
  double ch = fasttwosum(c[i][0], *l, l), cl = c[i][1] + *l;
  while(--i>=0){
    ch = muldd(xh,xl, ch,cl, &cl);
    ch = fastsum(c[i][0],c[i][1], ch,cl, &cl);
  }
  *l = cl;
  return ch;
}

static inline double polyddd(double x, int n, const double c[][2], double *l){
  int i = n-1;
  double ch = fasttwosum(c[i][0], *l, l), cl = c[i][1] + *l;
  while(--i>=0){
    ch = mulddd(x, ch,cl, &cl);
    ch = fastsum(c[i][0],c[i][1], ch,cl, &cl);
  }
  *l = cl;
  return ch;
}

static double __attribute__((noinline)) as_log1p_refine(double, double);

static const double rt[4][16] = {
  {0x1.6a09e68p+0, 0x1.5ab07ep+0, 0x1.4bfdad8p+0, 0x1.3dea65p+0,
   0x1.306fe08p+0, 0x1.2387a7p+0, 0x1.172b84p+0, 0x1.0b5587p+0,
   0x1p+0, 0x1.ea4afap-1, 0x1.d5818ep-1, 0x1.c199bep-1,
   0x1.ae89f98p-1, 0x1.9c4918p-1, 0x1.8ace54p-1, 0x1.7a1147p-1},
  {0x1.059b0dp+0, 0x1.04e5f7p+0, 0x1.04315e8p+0, 0x1.037d43p+0,
   0x1.02c9a4p+0, 0x1.0216818p+0, 0x1.0163da8p+0, 0x1.00b1af8p+0,
   0x1p+0, 0x1.fe9d968p-1, 0x1.fd3c228p-1, 0x1.fbdba38p-1,
   0x1.fa7c18p-1, 0x1.f91d8p-1, 0x1.f7bfdbp-1, 0x1.f663278p-1},
  {0x1.0058c88p+0, 0x1.004dad8p+0, 0x1.0042938p+0, 0x1.0037798p+0,
   0x1.002c608p+0, 0x1.0021478p+0, 0x1.00162fp+0, 0x1.000b178p+0,
   0x1p+0, 0x1.ffe9d2p-1, 0x1.ffd3a58p-1, 0x1.ffbd798p-1,
   0x1.ffa74e8p-1, 0x1.ff91248p-1, 0x1.ff7afb8p-1, 0x1.ff64d38p-1},
  {0x1.00058b8p+0, 0x1.0004dap+0, 0x1.0004288p+0, 0x1.0003778p+0,
   0x1.0002c6p+0,  0x1.0002148p+0, 0x1.000163p+0, 0x1.0000b18p+0,
   0x1p+0, 0x1.fffe9dp-1, 0x1.fffd3ap-1, 0x1.fffbd78p-1,
   0x1.fffa748p-1, 0x1.fff9118p-1, 0x1.fff7ae8p-1, 0x1.fff64cp-1}
};

static const double ln[4][16][3] = {
  {{0x1.61fa45b636ea2p-95, -0x1.497697cb3134p-45, -0x1.62e43033a8p-2},
   {-0x1.42b671e5e6d5ep-93, 0x1.fc50ee605bd8p-44, -0x1.3687aa721cp-2},
   {0x1.0cc2ea968d329p-95, 0x1.f48666c24eacp-43, -0x1.0a2b247d58p-2},
   {-0x1.963904db0a34ep-93, -0x1.58495e50d36ap-43, -0x1.bb9d3d808p-3},
   {0x1.337e0d9f81652p-95, 0x1.c89a0bede978p-45, -0x1.62e42f0378p-3},
   {-0x1.1bdf53871d50bp-95, 0x1.2d67e8a0c1d7p-42, -0x1.0a2b24a1p-3},
   {-0x1.781dd6fc4053p-98, 0x1.07a22a9e81ed8p-42, -0x1.62e43327cp-4},
   {-0x1.84e4562b8f2f3p-96, -0x1.868afe14895p-43, -0x1.62e435baep-5},
   {0x0p+0, 0x0p+0, 0x0p+0},
   {-0x1.17e3ec05cde7p-96, 0x1.174a19689569p-42, 0x1.62e432b22p-5},
   {-0x1.8067ec2a27738p-95, 0x1.111a4eadf313p-43, 0x1.62e42e4a8p-4},
   {-0x1.60bef06b14405p-93, 0x1.dd4ec4e1d421p-43, 0x1.0a2b233f08p-3},
   {-0x1.b18e160362c24p-94, 0x1.6bd65e8b0b7p-45, 0x1.62e43056cp-3},
   {-0x1.c6ac3f1862a6bp-93, -0x1.0f5aad509ea8p-44, 0x1.bb9d3cbd68p-3},
   {-0x1.dead1a4581acfp-93, -0x1.ac842de00564p-43, 0x1.0a2b244da4p-2},
   {-0x1.96b1f2f60e3e4p-93, 0x1.9cec9a50db228p-42, 0x1.3687aa9b78p-2}},
  {{0x1.dbd93cdc08613p-96, 0x1.96f39c49fa8ap-44, -0x1.62e423dd4p-6},
   {-0x1.3c008f058a015p-94, -0x1.f26429af65b7p-43, -0x1.36879e514p-6},
   {0x1.f96a86136c1d1p-93, -0x1.bfc5afb2bd378p-42, -0x1.0a2b2240cp-6},
   {0x1.3917cea885bep-93, -0x1.421e0c297638p-46, -0x1.bb9d4b288p-7},
   {0x1.f0f1b5d989675p-93, 0x1.92a7506c90a18p-42, -0x1.62e43c12p-7},
   {-0x1.f7f8d6a190051p-93, -0x1.94925ae0e524p-42, -0x1.0a2b41ddp-7},
   {0x1.ccb0919263869p-94, -0x1.0fb337ff2cf28p-42, -0x1.62e41068p-8},
   {0x1.c9a5efe4b7c6cp-93, -0x1.966c27788d1cp-44, -0x1.62e3e4ccp-9},
   {0x0p+0, 0x0p+0, 0x0p+0},
   {-0x1.40288ccae8f0fp-95, -0x1.f3c5ff88d19c8p-42, 0x1.62e462b6p-9},
   {-0x1.f7e788a87135p-94, -0x1.a9610028771cp-43, 0x1.62e44c93p-8},
   {0x1.574cc6d3f577dp-93, 0x1.2a3a1a65aa398p-42, 0x1.0a2b1e33p-7},
   {-0x1.56bb79b254f33p-99, -0x1.4a995b6d9ddcp-44, 0x1.62e4367cp-7},
   {0x1.c4c209ca6783dp-94, -0x1.d7b98ef45911p-43, 0x1.bb9d449a8p-7},
   {0x1.e9731de7f0155p-93, -0x1.de722390bbd6p-43, 0x1.0a2b1f194p-6},
   {0x1.fadc62522444dp-96, -0x1.fd95cb835e38p-45, 0x1.3687ad114p-6}},
  {{-0x1.1f72d2a6a460ep-95, 0x1.a25045c37c33p-43, -0x1.62e4795p-10},
   {0x1.2847e318fd3ffp-93, -0x1.c633239e1cdc8p-42, -0x1.3686e5dcp-10},
   {0x1.63f5750b9d826p-93, -0x1.cb15590f1cd78p-42, -0x1.0a2b6538p-10},
   {-0x1.1382d6395c24cp-94, 0x1.7636a5400e2bp-43, -0x1.bb9bf138p-11},
   {0x1.be4491ec20322p-94, 0x1.d0c5e2c9b6be8p-42, -0x1.62e53e5p-11},
   {-0x1.706f04cc2c9f6p-93, -0x1.f3bc0ce9b9a08p-42, -0x1.0a2ab37p-11},
   {0x1.ab6f41df01d7fp-93, 0x1.a7cbc9a97ba4p-44, -0x1.62e0ap-12},
   {-0x1.67e487663ca3ep-100, -0x1.957976dc5f36p-43, -0x1.62e84fcp-13},
   {0x0p+0, 0x0p+0, 0x0p+0},     
   {0x1.9681e48dde135p-93, -0x1.868625640a69p-43, 0x1.62e7bp-13},
   {0x1.a2948cd558655p-93, -0x1.2ee3d96b696ap-42, 0x1.62e35f6p-12},
   {-0x1.cfc26ccf6d0e4p-96, 0x1.53edbcf1165p-46, 0x1.0a2b4b2p-11},
   {0x1.f68d24b9e338dp-93, 0x1.783e334613p-51, 0x1.62e4be1p-11},
   {-0x1.f33369bf7dff1p-95, -0x1.60785f20acb2p-42, 0x1.bb9e085p-11},
   {-0x1.685a35575eff1p-95, -0x1.5a62ec66568p-48, 0x1.0a2b94d4p-10},
   {-0x1.4c4d1abca79bfp-95, 0x1.7ded26dc813p-46, 0x1.368810f8p-10}},
  {{0x1.40676dea39b19p-94, 0x1.25337681fa9p-42, -0x1.62dc284p-14},
   {-0x1.ec76c964ac8d6p-94, -0x1.504032e2ed388p-42, -0x1.367d0ecp-14},
   {0x1.8ecaeb6d9523cp-95, 0x1.c28fde83047ep-44, -0x1.0a1dd6cp-14},
   {-0x1.e4069214576bep-93, 0x1.58f2757976ef8p-42, -0x1.bbbcffp-15},
   {0x1.86146d1f9b91fp-93, -0x1.d4717ca0a323p-42, -0x1.62fe138p-15},
   {0x1.1eb44b07102c5p-93, -0x1.86ffcda25278p-43, -0x1.0a3eebp-15},
   {-0x1.95bbb5fd5a466p-93, 0x1.1ce399729e7bp-43, -0x1.62ff0ap-16},
   {-0x1.79059023f8767p-93, -0x1.dbf1c6a400408p-42, -0x1.62ff84p-17},
   {0x0p+0, 0x0p+0, 0x0p+0},
   {-0x1.b4cb08d6fd48fp-93, -0x1.db0e38e5aaa98p-42, 0x1.63007cp-17},
   {-0x1.0e6a88bfc3838p-93, 0x1.2b1c75580439p-43, 0x1.6300f6p-16},
   {0x1.8074feacfe49dp-94, 0x1.401bb919f14ep-42, 0x1.0a21148p-15},
   {-0x1.f56f5f168db15p-93, -0x1.85d6f6487ce2p-44, 0x1.62e1ecp-15},
   {0x1.9eae7e05a0143p-93, -0x1.af5d58a7c9218p-42, 0x1.bba301p-15},
   {-0x1.a859095999ae3p-94, 0x1.590faa0883bd8p-42, 0x1.0a32298p-14},
   {0x1.4f787e495e5ep-93, -0x1.fbb791220a18p-46, 0x1.3682f14p-14}}
};

double cr_log1p(double x){
  static const struct {ushort c0; short c1;} B[] = {
    {301, 27565}, {7189, 24786}, {13383, 22167}, {18923, 19696},
    {23845, 17361}, {28184, 15150}, {31969, 13054}, {35231, 11064},
    {37996, 9173}, {40288, 7372}, {42129, 5657}, {43542, 4020},
    {44546, 2457}, {45160, 962}, {45399, -468}, {45281, -1838},
    {44821, -3151}, {44032, -4412}, {42929, -5622}, {41522, -6786},
    {39825, -7905}, {37848, -8982}, {35602, -10020}, {33097, -11020},
    {30341, -11985}, {27345, -12916}, {24115, -13816}, {20661, -14685},
    {16989, -15526}, {13107, -16339}, {9022, -17126}, {4740, -17889}
  };

  b64u64_u ix = {.f = x};
  u64 ax = ix.u<<1;
  double ln1, ln0, eps;
  if(__builtin_expect(ax<0x7f60000000000000ul, 1)){
    double x2 = x*x;
    if(__builtin_expect(ax<0x7e60000000000000ul, 1)){
      ln1 = x;
      eps = 0x1.6p-64*x;
      if(__builtin_expect(ax<0x7d43360000000000ul, 1)){
	if(ax<0x7940000000000000ul){
	  if(!ax) return x;
	  return __builtin_fma(__builtin_fabs(x), -0x1p-54, x);
	}
	static const double c[] = {-0x1.00000000001d1p-1, 0x1.55555555558f7p-2};
	ln0 = x2*(c[0] + x*c[1]);
      } else {
	static const double c[] =
	  {-0x1.ffffffffffffdp-2, 0x1.5555555555551p-2, -0x1.000000d5555e1p-2, 0x1.99999b442f73fp-3};
	ln0 = x2*((c[0]+x*c[1])+x2*(c[2]+x*c[3]));
      }
    } else {
      static const double c[] =
	{0x1.5555555555555p-2, -0x1p-2, 0x1.9999999999b41p-3, -0x1.555555555583bp-3,
	 0x1.24924923f39ep-3, -0x1.fffffffe42e43p-4, 0x1.c71c75511d70bp-4, -0x1.99999de10510fp-4,
	 0x1.7457e81b175f6p-4, -0x1.554fb43e54e0fp-4, 0x1.3ed68744f3d18p-4, -0x1.28558ad5a7ac4p-4};
      double x3 = x2*x, x4 = x2*x2, hx = -0.5*x;
      ln1 = __builtin_fma(hx,x,x);
      ln0 = __builtin_fma(hx,x,x-ln1);
      double f = ((c[0]+x*c[1])+x2*(c[2]+x*c[3])) +
	x4*(((c[4]+x*c[5])+x2*(c[6]+x*c[7])) + x4*((c[8]+x*c[9])+x2*(c[10]+x*c[11])));
      ln0 += x3*f;
      eps = x3*0x1.94p-52;
    }
  } else {
    static const double c[] =
      {-0x1p-1, 0x1.5555555555555p-2, -0x1.fffffffff2p-3, 0x1.99999999887e6p-3,
       -0x1.5555c55570f25p-3, 0x1.2492b9217dbfep-3};
    b64u64_u t, dt;
    if(__builtin_expect((i64)ix.u<0x4340000000000000l && ix.u<0xbff0000000000000ul, 1)){
      t.f = fasttwosum(1.0, x, &dt.f);
    } else {
      if(__builtin_expect(ix.u<0x4690000000000000ul, 1)){
	t.f = x; dt.f = 1;
      } else {
	if(__builtin_expect(ix.u<0x7ff0000000000000ul, 1)){
	  t.f = x; dt.f = 0;
	} else {
	  if((ix.u<<1)>0xffe0000000000000ul) return x; // nan
	  if(ix.u==0x7ff0000000000000ul) return x; // +inf
	  if(ix.u==0xbff0000000000000ul) return -1./0.0; // -1
	  return 0.0/0.0; // <-1
	}
      }
    }  
    u64 i = (t.u>>(52-5))&0x1f;
    long d = t.u & (~0ul>>17);
    i64 j = t.u + 13841960651842387968ul + ((u64)B[i].c0<<33) + ((long)B[i].c1*(d>>16));
    i64 j1 = (j>>48)&0xf, j2 = (j>>44)&0xf, je = j>>52;
    i64 eoff = je<<52;
    t.u -= eoff;
    if(__builtin_expect(dt.u<<1, 1)) dt.u -= eoff;
    double r = rt[0][j1]*rt[1][j2];
    double dh = r*t.f, dl = __builtin_fma(r,t.f,-dh) + r*dt.f;
    double xl, xh = fasttwosum(dh-1.0, dl, &xl), x2 = xh*xh;
    xl += x2*((c[0] + xh*c[1]) + x2*((c[2] + xh*c[3]) + x2*(c[4] + xh*c[5])));
    double L1 = 0x1.62e42fefa4p-1*je, L0 = -0x1.8432a1b0e2634p-43*je;
    ln1 = ln[0][j1][2] + ln[1][j2][2] + L1;
    ln0 = ln[0][j1][1] + ln[1][j2][1] + L0;
    ln1 = fastsum(ln1, ln0, xh, xl, &ln0);
    eps = 0x1p-68;
  }
  double lb = ln1 + (ln0 - eps), ub = ln1 + (ln0 + eps);
  if(__builtin_expect(lb != ub, 0)) return as_log1p_refine(x, lb);
  return lb;
}

static double __attribute__((noinline)) as_log1p_refine(double x, double a){
  static const double cz[][2] = {
    {0x1.5555555555555p-2, 0x1.5555555555556p-56}, {-0x1p-2, 0x1.25558eff3c1efp-86},
    {0x1.999999999999ap-3, -0x1.999999a91d6cap-57}, {-0x1.5555555555555p-3, -0x1.588aab185593ap-57},
    {0x1.2492492492492p-3, 0x1.28c490abc528fp-57}
  };
  static const double czl[] = {
    -0x1.fffffffffc555p-4, 0x1.c71c71c7185aap-4, -0x1.9999d44449a31p-4, 0x1.745d51f1817c7p-4
  };
  static const double cy[][2] = {
    {1,0}, {-0.5, 0},
    {0x1.5555555555555p-2, 0x1.5555555555555p-56}, {-0x1p-2, -0x1.80007ed9858b7p-107},
    {0x1.999999999999ap-3, -0x1.9999999999991p-57}
  };
  static const double cl[4] = {
    -0x1.5555555555555p-3, 0x1.2492492492492p-3, -0x1.0000000073334p-3, 0x1.c71c71c802f68p-4
  };
  b64u64_u ix = {.f = x};
  double ln22, ln21, ln20;
  u64 ax = ix.u<<1;
  if(ax<0x7ea0000000000000ul){
    if(ax<0x7940000000000000ul){
      if(!ax) return x;
      return __builtin_fma(__builtin_fabs(x), -0x1p-54, x);
    }
    double x2h = x*x, x2l = __builtin_fma(x,x,-x2h);
    double x3l, x3h = mulddd(x,x2h, x2l, &x3l);
    double sl = x*((czl[0] + x*czl[1]) + x2h*(czl[2] + x*czl[3]));
    double sh = polyddd(x, 5,cz, &sl);
    sh = muldd(sh,sl,x3h,x3l, &sl);
    x2h *= -0.5;
    x2l *= -0.5;
    ln22 = x;
    ln21 = fastsum(x2h,x2l, sh,sl, &ln20);
  } else {
    i64 j = __builtin_roundeven(a*0x1.71547652b82fep+16);
    i64 i = j + 34952;
    int j1 = (i>>12)&0xf, j2 = (i>>8)&0xf, j3 = (i>>4)&0xf, j4 = i&0xf, je = i>>16;
    double L[3];
    L[0] = (ln[0][j1][0] + ln[1][j2][0]) + (ln[2][j3][0] + ln[3][j4][0]);
    L[1] = (ln[0][j1][1] + ln[1][j2][1]) + (ln[2][j3][1] + ln[3][j4][1]);
    L[2] = (ln[0][j1][2] + ln[1][j2][2]) + (ln[2][j3][2] + ln[3][j4][2]);

    b64u64_u t, dt;
    if((i64)ix.u<0x4690000000000000l && ix.u<0xbfe0000000000000ul ){
      t.f = twosum(1.0, x, &dt.f);
      if(__builtin_expect(!(dt.u<<1), 0)) dt.u = 0;
    } else {
      if(__builtin_expect((i64)ix.u>=0x4690000000000000l, 0))
	t.f = x;
      else 
	t.f = 1 + x;
      dt.u = 0;
    }
    t.u -= (long)je<<52;

    double t12 = rt[0][j1]*rt[1][j2], t34 = rt[2][j3]*rt[3][j4];
    double th = t12*t34, tl = __builtin_fma(t12,t34,-th);
    double dh = th*t.f, dl = __builtin_fma(th,t.f,-dh);
    double sh = tl*t.f, sl = __builtin_fma(tl,t.f,-sh);
    double xl, xh = fasttwosum(dh-1, dl, &xl);
    xh = fastsum(xh, xl, sh, sl, &xl);
    if(dt.u){
      dt.u -= (long)je<<52;
      double ddh = th*dt.f, ddl = __builtin_fma(th,dt.f,-ddh) + tl*dt.f;
      xh = fastsum(xh, xl, ddh, ddl, &xl);
    }
    sl = xh*((cl[0] + xh*cl[1]) + (xh*xh)*(cl[2] + xh*cl[3]));
    sh = polydd(xh,xl, 5,cy, &sl);
    sh = muldd(xh,xl, sh,sl, &sl);
    if(je){
      ln22 = 0x1.62e42fefa4p-1*je; ln21 = -0x1.8432a1b0e28p-43*je; ln20 = 0x1.cc01f97b57a08p-87*je;
      ln22 += L[2];
      ln21 = fastsum(ln21, ln20, L[1], L[0], &ln20);
    } else {
      ln22 = L[2];
      ln21 = L[1];
      ln20 = L[0];
    }
    ln21 = sum(ln21, ln20, sh, sl, &ln20);
  }
  ln22 = fasttwosum(ln22,ln21, &ln21);
  ln21 = fasttwosum(ln21,ln20, &ln20);
  
  b64u64_u t = {.f = ln21};
  if(__builtin_expect(!(t.u&(~0ul>>12)), 0)){
    b64u64_u w = {.f = ln20};
    if((w.u^t.u)>>63)
      t.u--;
    else
      t.u++;
    ln21 = t.f;
  }
  return ln22 + ln21;
}
