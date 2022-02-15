/* Correctly rounded exponential function for binary64 values.

Copyright (c) 2021-2022 Paul Zimmermann, Inria.

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

#include <math.h>
#include <fenv.h>

static int rnd = FE_TONEAREST;

/* Implements the algorithm from "Algorithms for Calculating Correctly
   Rounded Exponential Function in Double-Precision Arithmetic", Alexander
   Godunov, IEEE Transactions on Computers, vol. 69, no. 9, September 2020,
   pages 1388-1400 (+ appendix of 10 pages).

   Warning: it is correctly rounded only for rounding to nearest. For example
   for x=-0x1.62c9fc0dcbdb5p+9 and rounding towards zero, it returns
   0x0.4e89e7c1b34c9p-1022, but the correct rounding is
   0x0.4e89e7c1b34cap-1022.
*/

/* Add a + b exactly, such that *hi + *lo = a + b.
   Assumes |a| >= |b| and rounding to nearest.  */
static void
fast_two_sum (double *hi, double *lo, double a, double b)
{
  double e;

  *hi = a + b;
  e = *hi - a; /* exact */
  *lo = b - e; /* exact */
  /* Now *hi + *lo = a + b exactly.  */
}

/* Algorithm 4.1. Case |x| < 2^(-36).  */
static double
algo41 (double x)
{
  double v, Y, Y1;
  {
    v = (x + 0.5) - 0.5;
    fast_two_sum (&Y, &Y1, 1.0, v);
  }

  /* rounding to nearest */
  if (rnd == FE_TONEAREST)
  {
    if (Y1 == 0)
      return Y;
    double t1 = ((x - v) + 0.5 * v * v);
    if (t1 > 0)
      return Y + (Y1 + 0x1p-55);
    else
      return Y + (Y1 - 0x1p-55);
  }
  
  /* directed roundings */
  if (Y1 == 0)
  {
    Y1 = (x - v) + 0.5 * v * v;
    if (Y1 == 0)
    {
      if (x > 0)
        Y1 = -0x1p-55;
      if (x < 0)
        Y1 = 0x1p-55;
    }
  }

  /* rounding to +Inf */
  if (rnd == FE_UPWARD)
  {
    if (Y1 > 0)
    {
      if (Y >= 1)
        Y = Y + 0x1p-52;
      else
        Y = Y + 0x1p-53;
    }
    return Y;
  }

  /* rounding to -Inf or zero */
  if (Y1 < 0)
  {
    if (Y > 1)
      Y = Y - 0x1p-52;
    else
      Y = Y - 0x1p-53;
  }
  return Y;
}

/* Algorithm 5.1. Case -2^(-28) < x < 2^(-27).  */
static double
algo51 (double x)
{
  double v, Y, Y1, t1, t2, v2;
  {
    v = (x + 0.5) - 0.5;
    fast_two_sum (&Y, &Y1, 1.0, v);
    v2 = v * v;
  }
  
  /* rounding to nearest */
  if (rnd == FE_TONEAREST)
  {
    if (Y == 0)
      return Y;
    t1 = ((x - v) + 0.5 * v2);
    t2 = - (v2 * v) / 3.0;
    if (t1 + t2 > 0)
      return Y + (Y1 + 0x1p-55);
    else
      return Y + (Y1 - 0x1p-55);
  }

  /* directed roundings */
  if (Y1 == 0)
  {
    t1 = ((x - v) + 0.5 * v2);
    t2 = - (v2 * v) / 3.0;
    Y1 = t1 + t2;
  }

  /* the end of the algorithm is the same as in Algorithm 4.1 */
  
  /* rounding to +Inf */
  if (rnd == FE_UPWARD)
  {
    if (Y1 > 0)
    {
      if (Y >= 1)
        Y = Y + 0x1p-52;
      else
        Y = Y + 0x1p-53;
    }
    return Y;
  }

  /* rounding to -Inf or zero */
  if (Y1 < 0)
  {
    if (Y > 1)
      Y = Y - 0x1p-52;
    else
      Y = Y - 0x1p-53;
  }
  return Y;
}

/* Return R_0(x).
   Warning: this assumes rounding to nearest, since for example with rounding
   toward +Inf, if x = 0.25, then C+x is rounded to C+1, and we get 1 instead
   of the expected value of 0. */
static inline double
R0 (double x)
{
  double C = 0x1.8p+52;
  return (C + x) - C;
}

/* return R_{-16}(x) */
static inline double
Rm16 (double x)
{
  double C = 0x1.8p+36;
  return (C + x) - C;
}

/* return R_{-34}(x) */
static inline double
Rm34 (double x)
{
  double C = 0x1.8p+18;
  return (C + x) - C;
}

/* return R_{-35}(x) */
static inline double
Rm35 (double x)
{
  double C = 0x1.8p+17;
  return (C + x) - C;
}

/* return R_{-66}(x) */
static inline double
Rm66 (double x)
{
  double C = 0x1.8p-14;
  return (C + x) - C;
}

static void
algo71 (double x, double *m0, double *x01, double *x02, double *x03)
{
  /* l1, l2, l3 approximate log(2) */
  /* Sage: l1=round(log(2)*2^42)/2^42 (l1 has 42 significant bits) */
  double l1 = 0x1.62e42fefa38p-1;
  /* Sage: l2=round((log(2)-l1)*2^85)/2^85 (l2 has 40 significant bits) */
  double l2 = 0x1.ef35793c76p-45;
  /* Sage: l3=round((log(2)-l1-l2)*2^138)/2^138 (50 significant bits) */
  double l3 = 0x1.cc01f97b57a08p-87;
  double l_tilde = 0x1.71547652b82fep+0;
  *m0 = R0 (x * l_tilde);
  *x01 = x - *m0 * l1;
  *x02 = -*m0 * l2;
  *x03 = -*m0 * l3;
}

static const double table1[45][4] = {
  {-19, -0x1.68ac83e9c6a14p-2, -0x1.a64p-58, -0x1.d5bae802f0ce4p-71},
  {-18, -0x1.522ae0738a3d8p-2, 0x1.8f8p-57, -0x1.64c759686a22p-73},
  {-17, -0x1.3c25277333184p-2, 0x1.2adp-56, 0x1.3f285476315ecp-71},
  {-16, -0x1.269621134db92p-2, -0x1.e0fp-56, 0x1.489893f5563ep-74},
  {-16, -0x1.269621134db92p-2, -0x1.e0fp-56, 0x1.489893f5563ep-74},
  {-15, -0x1.1178e8227e47cp-2, 0x1.0e6p-57, 0x1.d2f80e3485f88p-72},
  {-14, -0x1.f991c6cb3b37ap-3, 0x1.04dp-56, -0x1.419be6028636cp-71},
  {-13, -0x1.d1037f2655e7cp-3, 0x1.4fdp-56, -0x1.4921238d10f8p-72},
  {-13, -0x1.d1037f2655e7cp-3, 0x1.4fdp-56, -0x1.4921238d10f8p-72},
  {-12, -0x1.a93ed3c8ad9e4p-3, 0x1.21bp-56, -0x1.f53bd2e406e66p-70},
  {-11, -0x1.823c16551a3c2p-3, 0x1.124p-57, -0x1.a631e830fd30ap-70},
  {-10, -0x1.5bf406b543db2p-3, 0x1.2p-61, -0x1.49767e410316ep-70},
  {-9, -0x1.365fcb0159016p-3, -0x1.7d4p-58, -0x1.1a5b944aca88p-74},
  {-8, -0x1.1178e8227e47cp-3, 0x1.0e8p-58, -0x1.c5a0fe396f41p-70},
  {-8, -0x1.1178e8227e47cp-3, 0x1.0e8p-58, -0x1.c5a0fe396f41p-70},
  {-7, -0x1.da727638446a4p-4, 0x1.bp-56, -0x1.fa71733018beap-70},
  {-6, -0x1.9335e5d594988p-4, -0x1.5c4p-57, 0x1.50ae09996d7bcp-71},
  {-5, -0x1.4d3115d207eacp-4, -0x1.768p-58, -0x1.f42c7842cbd9ep-70},
  {-4, -0x1.08598b59e3a08p-4, 0x1.776p-56, -0x1.fecdfa819b96p-71},
  {-3, -0x1.894aa149fb34p-5, -0x1.9a9p-56, 0x1.05a267d770ceap-70},
  {-2, -0x1.0415d89e74448p-5, 0x1.c7fp-56, 0x1.18714564eeec4p-70},
  {-1, -0x1.020565893584p-6, -0x1.d28p-56, 0x1.b8bdf48c7088cp-71},
  {0, 0x0p3, 0x0p3, 0x0p3},
  {1, 0x1.fc0a8b0fc03ep-7, 0x1.e7cp-58, 0x1.eda74d37abd58p-71},
  {2, 0x1.f829b0e7833p-6, 0x1.34p-60, -0x1.c0fb0e10dd608p-72},
  {3, 0x1.77458f632ddp-5, -0x1.dcep-56, -0x1.61abc236b4fa8p-70},
  {4, 0x1.f0a30c01162a8p-5, -0x1.9e8p-57, -0x1.9b47488a66878p-72},
  {5, 0x1.341d7961bd1dp-4, 0x1.254p-57, -0x1.9f227becbb22cp-70},
  {6, 0x1.6f0d28ae56b4cp-4, -0x1.908p-58, 0x1.266e7b466d904p-70},
  {7, 0x1.a926d3a4ad564p-4, -0x1.35ep-57, -0x1.0b7558f156ce2p-70},
  {9, 0x1.0d77e7cd08e5ap-3, -0x1.32dp-56, -0x1.1d0b7e7aa2e38p-72},
  {10, 0x1.29552f81ff524p-3, -0x1.67fp-56, -0x1.11c77f0482cbcp-70},
  {11, 0x1.44d2b6ccb7d1ep-3, 0x1.9f4p-57, 0x1.eca87c3f0f062p-70},
  {12, 0x1.5ff3070a793d4p-3, -0x1.bc8p-58, 0x1.f105039091dd8p-70},
  {13, 0x1.7ab890210d90ap-3, -0x1.c84p-56, 0x1.b596b50304032p-70},
  {14, 0x1.9525a9cf456b4p-3, 0x1.d9p-57, 0x1.307538b896764p-71},
  {16, 0x1.c8ff7c79a9a22p-3, -0x1.4f6p-57, -0x1.13f08680232eep-70},
  {17, 0x1.e27076e2af2e6p-3, -0x1.618p-59, 0x1.43fff0ff4f0aap-70},
  {18, 0x1.fb9186d5e3e2ap-3, 0x1.1abp-56, -0x1.5cc9e435969bp-70},
  {19, 0x1.0a324e27390e3p-2, 0x1.7ddp-56, -0x1.0bfcf1fe78ecp-75},
  {21, 0x1.22941fbcf7966p-2, -0x1.76fp-56, -0x1.7ac258a2bcd0cp-70},
  {22, 0x1.2e8e2bae11d31p-2, -0x1.8f5p-56, 0x1.92350a103679cp-71},
  {23, 0x1.3a64c556945eap-2, -0x1.c68p-57, -0x1.946517e5ee414p-71},
  {25, 0x1.51aad872df82dp-2, 0x1.39p-59, 0x1.3d60cfaaf188ep-70},
  {26, 0x1.5d1bdbf5809cap-2, 0x1.423p-56, 0x1.8e0f71ff84568p-70}
};

static double
algo72 (double x01, double *x11, double *x12, double *x13)
{
  int n1 = R0 (x01 * 0x1p6);
  /* We should have -22 <= n1 <= 22.  */
  const double *p = table1[22 + n1];
  *x11 = x01 - p[1];
  *x12 = -p[2];
  *x13 = -p[3];
  return p[0]; /* m1 */
}

static const double table2[71][4] = {
  {-35, -0x1.1a6b91ac7338p-6, -0x1.875p-56, 0x1.2ac748a672ed4p-71},
  {-34, -0x1.12487a5507f7p-6, 0x1.ap-60, 0x1.2b4c6d7e52b38p-72},
  {-33, -0x1.0a266bb50869p-6, -0x1.be1p-56, -0x1.f0bb13ff3b29p-70},
  {-32, -0x1.020565893584p-6, -0x1.d28p-56, 0x1.b8bdf48c7088cp-71},
  {-31, -0x1.f3cacf1cd3b2p-7, -0x1.b08p-57, 0x1.b2cbc51419bap-71},
  {-30, -0x1.e38ce303331p-7, -0x1.776p-56, 0x1.ed140cf86c7bp-71},
  {-29, -0x1.d351063fa46ap-7, -0x1.b24p-58, -0x1.be0017d6367fp-70},
  {-28, -0x1.c317384c75fp-7, -0x1.98p-57, -0x1.8823013087f04p-71},
  {-27, -0x1.b2df78a428a6p-7, -0x1.92p-58, 0x1.68aeaf4333ec4p-70},
  {-26, -0x1.a2a9c6c17046p-7, -0x1.0ecp-58, -0x1.4965c6f72055cp-71},
  {-25, -0x1.9276221f3328p-7, -0x1.5afp-56, 0x1.6cbb5d815d7b6p-70},
  {-24, -0x1.82448a388a2ap-7, -0x1.441p-56, -0x1.62c26fe134004p-71},
  {-23, -0x1.7214fe88c094p-7, 0x1.b14p-57, 0x1.4663fbfc17de8p-71},
  {-22, -0x1.61e77e8b53fcp-7, -0x1.822p-57, -0x1.41c5432d00cecp-70},
  {-21, -0x1.51bc09bbf436p-7, 0x1.e14p-58, 0x1.ada75c3992c36p-70},
  {-20, -0x1.41929f96832ep-7, -0x1.f1dp-56, -0x1.5d013687bbab6p-70},
  {-19, -0x1.316b3f9714dcp-7, -0x1.b9ap-57, 0x1.d16aa70c6a982p-70},
  {-18, -0x1.2145e939ef1ep-7, -0x1.124p-56, 0x1.cec58013244ep-70},
  {-17, -0x1.11229bfb89a8p-7, -0x1.dc8p-57, -0x1.7f414ec3e7e5cp-70},
  {-16, -0x1.010157588de8p-7, 0x1.daep-56, 0x1.99d2be8312ff8p-70},
  {-15, -0x1.e1c4359baddp-8, 0x1.5f8p-58, -0x1.5b76c04295738p-70},
  {-14, -0x1.c189cbb0e28p-8, 0x1.544p-58, 0x1.2c42b02675e18p-71},
  {-13, -0x1.a1536feb35e4p-8, -0x1.c74p-56, -0x1.236d38492f742p-70},
  {-12, -0x1.8121214586b4p-8, -0x1.40ep-56, -0x1.4b9f9377a1d3p-73},
  {-11, -0x1.60f2debb161cp-8, 0x1.2fap-57, -0x1.612245d4d67e2p-70},
  {-10, -0x1.40c8a747879p-8, 0x1.e39p-56, -0x1.f070001892f1ep-70},
  {-9, -0x1.20a279e6e09cp-8, 0x1.935p-56, -0x1.244e4f5377494p-70},
  {-8, -0x1.0080559588b4p-8, 0x1.504p-57, -0x1.6638cf63676cep-70},
  {-7, -0x1.c0c472a092p-9, -0x1.57ep-57, 0x1.939c6a68bf2d8p-70},
  {-6, -0x1.80904828986p-9, 0x1.fcbp-56, 0x1.2b1e7e76024dep-70},
  {-5, -0x1.406429be3c7p-9, -0x1.e1cp-56, 0x1.3d91690c8e138p-71},
  {-4, -0x1.0040155d5888p-9, -0x1.de8p-57, 0x1.f31c2227e8062p-70},
  {-3, -0x1.80481205118p-10, -0x1.51ap-56, 0x1.8314f21a7dfecp-71},
  {-2, -0x1.00200556559p-10, 0x1.db3p-56, 0x1.506503c16f284p-71},
  {-1, -0x1.00100155756p-11, 0x1.ddcp-57, 0x1.10c7e47de444ap-70},
  {0, 0x0p3, 0x0p3, 0x0p3},
  {1, 0x1.ffe002aa6acp-12, -0x1.dep-57, 0x1.99e2b62cc632cp-70},
  {2, 0x1.ffc00aa8ab2p-11, -0x1.e08p-56, -0x1.fd97d736d9b66p-70},
  {3, 0x1.7fb811faf18p-10, 0x1.14ep-56, -0x1.86491276c1a3p-71},
  {4, 0x1.ff802a9ab11p-10, -0x1.988p-58, 0x1.4f1d0a9f1d8d4p-71},
  {5, 0x1.3f9c29972c68p-9, 0x1.974p-57, 0x1.ac83c856f4ed6p-70},
  {6, 0x1.7f7047d7984p-9, -0x1.2ccp-56, 0x1.3ad0c90273ff6p-70},
  {7, 0x1.bf3c720a81b8p-9, -0x1.a04p-56, -0x1.a184b73011478p-72},
  {8, 0x1.ff00aa2b10cp-9, -0x1.fdbp-56, 0x1.0d6ad369a96ep-74},
  {9, 0x1.1f5e7919d7ecp-8, 0x1.6d5p-56, 0x1.6b0ae555028ap-74},
  {10, 0x1.3f38a60f0648p-8, 0x1.2b4p-57, 0x1.3c937494045eap-70},
  {11, 0x1.5f0edcf18bdcp-8, -0x1.491p-56, -0x1.8c7950928041p-72},
  {12, 0x1.7ee11ebd82e8p-8, 0x1.3a8p-56, -0x1.e96e2fc5d8ff8p-70},
  {13, 0x1.9eaf6c6ea7c4p-8, -0x1.7fp-56, 0x1.d6aa5e18a19eep-70},
  {14, 0x1.be79c70058ecp-8, 0x1.1f4p-57, -0x1.64fefef02b628p-70},
  {15, 0x1.de402f6d9754p-8, -0x1.834p-58, 0x1.16e42c40b3202p-70},
  {16, 0x1.fe02a6b10678p-8, 0x1.1f8p-57, 0x1.bb481c8ee1418p-71},
  {17, 0x1.0ee096e2765p-7, -0x1.838p-57, -0x1.9b7bf9f8bc934p-71},
  {18, 0x1.1ebde2d1997ep-7, 0x1.7cp-57, 0x1.6e479344dfc4p-74},
  {19, 0x1.2e9937a2b2f2p-7, 0x1.75ap-56, -0x1.0df74bd4b1f58p-70},
  {20, 0x1.3e7295d25a7ep-7, -0x1.c2p-57, 0x1.acbdd778bf36p-74},
  {21, 0x1.4e49fddcf9aep-7, -0x1.54cp-56, 0x1.3212f4753101p-71},
  {22, 0x1.5e1f703ecbe6p-7, -0x1.f61p-56, 0x1.617686d277124p-70},
  {23, 0x1.6df2ed73de72p-7, 0x1.d44p-56, -0x1.5d6396fc72c74p-71},
  {24, 0x1.7dc475f810a8p-7, -0x1.246p-56, 0x1.289782c20df36p-70},
  {25, 0x1.8d940a4713ecp-7, -0x1.8ccp-58, 0x1.7b0357cf5e348p-71},
  {26, 0x1.9d61aadc6bd8p-7, 0x1.963p-56, 0x1.3da12b69ae59cp-71},
  {27, 0x1.ad2d58336e4ep-7, -0x1.7d6p-56, -0x1.525bb5226dc58p-71},
  {28, 0x1.bcf712c74384p-7, 0x1.782p-56, 0x1.7bd9770b665bp-70},
  {29, 0x1.ccbedb12e62ep-7, -0x1.4fp-60, -0x1.24e9f34718388p-72},
  {30, 0x1.dc84b1912382p-7, -0x1.6bp-56, 0x1.e3a5b304fbe58p-70},
  {31, 0x1.ec4896bc9b58p-7, -0x1.bf4p-56, 0x1.b080db0803b18p-72},
  {32, 0x1.fc0a8b0fc03ep-7, 0x1.e7cp-58, 0x1.eda74d37abd58p-71},
  {33, 0x1.05e547826bc9p-6, 0x1.0fp-59, -0x1.f5e1f7e49e22p-72},
  {34, 0x1.0dc4518afcc8p-6, 0x1.e86p-56, 0x1.ff3b0100c426p-71},
  {35, 0x1.15a263de88b9p-6, 0x1.53ap-57, 0x1.d76157c5a14b8p-72}
};

static double
algo73 (double x11, double *x21, double *x22, double *x23)
{
  int n2 = R0 (x11 * 0x1p11);
  /* We should have -35 <= n2 <= 35.  */
  const double *p = table2[35 + n2];
  *x21 = x11 - p[1];
  *x22 = -p[2];
  *x23 = -p[3];
  return p[0]; /* m2 */
}

static const double table3[53][4] = {
  {-26, -0x1.a015216e468p-12, -0x1.12fp-56, -0x1.f7ea523c98dp-77},
  {-25, -0x1.901389459d4p-12, 0x1.14ap-56, -0x1.49dac2a33d262p-70},
  {-24, -0x1.80120120144p-12, -0x1.85p-60, 0x1.4d0c2247f3678p-72},
  {-23, -0x1.701088fd8bcp-12, -0x1.a4p-62, -0x1.0fb599498e5cp-73},
  {-22, -0x1.600f20dde3cp-12, 0x1.d6fp-56, -0x1.78500a88546ap-73},
  {-21, -0x1.500dc8c0fbcp-12, -0x1.f8bp-56, -0x1.e397f6f734088p-70},
  {-20, -0x1.400c80a6b48p-12, 0x1.0b9p-56, 0x1.650f877c54dfp-73},
  {-19, -0x1.300b488eed4p-12, -0x1.424p-57, -0x1.1033a4d1ee5c2p-70},
  {-18, -0x1.200a2079868p-12, 0x1.764p-56, -0x1.24078c1d178d8p-70},
  {-17, -0x1.100908665fcp-12, -0x1.fap-59, -0x1.5bb8ff3ad374p-75},
  {-16, -0x1.00080055594p-12, -0x1.589p-56, 0x1.d332a0e20e2f2p-70},
  {-15, -0x1.e00e108ca6p-13, -0x1.729p-56, -0x1.5ff93a04bb768p-72},
  {-14, -0x1.c00c40725ap-13, -0x1.828p-59, 0x1.1561980a5e7dp-72},
  {-13, -0x1.a00a905b8ep-13, -0x1.3abp-56, -0x1.dd0d7af9e5e76p-70},
  {-12, -0x1.80090048028p-13, -0x1.03p-58, -0x1.36ff5a42842fcp-71},
  {-11, -0x1.6007903777p-13, -0x1.eeep-57, 0x1.dfa826ded53e8p-70},
  {-10, -0x1.40064029acp-13, 0x1.cccp-57, -0x1.bbf03ae68d184p-71},
  {-9, -0x1.2005101e61p-13, 0x1.979p-56, 0x1.dcd92b45cfcc8p-72},
  {-8, -0x1.0004001556p-13, 0x1.554p-56, -0x1.113bbce056052p-70},
  {-7, -0x1.c006201c96p-14, 0x1.498p-58, -0x1.3b842ad90132cp-71},
  {-6, -0x1.80048012p-14, -0x1.44p-56, -0x1.84d464f3db50cp-70},
  {-5, -0x1.4003200a6bp-14, 0x1.722p-57, 0x1.7225947f889a8p-71},
  {-4, -0x1.0002000555p-14, -0x1.955p-56, -0x1.888933357c5fcp-70},
  {-3, -0x1.800240048p-15, -0x1.44p-60, -0x1.84d098d69054p-75},
  {-2, -0x1.0001000156p-15, 0x1.515p-56, 0x1.53bbb9110c7ecp-70},
  {-1, -0x1.0000800054p-16, -0x1.559p-56, -0x1.5562222cccd6p-70},
  {0, 0x0p3, 0x0p3, 0x0p3},
  {1, 0x1.ffff0000a8p-17, 0x1.551p-56, 0x1.556222177780ap-70},
  {2, 0x1.fffe0002acp-16, -0x1.595p-56, -0x1.53bbbe6661d42p-70},
  {3, 0x1.7ffdc0048p-15, -0x1.44p-60, 0x1.84c900d6902p-75},
  {4, 0x1.fffc000aaap-15, 0x1.155p-56, 0x1.8887dde026fa8p-70},
  {5, 0x1.3ffce00a6bp-14, -0x1.f19p-56, -0x1.722fc0aa3404p-71},
  {6, 0x1.7ffb8012p-14, -0x1.44p-56, 0x1.84c534f3d9b6ap-70},
  {7, 0x1.bff9e01c95p-14, -0x1.02ep-56, 0x1.3b3792ae4b894p-71},
  {8, 0x1.fff8002aaap-14, -0x1.554p-56, 0x1.10e6678af0afcp-70},
  {9, 0x1.1ffaf01e5fp-13, 0x1.97fp-56, -0x1.df8d264674028p-72},
  {10, 0x1.3ff9c029a98p-13, -0x1.b98p-58, 0x1.b965303b23b18p-71},
  {11, 0x1.5ff87037738p-13, 0x1.7bcp-58, -0x1.e1e8d4f4f701ap-70},
  {12, 0x1.7ff70047fd8p-13, -0x1.fap-59, 0x1.2f675a3f500fcp-71},
  {13, 0x1.9ff5705b87p-13, 0x1.723p-56, 0x1.d6ea40e186a1p-70},
  {14, 0x1.bff3c072508p-13, 0x1.285p-56, -0x1.3badad75b1828p-72},
  {15, 0x1.dff1f08c9ap-13, -0x1.6dfp-56, 0x1.2609aede8acf8p-72},
  {16, 0x1.fff000aaa28p-13, 0x1.589p-56, -0x1.e887f64763848p-70},
  {17, 0x1.0ff6f866558p-12, 0x1.1ecp-56, -0x1.3d380357b53ap-74},
  {18, 0x1.1ff5e079798p-12, 0x1.81cp-56, 0x1.f18fb7e814328p-71},
  {19, 0x1.2ff4b88edd8p-12, -0x1.e76p-56, 0x1.a8c274faa879cp-71},
  {20, 0x1.3ff380a6a1p-12, -0x1.8b9p-56, -0x1.f80d1a90f8068p-72},
  {21, 0x1.4ff238c0e44p-12, -0x1.dfdp-56, 0x1.76894daa4d574p-70},
  {22, 0x1.5ff0e0ddc7p-12, 0x1.422p-57, -0x1.84861342e1fd8p-72},
  {23, 0x1.6fef78fd698p-12, 0x1.761p-56, -0x1.348c241fbd71p-71},
  {24, 0x1.7fee011febcp-12, 0x1.85p-60, -0x1.46430a2c0cdcep-70},
  {25, 0x1.8fec79456d8p-12, -0x1.48p-62, 0x1.369e45578a58p-74},
  {26, 0x1.9feae16e0ecp-12, 0x1.15ep-57, -0x1.84deb3bacd2c2p-70}
};

static double
algo74 (double x21, double *x31, double *x32, double *x33)
{
  int n3 = R0 (x21 * 0x1p16);
  /* We should have -26 <= n3 <= 26.  */
  const double *p = table3[26 + n3];
  *x31 = x21 - p[1];
  *x32 = -p[2];
  *x33 = -p[3];
  return p[0]; /* m3 */
}

/* page 1398, left column */
static inline double
P (double v)
{
  double a0 = -0x1p-2;
  double a1 = 0x1.99999999d2ac6p-3;
  double a2 = -0x1.55555555a890cp-3;
  return a0 + v * (a1 + a2 * v);
}

/* Algorithms 7.5 and 7.6 */
static double
algo75 (double x)
{
  double d, x3, Y, Y1, Y3, Y5, S1, Sprime1, Sprime2, v4, v42, M1, M2, M3, m0;
  int n;
  {
    /* We should have |x| < 745.2.  */
    double t;
    double x01, x02, x03;
    algo71 (x, &m0, &x01, &x02, &x03);
    double m1, x11, x12, x13;
    m1 = algo72 (x01, &x11, &x12, &x13);
    double m2, x21, x22, x23;
    m2 = algo73 (x11, &x21, &x22, &x23);
    double m3, x31, x32, x33;
    m3 = algo74 (x21, &x31, &x32, &x33);
    double x1 = x31 + (x12 + (x22 + x32)); /* defined in Corollary 7.1 */
    double x2 = x02;
    x3 = ((x03 + x13) + x23) + x33;
    /* v1 = m1/2^6, v2 = m2/2^11, v3 = m3/2^16 */
    double Mprime0 = (1 + m1 * 0x1p-6) * (1 + m2 * 0x1p-11) * (1 + m3 * 0x1p-16);
    double Mprime1 = Rm16 (Mprime0);
    double vprime4 = Rm34 (x1 * (1.0 / 3.0));
    v4 = 3.0 * vprime4;
    t = (Mprime0 - Mprime1) + Mprime0 * v4;
    double Mprime2 = Rm35 (t);
    double Mprime3 = t - Mprime2;
    if (m0 < -1022)
      n = m0 + 1022;
    else
      n = 0;
    M1 = ldexp (Mprime1, n);
    M2 = ldexp (Mprime2, n);
    M3 = ldexp (Mprime3, n);
    v42 = v4 * v4;
    fast_two_sum (&Sprime1, &Sprime2, (x1-v4) + x2 + 0.5 * v42, -v42 * vprime4);
    S1 = Rm66 (Sprime1);
    double Y4;
    fast_two_sum (&Y4, &Y5, M3 + M1 * S1, M2 * S1);
    double Y2;
    fast_two_sum (&Y2, &Y3, M1 + M2, Y4);
    if (m0 < -1022 || (m0 == -1022 && Y2 < 1))
    {
      d = 1;
      fast_two_sum (&Y, &Y1, 1, Y2);
    }
    else
    {
      d = 0;
      Y = Y2;
      Y1 = Y3;
      Y3 = 0;
    }
  }
  double D0 = 0;
  if (rnd == FE_TONEAREST)
  {
    if (Y1 >= 0)
      D0 = (Y < 1) ? 0x1p-54 : 0x1p-53;
    else
      D0 = (Y < 1) ? -0x1p-54 : -0x1p-53;
  }
  double D1 = (Y1 - D0) + Y3;
  if (fabs (D1) < ldexp (1.567, -65))
  {
    /* FIXME: we could compute Sprime1/Sprime2 only in this case */
    double S2 = (Sprime1 - S1) + (x3 + (Sprime2 - v42 * v42 * P (v4)));
    double M = M1 + M2 + M3;
    double S0 = S2 + (S1 + S2) * (S1 + S2) * (0.5 + (S1 + S2) / 6);
    D1 = D1 + M * S0;
    if (fabs (D1) < ldexp (1.523, -69))
      D1 = D1 + (M3 * S1 + Y5);
  }
  if (rnd == FE_TONEAREST)
  {
    if ((D0 > 0 && D1 > 0) || (D0 < 0 && D1 < 0))
      Y = Y + (D0 + D0);
  }
  else if (rnd == FE_UPWARD)
  {
    if (D1 > 0)
      Y = Y + ((Y >= 1) ? 0x1p-52 : 0x1p-53);
  }
  else /* FE_DOWNWARD */
  {
    if (D1 < 0)
      Y = Y - ((Y > 1) ? 0x1p-52 : 0x1p-53);
  }
  return ldexp (Y - d, m0 - n);
}

double
cr_exp (double x)
{
  if (isnan (x))
    return x + x; /* always return qNaN, even for sNaN input */

  /* Other special cases, cf page 1390, left column.  */

  /* Deal with overflow.  */
#define X_MAX 0x1.62e42fefa39efp+9
  if (x > X_MAX)
  {
    /* The +0x1p970 returns 0x1.fffffffffffffp1023 or +Inf according to the
       rounding mode.  */
    double x0 = 0x1.fffffffffffffp1023; /* largest representable number */
    /* If x > x0, then x is +Inf.  */
    return (x > x0) ? x : x0 + 0x1p970;
  }
  
  /* deal with underflow */
#define X_MIN -0x1.6232bdd7abcd2p+9
#define X_DMIN -0x1.74385446d71c3p+9
#define X_DMIN2 -0x1.74910d52d3051p+9
  if (x < X_MIN)
  {
    if (x < X_DMIN && rnd != FE_TONEAREST)
    {
      /* Warning: 0x1p-1074 * 0.5 is optimized by gcc to 0, whatever the
         rounding mode: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=34678.
         We add 'volatile' which seems to solve this issue. */
      volatile double xmin = 0x1p-1074;
      return xmin * 0.5;
    }
    if (x < X_DMIN2 && rnd == FE_TONEAREST)
      return 0;
  }

  if (fabs (x) <= 0x1p-36)
    return algo41 (x);

  if (-0x1p-28 <= x && x <= 0x1p-27)
    return algo51 (x);
  
  return algo75 (x);
}
