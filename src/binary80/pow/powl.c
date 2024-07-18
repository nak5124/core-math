/* Correctly rounded powl function for binary80 values.

Copyright (c) 2024 Sélène Corbineau and Paul Zimmermann

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
#include <fenv.h>
#include <stdbool.h>
#include <stdio.h>

#ifdef __x86_64__
#include <x86intrin.h>
#endif

#ifdef POWL_DEBUG
#include <stdio.h>
#define POWL_DPRINTF(...) printf(__VA_ARGS__)
#define SAGE_RR "R(\"%a\",16)"
#define SAGE_RE "R(\"%La\",16)"
#define SAGE_DD "(R(\"%a\",16)+R(\"%a\",16))"
#else
#define POWL_DPRINTF(...)
#endif

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {long double f; struct {uint64_t m; uint16_t e;};} b80u80_t;
typedef union {
	double f;
	struct __attribute__((packed)) {uint64_t m:52;uint32_t e:11;uint32_t s:1;};
	uint64_t u;
} b64u64_t;

static inline int get_rounding_mode (void)
{
#ifdef __x86_64__
  const unsigned flagp = _mm_getcsr ();
  return (flagp&(3<<13))>>3;
#else
  return fegetround ();
#endif
}

/* Split a number of exponent 0 (1 <= |x| < 2)
   into a high part fitting in 33 bits and a low part fitting in 31 bits. */
static inline
void split(double* rh, double* rl, long double x) {
	static long double C = 0x1.8p+31L; // ulp(C)=2^-32
	long double y = (x + C) - C;
	/* Given the relative sizes of C and x, x + C has the same binade as C.
           Therefore, the difference is exact. Furthermore,
           ulp(x + C) = ulp(C) = 2^-32.
           The rounding error in x + C is therefore less than 2^-32.
	   Thus, |x - y| < 2^-32. Note that since 2^31 <= x + C < 2^32 and the
	   difference is exact, y is a multiple of ulp(x + C) = 2^-32.
           Since |x| < 2, and the roundings are monotonous, x + C is bounded
           by the values obtained with |x| = 2, namely 0x1.7ffffffcp+31 and
           0x1.80000004p+31, and likely for y, namely -2 and 2.
           Since y is a multiple of 2^-32, this ensures y = k*2^-32
           with |k| <= 2^-33, thus y fits in 33 bits.
           (If |y| = 2, it trivially fits.) */
	*rh = y; // This conversion is exact by the argument above.
	*rl = x - y;
	/* 
           |x - y| < 2^-32. Note that x and y are both multiples of
           ulp_64(1) = 2^-63; therefore x - y too. This implies that
           x - y = l*2^-63 with |l| < 2^31, thus rl fits in 31 bits,
           and the difference is exact. */
}


static inline
void add22(double* zh, double* zl, double xh, double xl, double yh, double yl) {
	double r,s;
	r = xh+yh;
	s = ((xh-r)+yh)+yl+xl;
	*zh = r+s;
	*zl = (r - (*zh)) + s;
}

static inline
void fast_two_sum(double* rh, double* rl, double a, double b) {
	*rh = a + b;
	double e = *rh - a;
	*rl = b - e;
}

static inline
void two_sum(double* rh, double* rl, double a, double b) {
	*rh = a + b;
	double ap = *rh - b;
	double bp = *rh - ap;
	*rl = (a - ap) + (b - bp);
}

/* Computes an approximation of a + bh + bl assuming a = 0 or |a| >= |bh|*/
static inline
void high_sum(double* rh, double* rl, double a, double bh, double bl) {
	double e;
	fast_two_sum(rh, &e, a, bh);
	*rl = bl + e;
}

/* Computes rh + rl = a * b exactly */
static inline
void a_mul(double* rh, double* rl, double a, double b) {
	*rh = a*b;
	*rl = __builtin_fma(a,b,-*rh);
}

/* Computes an approximation of (ah+al)(bh+bl)-al*bl */
static inline
void d_mul(double* rh, double* rl, double ah, double al,
                                   double bh, double bl)
{ double p;
	a_mul(rh, &p, ah, bh);
	double q = __builtin_fma(al, bh, p);
	*rl = __builtin_fma(ah, bl, q);
}

#include "powl_tables.h"

/* Let x = xh + xl. Assume |x| <= 2^-12
Then polyeval(&rh, &rl, xh, xl) returns in 
rh + rl an estimate of log2(1 + x) with relative error at most 2^-98.429.
*/
static inline
void polyeval(double* rh, double* rl, double xh, double xl) {
	/* We approximate log2(1 + x) by x/ln(2) * (c0 + c1*x + ... + c7*x^7)
	  This polynomial has intrinsic relative error 2^-105.879
	  (the biggest error we make is in x^4/5 ~> 2^(-52)*2^(-42) = 2^-94
	*/
	double ln2invh = 0x1.71547652b82fep+0, ln2invl = 0x1.777d0ffda0d24p-56;
	double scaleh, scalel; d_mul(&scaleh, &scalel, ln2invh, ln2invl, xh, xl);
	/* Expanding the d_mul call we get that :
	   - |p| <= 2^-52|ah*bh|
	   - |al*bh+p| <= 2^-55.976|ah*bh| + 2^-104|ah*bh| so that q is such that
	     |q| <= 2^-55.975|ah*bh| and q's rounding error is
	     at most 2^-107.975|ah*bh|.
	   - |ah*bl+q| <= 2^-52|ah*bh| + 2^-55.975|ah*bh| so that
	     |scalel| <= 2^-51.911|ah*bh| and the final rounding error is at most
	     2^-103.911|ah*bh|.
	   - the error neglecting |al*bl| is at most 2^-52*2^-55.976|ah*bh|
	   The total relative error in terms of |ah*bh| is thus at most
	     2^-107.976 + 2^-103.91 + 2^-107.975 <= 2^-103.747.
	   Expressing the error relative to ln2inv*x, we get a relative error at most
	     2^-103.746.
	*/

	double ord01h, ord01l;
	ord01h = -xh/2; ord01l = -xl/2; /*c1 = 1/2*/ // Exact.
	high_sum(&ord01h,&ord01l, 1, ord01h, ord01l);
	/* Expanding the high_sum call we get that :
	   - the fast_two_sum incurs an error of 2^-105(1 + 2^-12.999) <= 2^-104.999
	     and its' low value is at most 2^-52.
	   - the final sum has value at most 2^-52+2^-12.999*2^-52 <= 2^-51.999; this
	     implies a rounding error of at most 2^-104.
	  Therefore, at output we have |ord01h| <= 1 + 2^-12.998,
	  |ord01l| <= 2^-51.999 and the total error on 1 - x/2 is at most 2^-103.414.
	*/


	double ord23h, ord23l;
	ord23h = -xh/4; /*c3 = -1/4*/
	ord23l = __builtin_fma(-xl, 1./4, 0x1.55555a5b705aap-56); /*c2l*/
	high_sum(&ord23h, &ord23l, 0x1.5555555555555p-2, ord23h, ord23l); /*c2h*/
	/* We compute that
	    |-xl/4 + c2l| <= 2^-11.999*2^(-54) + 2^-55.584 <= 2^-55.582.
	   Therefore after the fma |ord23l| <= 2^-55.581 and the associated rounding
	   error is at most 2^-56-52 = 2^-108.

	   Expanding the high_sum call we get that :
	   - the fast_two_sum incurs an error of 2^-105(2^-1.584 + 2^-13.999)
	     <= 2^-106.583 and its' low value is at most 2^-54.
	   - the final sum has value at most 2^-54 + 2^-55.581 <= 2^-53.583; this
	     implies a rounding error of at most 2^-106.
	   At output we thus have |ord23h| <= 2^-1.583, |ord23l| <= 2^-53.582
	   and the total error on c2 + c3x is at most
	      2^-106+2^-108+2^-106.583 <= 2^-105.060.
	*/

	double xsqh, xsql; d_mul(&xsqh, &xsql, xh, xl, xh, xl);
	// FIXME directly analyze the absolute errors ?
	/* Expanding the d_mul call we get that :
	   - |p| <= 2^-52|xh^2|
	   - |xh*xl+p| <= 2^-51|xh^2| so that |q| <= 2^-50.999|xh^2| and q's rounding
	     error is at most 2^-102.999|xh^2|.
	   - |xl*xh+q| <= 2^-50.414|xh|^2 so that |xsql| <= 2^-50.413|xh^2| and the
	     rounding error is at most 2^-102.413|xh^2|
	   - the error made by neglecting xl^2 is at most 2^-104|xh|.
	   The total relative error in terms of |xh^2| is thus at most
	     2^-102.413 + 2^-102.999 + 2^-104 <= 2^-101.413
	   Expressing the error relative to x^2 we get a relative error at most
	     2^-101.412.
	   This translates to an absolute error less than 2^-125.41 
	   Also, at output we have |xsqh| <= 2^-23.997 and
	   |xsql| <= 2^-50.413|xsqh| <= 2^-74.410 .
	*/

	d_mul(&ord23h, &ord23l, ord23h, ord23l, xsqh, xsql);
	/* Expanding the d_mul call we get that:
	   - |p| <= 2^-52 * (2^-23.997*2^-1.583) <= 2^-77.579.
	   - |al*bh + p| <= 2^-53.582*2^-23.997 + 2^-77.579 <= 2^-76.578, so that 
	     |q| <= 2^-76.577 and the associated rounding error is at most 2^-129.
	   - |ah*bl + q| <= 2^-1.583*2^-74.410 + 2^-76.578 <= 2^-75.256. Therefore
	     at output |ord23l| <= 2^-75.255 and the associated rounding error is at
	     most 2^-128.
	   - the error made by neglecting xsql*ord23l is at most 2^-74.410*2^-53.582
	     <= 2^-127.991
	   Propagating the errors on x^2 and c2 + c3x gives an intrinsic error of
	   at most :
	     2^-125.41 * (c2 + 2^-11.999*|c3|) + 2^-105.060*2^(-11.999*2)
	     + 2^-125.41 * 2^(-11.999*2) <= 2^-126.685
	   The total absolute error computing x^2(c2 + c3x) is thus bounded by
	     2^-128+2^-129+2^-127.991+2^-126.685 <= 2^-125.679.
	*/

	double x4 = xsqh*xsqh;
	/* Neglecting 2*xsqh*xsql + xsql^2 creates an error of at most
	   2*2^-23.997*2^-74.410 + 2^(-74.410*2) <= 2^-97.406.
	   Also, |xsqh*xsqh| <= 2^(-23.997*2) <= 2^-47.994. Therefore, the rounding
	   error of the product is at most 2^-100.

	   We obtain that |x4| <= 2^-47.993
	   and that |x4 - x^4| <= 2^-97.406 + 2^-100 + |(xsqh+xsql)^2 - x^4|
	                       <= 2^-97.184 + |(xsqh+xsql)+x^2||(xsqh+xsql) - x^2|
	                       <= 2^-97.184 + (2*(2^-11.999)^2 + 2^-125.41)*2^-125.41
	                       <= 2^-97.183
	*/

	double acc = __builtin_fma(xh, -0x1.555555555554dp-3/*c5*/,
	                                0x1.999999999998ap-3/*c4*/);
	double bcc = __builtin_fma(xh, -0x1.0000014f8ec21p-3/*c7*/,
	                                0x1.24924ad7557bep-3/*c6*/);
	/* We compute that |xh*c5 + c4| <= 2^-2.321. This implies that
	   |acc| <= 2^-2.321 and that the rounding error is at most 2^-55.
	   Neglecting xl*c5 incurs an error of at most 2^-64*2^-2.584 <= 2^-66.584
	   The total error computing c4 + c5x is thus at most
	     2^-55+2^-66.584 <= 2^-54.999.

	   In the same way, we compute that |xh*c7 + c6| <= 2^-2.807, which implies
	   |bcc| <= 2^-2.806 and that the rounding error is at most 2^-55.
	   Neglecting xl*c7 incurs an error of at most 2^-64*2^-2.999 <= 2^-66.999.
	   The total error computing c6 + c7x is thus at most
	     2^-55+2^-66.999 <= 2^-54.999.
	*/

	acc = __builtin_fma(xsqh, bcc, acc);
	/* We compute that |xsqh*bcc + acc| <= 2^-23.997*2^-2.806+2^-2.321<= 2^-2.320.
	   This ensures that at output, |acc| <= 2^-2.319. Also, the rounding error is
	   at most 2^-55.
	   Since |xsqh+xsql - x^2| <= 2^-125.41, we have
	   |xsqh - x^2| <= 2^-125.41+2^-74.410 <= 2^-74.409. Propagating the errors
	   on acc and bcc we get an intrinsic error of 
	     2^-74.409*2^-2.806 + 2^-54.999*2^(-11.999*2) + 2^-54.999*2^(-11.999*2)
	     + 2^-54.999 <= 2^-54.998.
	   The total error computing c4 + c5x + x^2(c6 + c7x) is thus at most
	    2^-54.998 + 2^-55 <= 2^-53.998 
	*/

	ord01l = __builtin_fma(x4, acc, ord01l);
	/* Propagating the errors on x4 and acc yields an intrinsic error of
	     2^-54.999*2^-47.996 + 2^-97.183*2^-2.319 + 2^-97.183*2^-54.999
	     <= 2^-99.379.
	   We compute that |x4*acc+ord01l|<=2^-47.993*2^-2.319+2^-51.999 <= 2^-49.921.
	   This implies that at output |ord01l| <= 2^-49.920 and that the rounding
	   error is 2^-102 at most.
	   The error in this step is thus 2^-102+2^-99.379 = 2^-99.161 
	*/

	high_sum(&ord23h, &ord23l, ord01h, ord23h, ord23l);
	ord23l += ord01l; // Rewrite ? 
	//add22(&ord01h, &ord01l, ord01h, ord01l, ord23h, ord23l);
	/* Propagating the errors on subterms gives a total intrinsic error of
	   2^-99.161 + 2^-105.060 + 2^-103.414 <= 2^-99.064.

	   Expanding high_sum, we see that the fast_two_sum incurs an error at most
	   2^-105(1 + 2^-12.998 + 2^-1.583) <= 2^-104.584. We also see that
	   |ord23h| <= 2^0.416 at output. Furthermore the fast_two_sum's low value
	   is at most 2^-52. The final sum has absolute value at most
	   2^-52+2^-53.582 <= 2^-51.584. This ensures that |ord23l| <= 2^-51.583
	   after the high_sum and that the associated rounding error is at most
	   2^-104.

	   Given that |ord01l| <= 2^-49.920 we get that the final sum is less than
	   2^-49.920+2^-51.583 <= 2^-49.524. This ensures that at output
	   |ord23l| <= 2^-49.523 and that the sum's rounding error is at most 2^-102.

	   The total accumulated error is therefore at most
	   2^-102 + 2^-104 + 2^-104.584 + 2^-99.064 <= 2^-98.818.
	   We have tried to compute log2(1 + xr)/(ln2inv*xr), which 
	   is positive and decreasing by concavity. We can check that it's value
	   at xr = 2^-11.999 is at least .9998. Therefore, the absolute error above
	   translates to a relative error of 2^-98.818/.9998 <= 2^-98.817. Added
	   to the polynomial's intrinsic relative error of 2^-95, we ensure
	   that we have computed log2(1+xr)/(ln2inv*xr) with relative error at most
	   2^-105.879+2^-98.817 <= 2^-98.806
	*/

	d_mul(rh, rl, scaleh, scalel, ord23h, ord23l);
	/* Propagating each term's relative errors we get a total intrinsic relative
	   error of (1 + 2^-103.746)*(1 + 2^-98.806) - 1 <= 2^-98.759
	   Expanding the d_mul call we get :
	     - |p| <= 2^-52*|scaleh|*1.0002
	     - |scalel*ord23h + p| <= (2^-51.910+2^-52)*1.0002*|scaleh|. This ensures
	       that |q| <= 2^-50.954|scaleh| and that the rounding error is bounded by
	       2^-102.954|scaleh|.
	     - |scaleh*ord23l + q| <= (2^-49.523 + 2^-50.954)|scaleh|. This ensures
	       that |rl| <= 2^-49.067|scaleh| and that the final rounding error is
	       at most 2^-101.067|scaleh|.
	  The total rounding error is thus at most
	    (2^-101.067+2^-102.954)|scaleh| <= 2^-100.721|scaleh|.
	  Writing
	    |scaleh|/|log2(1+x)| = (|scaleh|/|scale*ord23|)*|scale*ord23|/|log2(1+x)|
		we see that the second factor is off from 1 by at most 2^-98.756 and that
	  the first factor is off from 1 by at most 0.0002. This ensures that the
	  rounding error can be expressed as 2^-100.720|log2(1+x)|.

	  The total relative error of this routine is therefore at most
	   2^-98.759 + 2^-100.721 <= 2^-98.429.
	  We have the postcondition |rl| <= 2^-49.066|rh|.
	*/
}

/* Computes an approximation of ylog2|x| under the following conditions :

- 2^-80 <= |y| < 2^78
- x is normal, i.e., |x| >= 2^-16382
*/
static inline
void compute_log2pow(double* rh, double* rl, long double x, long double y) {
	b80u80_t cvt_x = {.f = x};
	int extra_int = (cvt_x.e&0x7fff) - 16383;
        // -16382 <= extra_int <= 16383
	cvt_x.e = 16383; // New wanted exponent
	x = cvt_x.f;
        // original x = 2^extra_int * x

	double xh, xl; // a (resp b) bits
	split(&xh, &xl, x);
        // x = xh + xl with xh on 33 bits and xl on 31 bits

	POWL_DPRINTF("sx = " SAGE_RE "\nei = %d\n", x, extra_int);
	// Uses the high 7 bits of x's mantissa.
	lut_t l = coarse[cvt_x.m>>56 & 0x7f];
	POWL_DPRINTF("key=0x%lx\n", (cvt_x.m>>56 & 0x7f));

	/* We always have |x*r1 - 1| <= 0x1p-7. The term l.z is chosen such that
	   l.z+mlogr1h + mlogr1 approximates -log2(r1) with
     relative error < 2^-107.22. Note that |mlogr1h+mlogr1l| < .5
	   The tables are constructed in such a way that r fits in 9 bits.
	*/
	double r1      = l.r;
	double mlogr1h = l.mlogrh;
	double mlogr1l = l.mlogrl;
	extra_int     += l.z;

	POWL_DPRINTF("r1 = " SAGE_RR "\n", r1);

	// Eliminated if POWL_DEBUG is not defined	
	if(l.z) {POWL_DPRINTF("sx = sx/2\nei+=1\n");}
	xh *= r1; xl *= r1;
        /* The above multiplications are exact.
           now xh fits in 42 bits at most, xl in 40.
           More precisely the initial xh was a multiple of 2^-32.
	         Since r1 is a multiple of 2^-9 then the new value of xh is a multiple
           of 2^-41.
        */

	POWL_DPRINTF("get_hex(R(abs("SAGE_RR" - 1)))\n", xh);
	/* Note that now |xh - 1| <= 1p-7
	   Therefore, xh's mantissa (seen as a 53-bit integer) is either
	   1.00000 00p or 1.11111 1q.
	   We skip the upper 6 bits of the mantissa and use the next 7 bits
           to index another lookup table. A quarter of the table is wasted!

           We're looking at 1 + 2^-12*k, 1 + 2^-12*(k+1) for 0 <= k < 32.
	   Else we're looking at 1 - 2^-6 + 2^-13*k', 1 - 2^-6 + 2^-13*(k'+1)
           for 64 <= k' < 128.
	*/

	b64u64_t cvt_xh = {.f = xh};
	lut_t l2 = fine[(cvt_xh.u>>40) & 0x7f]; // k' = (cvt_xh.u>>40) & 0x7f
	// bit 52 goes to 6+5 = 11. Bits 11 - 8
	POWL_DPRINTF("key2 = 0x%lx\n", (cvt_xh.u>>40 & 0x7f));
	double r2 = l2.r;
	double mlogr2h = l2.mlogrh;
	double mlogr2l = l2.mlogrl;
	POWL_DPRINTF("r2 = " SAGE_RR "\n", r2);
	/* The fine table is built in such a way that :
	   i)  |r2*xh - 1| <= 2^-12
	   ii) r2 fits in 13 bits
	   iii) mlogr2h + mlogr2l approximates -log2(r2) with
	        relative error at most 2^-107.27
	   iv) |mlogr2h+mlogr2l| < 2^-6.47
           Since r2 is a multiple of 2^-13, xh*r2 is a multiple of 2^-54,
           and since |r2*xh - 1| <= 2^-12, then r2*xh - 1 is representable
           exactly on 42 bits.
	*/

	POWL_DPRINTF("get_hex(R(-log2(r1)-" SAGE_DD"))\n", mlogr1h, mlogr1l);
	POWL_DPRINTF("get_hex(R(-log2(r2)-" SAGE_DD"))\n", mlogr2h, mlogr2l);
	
	double mlogrh, mlogrl;
	high_sum(&mlogrh, &mlogrl, extra_int, mlogr1h, mlogr1l);
	/* Since |mlogr1h + mlogr1l| < .5, we indeed have extra_int = 0 or
	   |extra_int| > |mlogr1h|. If extra_int=0 everything is exact,
           and we get mlogrh=mlogr1h, mlogrl=mlogr1l.
	   Assume |extra_int| >= 1. Expanding the high_sum call, this implies
	   that the fast_two_sum introduces an error <= 2^-105|mlogrh| and
           that the low part of its result is at most 2^-52|mlogrh|. Notice
           that |mlogrh| >= .5 > |mlogr1h| so that
	   2^-52|mlogrh| >= 2^-52|mlogr1h| >= |mlogr1l|.
	   This implies that the "rl" sum of high_sum (i.e., mlogrl)
           is at most 2^-51|mlogrh| and that its rounding error is
           at most 2^-103|mlogrh|.

	   The total rounding error is at most
	     (2^-103+2^-105)|mlogrh| <= 2^-102.678 |mlogrh|:
             
             |mlogrh + mlogrl - (extra_int + mlogr1h + mlogr1l)|
             <= 2^-102.678 |mlogrh|.

           If one performs an exhaustive search on all possible values
           of extra_int (-16382 to 16383), on all rounding modes, and
           on all values of mlogr1h/mlogr1l, we obtain that the maximal
           relative error is bounded by 2^-105.003 |mlogrh|
           (see function analyse_first_high_sum() in powl.sage).
	*/

        double mlogr12h, mlogr12l;
	high_sum(&mlogr12h, &mlogr12l, mlogrh, mlogr2h, mlogr2l);
	mlogr12l += mlogrl;
	/* Let us prove that unless it is zero, |mlogrh| is in the same binade
           of |mlogr2h| or in a larger binade (so that the fast_two_sum
           condition is fulfilled).
	   If extra_int != 0, this is obvious because |mlogrh+mlogrl| >= .5,
           and |mlogr2h| < 2^-6.47.
	   Assume extra_int = 0. Then mlogrh = mlogr1h and looking at the
           tables we see that mlogr1h = 0 or
           |mlogr1h| >= 0x1.6fe50b6ef0851p-7 >= |mlogr2h| which allows us to
           conclude.

	   Expanding high_sum(), let t the low part of the fast_two_sum() call.
           As above the fast_two_sum yields an error <= 2^-105 |mlogr12h|
           and |t| <= 2^-52 |mlogr12h|.
	   In the last sum of the fast_two_sum(), notice that
           |mlogr2l| <= 2^-53 |mlogr2h| <= 2^-53 * 124.6 |mlogr12h|
           < 2^-46.03 |mlogr12h|.
	   The factor 124.6 is because several bits might cancel in the
           addition mlogrh + mlogr2h.
           The largest cancellation is obtained for extra_int=0, i1=126
           (the index of r1), and i2=31, where we get
           mlogrh = -0x1.6fe50b6ef0851p-7, mlogr2h = 0x1.6cf6ddd2611d4p-7,
           thus mlogr12h = -0x1.7716ce47b3e8p-14, |mlogr2h/mlogr12h| ~ 124.545.

	   Therefore, |mlogr2l + t| <= (2^-46.03 + 2^-52) |mlogr12h|
                                    <= 2^-46.007 |mlogr12h|.
           This implies that |mlogr12l| <= 2^-46.007 |mlogr12h| and that the
           rounding error of the high_sum's final sum is
           at most 2^-98.007 |mlogr12h|.

	   Note that the following sum mlogr12l + mlogrl is at most
	     |mlogr12l| + 2^-51 |mlogrh|
             <= (2^-46.007 + 125.545 * 2^-51) |mlogr12h|
           (where |mlogrl| <= 2^-51|mlogrh| was proven in the analysis of
           the first high_sum, and |mlogrh/mlogr12h| < 125.6 from the above
           cancellation example), we thus get:
           |mlogr12l| <= 2^-43.701 |mlogr12h| and a rounding error
	   of at most 2^-95.701 |mlogr12h|.

	   These steps of computation created an error at most
	     (2^-105 + 2^-98.007 + 2^-95.701) |mlogr12h|
             < 2^-95.433 |mlogr12h|
	   Propagating the previous errors gives another error term at most
	     125.6 (2^-107.22 + 2^-107.27 + 2^-102.678) |mlogr12h|
             < 2^-95.599 |mlogr12h|
	   (the terms 2^-107.22 and 2^-107.27 come from finite precision of
           the table coarse[] and fine[] respectively, and the factor 125.6
           from the cancellation from |mlogrh/mlogr12h|).
           We thus get a total relative error of at most:
           (2^-95.433 + 2^-95.599) < 2^-94.513 |mlogr12h|.
	*/
	POWL_DPRINTF("get_hex(R(-log2(r1) - log2(r2)+ei- "SAGE_DD"))\n",
		mlogr12h, mlogr12l);
	fast_two_sum(&mlogr12h, &mlogr12l, mlogr12h, mlogr12l);
	/* This renormalization incurs a relative error at most 2^-105. The total
	   relative error becomes at most 2^-96.223. This ensures
	   |mlogr12l| <= 2^-52|mlogr12h|.
	*/

	// |xh| <= 1p-12
	xh = __builtin_fma(r2, xh, -1); xl *= r2;
	// exact since xl fits in 40 bits and at least 11 bits are cancelled
	// when computing the fma for xh.

	two_sum(&xh, &xl, xh, xl); // We probably cannot use Fast2Sum
	/* At input, we have |xh| <= 1p-12 and |xl| < 2*2^-32 = 2^-31. Therefore at
	   output we have |xh| <= 2^-11.999, |xl| <= ulp(xh) <= 2^-64 and 
	   xr = xh + xl is such that |xr| <= 2^-11.999.

	   Note that without the two_sum, we have no guarantees on the relative
	   sizes of xh and xl but both are small (approx. xh <= 2^-12, xl <= 2^-33).
	   this may be enough ? (probably not)
	*/
	POWL_DPRINTF("get_hex(R(r1*r2*sx - 1 - "SAGE_DD"))\n", xh, xl);
	POWL_DPRINTF("s = r1*r2*sx - 1\n");
	POWL_DPRINTF("get_hex(s)\n");

	polyeval(rh, rl, xh, xl);
	/* By polyeval's error analysis, rh + rl gets an estimate of log2(1+x)
	   with relative error at most 2^-98.429. Furthermore |rl|<=2^-49.066|rh|.
	*/

	high_sum(rh, rl, mlogr12h,*rh,*rl);
	*rl += mlogr12l;
	/* Let us call rh', rl' the results of the computation, rh and rl the inputs.
	   Note that if mlogr12h != 0, then |mlogr12h| >= 0x1.5p-12 (manual check).
	   Given the argument above, this implies |mlogr12h| has at least the
	   binade of |rh|.

	   Expanding the high_sum call and calling t the fast_two_sum result's low
	   part, the previous argument ensures that the fast_two_sum creates an error
	   at most 2^-105|rh'| and that |t| <= 2^-52|rh'|. In the last sum, notice
	   that |rl + t| <= 2^-52|rh'| + 2^-49.066|rh|. Now, since
	   2^-11.469/|0x1.5p-12 - 2^-11.469| <= 2^3.447, we have |rh| <= 2^3.447|rh'|
	   and thus |rl + t| <= (2^-52+2^-49.066*2^3.447)|rh'|. This shows that after
	   the high sum, |rl'| <= 2^-45.601|rh'| and that the associated rounding
	   error is at most 2^-97.601|rh'|.

	   We compute that 0x1.5p-12/|0x1.5p-12 - 2^-11.469| <= 2^3.309
	   The sum on the last line is thus at most
	     2^-45.601|rh'| + 2^3.309*2^-52|rh'| <= 2^-45.440|rh'|
	   which implies that in the end |rl'| <= 2^-45.439|rh'| and that the
	   associated rounding error is at most 2^-97.539|rh'|.

	   These steps of computation created an error at most
	     (2^-97.601 + 2^-105 + 2^-97.539)|rh'|.
	   Propagating the previous errors gives another error term at most
	     (2^3.309*2^-96.223 + 2^3.447*2^-98.429)|rh'|.
	   The total relative error computing log2(1 + x) is therefore at most
	   2^-92.515. Also, |rl'| <= 2^-45.439|rh'|.
	*/

	double yh = y; double yl = y - (long double)(yh);
	POWL_DPRINTF("get_hex(R(log2(x)) - "SAGE_DD")\n", *rh, *rl);
	d_mul(rh, rl, yh, yl, *rh, *rl);
	/* Let us call again rh', rl' the output values for rh and rl, and rh and rl
	   the input values. The relative error on log2(1 + x) propagates, creating
	   an intrinsic relative error of 2^-92.515.
	   Expanding the d_mul call, we see that |p| <= 2^-52|yh*rh|; then
	     - |yl*rh + p| <= 2^-52|yh*rh| + 2^-52|yh*rh| <= 2^-51|yh*rh|.
	       This ensures |q| <= 2^-51|yh*rh|. Also, the associated
	       rounding error is at most 2^-103|yh*rh|.
	     - |yh*rl + q| <= (2^-45.439 + 2^-51)|yh*rh|. This ensures that 
	       |rl'| <= 2^-45.408|rh'| and that the associated rounding error is
	       at most 2^-97.408|rh'|.
	     - the error produced by neglecting |yl*rl| is at most
	       2^-52*2^-45.439|yh*rh|.
	   The total error up to this step is therefore at most
	     (2^-52*2^-45.439 + 2^-97.408 + 2^-103 + 2^-92.515)|yh*rh|
		 or at most 2^-92.421|rh'|. We also have |rl'| <= 2^-45.408|rh|.
	*/
}


/* computes 2^(xh + xl), assuming |xl| <= ulp(xh) and |xh| < 2^31 */
static inline
int exp2d(double* resh, double* resl, double xh, double xl) {
	b64u64_t cvt = {.f = xh};
	bool do_red	= cvt.e >= -20 + 0x3ff;

	static const double C = 0x1.8p+32; // ulp is 2^-20
	b64u64_t y = {.f = xh + C};
	uint64_t fracpart = y.u;
	int16_t extra_exponent = y.u>>20;

	if(__builtin_expect(do_red, 1)) {
		double rem = xh - (y.f - C);
		fast_two_sum(&xh,&xl,rem,xl);
	}
	/* Let xl_old/xh_old be the old values to xh and xl, and xh/xl the values
	   afterwards. We have |rem| < 2^-20 and |xl| <= 2^(14.1-45.408) <= 2^-31.307.
	   Therefore |rem+xl| < 2^-20 + 2^(14.1-45.408) which ensures
	   |xh| <= 2^-19.9994.
	   If rem = 0 or rem's exponent is at least that of xl_old, then the rounding
	   error is at most 2^-105|xh| <= 2^-124.99942 and |xl| < ulp(xh) <= 2^-72.
	   If this is not the case, since |xl_old| <= 2^-31.307 we must have
	   |rem| < 2^-32 which implies |xh| < 2^-30.612. (Theorem 2 from reference [])
	   We also have |xl| <= 2^-72 in that case too.
	*/

	int i0 = fracpart & 0x1f;
	int i1 = (fracpart >> 5) & 0x1f;
	int i2 = (fracpart >> 10) & 0x1f;
	int i3 = (fracpart >> 15) & 0x1f;

	double frcp_acc0_l, frcp_acc0_h, frcp_acc2_h, frcp_acc2_l;
	double xs_pow2_h, xs_pow2_l;

	d_mul(&frcp_acc0_h, &frcp_acc0_l,
		t0[i0][0], t0[i0][1],   // 2^(i0/2^20)
		t1[i1][0], t1[i1][1]);  // 2^(i1/2^15)
	d_mul(&frcp_acc2_h, &frcp_acc2_l,
		t2[i2][0], t2[i2][1],   // 2^(i2/2^10)
		t3[i3][0], t3[i3][1]);  // 2^(i3/2^5)
	d_mul(&xs_pow2_h, &xs_pow2_l, frcp_acc0_h, frcp_acc0_l,
		frcp_acc2_h, frcp_acc2_l);
	/* This step introduces relative error |rho2| at most 2^-99.1, see sc_expl*/

	/* Evaluating the Taylor polynomial for 2^xr where xr = xh + xl.
	   If do_red is true, then |xh| <= 2^-19.9994 and |xl| <= 2^-72 so that
	   |xr| <= 2^-19.999.
	   If do_red is false, we have |xr| < 2^-20 + 2^(-20-45.408) <= 2^-19.999,
	   and |xl| <= 2^(-20-45.408) = 2^-65.408.

	   Over the interval [-2^-19.999, 2^-19.999] the polynomial used has 2^-89.218
	   absolute error.
	*/

	double xsq = xh * xh;
	/* Neglecting 2*xl*xh + xl^2 brings an error of at most
           |2*xl*xh + xl^2| <= 2 * 2^-65.408 * 2^-19.999 + 2^(-65.408*2)
                            <= 2^-84.406.
	   Since |xh| <= 2^-19.999, we have |xh*xh| <= 2^-39.998. The rounding
           error on xsq is therefore at most ulp(2^-39.998) = 2^-92, and
	     |xsq| <= 2^-39.998 + 2^-92 <= 2^-39.997.
	   We have thus |xsq - xr^2| <= 2^-84.406 + 2^-92 <= 2^-84.398.
	*/
	double orders23 = xsq * __builtin_fma(xh,0x1.c6b08d704a1cdp-5,
		0x1.ebfbdff82c696p-3);

	/* We note A = 0x1.c6b08d704a1cdp-5 and B = 0x1.ebfbdff82c696p-3.
	   Analyzing the fma call:
	   Neglecting xl * A imparts an error bounded by
	     |A * xl| <= A*2^-65.408 <= 2^-69.579 
	   Since |xh| <= 2^-19.999, |A*xh+B| <= 2^-19.999*A + B < 1/4 - 0.008.
	   This implies that the fma's result is strictly less than 1/4.
	   The rounding error of the fma is therefore at most ulp(1/8) = 2^-55.
	   The total error evaluating Axr+B is thus at most
           2^-55 + 2^-69.579 <= 2^-54.999.

	   Given the errors on each factor, the product xsq * fma() carries an
           intrinsic error
	     |xr|^2*2^-54.999 + 2^-84.398*|Axr+B| + 2^-84.398*2^-54.999
	   Since |xr|<= 2^-19.999 we check that |Axr+B| <= 2^-2.049 and thus
           the error is at most 2^-86.443.

	   Since |xsq| <= 2^-39.997 and |A*xh+B| < 1/4 - 0.008, the product
           is bounded by 2^-39.997 * (1/4 - 0.008) < 2^-42.04. The rounding
           error on orders23 is thus at most ulp(2^-42.04) = 2^-95.
           Therefore |orders23| <= 2^-42.04 + 2^-95 < 2^-42.03.
	   The total error on orders23 is at most 2^-86.443+2^-95 <= 2^-86.439.
	*/

	double order1h, order1l;
	static const double coeff1h = 0x1.62e42fefa39efp-1;
	static const double coeff1l = 0x1.abc9e3b369936p-56;
        // 0x1.62e42fefa39ef35793c766d326cp-1 = coeff1h + coeff1l
	d_mul(&order1h, &order1l, coeff1h, coeff1l, xh, xl);
	/* Let's expand the d_mul call.
	   Since |coeff1h| < 2^-0.5287 and |xh| <= 2^-19.999 we have
           |ah*bh| <= 2^-20.527. This implies
	   |order1h| < 2^-20.5 and |s| < ulp(2^-20.5) = 2^-73.
	   Since |coeff1l| <= 2^-55.25, we compute that
	     |coeff1l*xh + s| <= 2^-55.25*2^-19.999 + 2^-73 <= 2^-72.724.
	   This ensures that the rounding error computing t is at most
           ulp(2^-72.724) = 2^-125 and that |t| <= 2^-72.7.
	   Since |coeff1h| <= 2^-0.528 and |xl| < 2^-65.408, we compute that
	     |coeff1h*xl + t| <= 2^-0.528 * 2^-65.408 + 2^-72.7 <= 2^-65.922.
	   The rounding error on order1l is therefore at most ulp(2^-65.922)
	   = 2^-118. We also get |order1l| < 2^-65.9.
	   The total rounding error is at most 2^-125 + 2^-118 <= 2^-117.988:

	   The error due to neglecting xl*coeff1l is at most
	     |xl*coeff1l| <= 2^-65.408 * 2^-55.25 <= 2^-120.658
	   The total error on order1 is at most 2^-117.988 + 2^-120.658 <= 2^-117.777:
           |order1h + order1l - (coeff1h + coeff1l) * (xh + xl)| < 2^-117.777.
	*/

	double finalh, finall;
	fast_two_sum(&finalh, &finall, 1, orders23);
	/* Since |orders23| < 2^-42.03, clearly |finalh| < 2, so that
           |finalh| <= 1 + 2^-42.03 + ulp(1) <= 2^0.001.
	   The arguments are in the right order, so this fast_two_sum
           introduces an error at most 2^-105*2^0.001 <= 2^-104.999.
	   Also since |finalh| < 2 we get |finall| < ulp(1) = 2^-52.
	*/

	double tmp;
	fast_two_sum(&finalh, &tmp, finalh, order1h);
	/* At input we have 1/2 < |finalh| < 2^0.001 and |order1h| < 2^-20.5,
           which ensures that the order is respected.
           Furthermore |finalh+order1h| <= 2^0.001 + 2^-20.5 < 2^0.002. This
           ensures that the error is bounded by 2^-105*2^0.002 <= 2^-104.998.
           Also, we get for the new value of finalh:
           |finalh| <= 2^0.002 and |tmp| <= ulp(2^0.002) = 2^-52.
	*/

	finall = tmp + (finall + order1l);
	/* At input, we have |finall| < 2^-52 and |order1l| <= 2^-71.4.
           Therefore the internal sum is strictly bounded by 2^-51.9 and has
           rounding error at most ulp(2^-51.9) = 2^-104.
	   Calling S the intermediate result, since |tmp| <= 2^-52 we have
		 |tmp + S| < 2^-52 + 2^-51.9 < 2^-50.9, which ensures a
                 rounding error of at most ulp(2^-50.9) = 2^-103 and that the
                 new value of finall satisfies:
                 |finall| < 2^-50.9 + 2^-103 <= 2^-50.8.
	   The total rounding error here is 2^-104 + 2^-103 <= 2^-102.415.
	*/

	/* Summing up the errors we get :
	   - 2^-86.439  computing orders23,
	   - 2^-117.777 computing order1,
	   - 2^-104.999 in the first fast_two_sum,
	   - 2^-104.998 in the second fast_two_sum,
	   - 2^-102.415 in the last sum.
	  The polynomial itself was only precise to 2^-89.218. Therefore,
          we have computed 2^xr with error at most :
	   2^-86.439 + 2^-117.777 + 2^-104.999 + 2^-104.998 + 2^-102.415 + 2^-89.218
	   <= 2^-86.2427
	   Since xr >= -2^-19.999, this gives a relative error rho3 less than
	   2^-86.2427/2^(-2^-19.999) < 2^-86.242:
	   finalh + finall = 2^xr * (1 + rho3) with |rho3| < 2^-86.242
	*/
	if(__builtin_expect(do_red,1)) {
		d_mul(&finalh, &finall, finalh, finall, xs_pow2_h, xs_pow2_l);
	  /* We have xs_pow2_h + xs_pow2_l = 2^frac(r) * (1 + rho2)
	     with |rho2| < 2^-99.1, and
	     finalh + finall = 2^xr * (1 + rho3) with |rho3| < 2^-86.242
	     The intrinsic relative error of the product is at most
	     (1 + rho2) * (1 + rho3) - 1 <= 2^-99.1 + 2^-86.242
	                                  + 2^-99.1*2^-86.242 <= 2^-86.241

	     Expanding the d_mul call, we see that:
	     Since |finalh*xs_pow2_h| < 2^0.002*2 <= 2^1.002,
             |s| <= ulp(2^1.002) = 2^-51. Then, since
	     |finall*xs_pow2_h + s| <= 2^-50.8 * 2 + 2^-51 <= 2^-49.278 we get
	     |t| <= 2^-49.277 and that the rounding error on t is at most
             ulp(2^-49.278) = 2^-102.
	     We compute |finalh*xs_pow2_l + t| <= 2^0.002*2^-48.2+2^-49.277
	     <= 2^-47.638. This ensures that the rounding error computing
             *lo is less than ulp(2^-47.638) = 2^-100.
             The total rounding error is therefore at most
	     2^-102 + 2^-100 <= 2^-99.678.
	     The error due to neglecting xs_pow2l * finall is at most
	       |xs_pow2l * finall| <= 2^-48.2 * 2^-50.8 <= 2^-99, thus adding
             it yields an arror < 2^-99.678 + 2^-99 <= 2^-98.299.

	     Since the product should be at least exp(-2^-19.999), this
             translates to an additional relative error
             rho4 <= 2^-98.299/2^(-2^-19.999), so rho4 <= 2^-98.298.
	     Taking into account rho4, the total relative error is thus at most
	      (1 + 2^-86.241)(1 + rho4) - 1 <= 2^-86.240:
             |finalh + finall - 2^frac(r) * 2^xr| < 2^-86.240.
	  */
	} else {
	  /* The only error made is rho3, the total relative error is 2^-86.242. */
		extra_exponent = 0;
	}
	*resh = finalh;
	*resl = finall;
	return extra_exponent;
}

/* Rounding and rounding test for the fastpath*/
inline static
long double fastpath_roundtest(double rh, double rl,int extra_exp,
                               bool invert, bool* fail) {
	unsigned rm = get_rounding_mode();
	b64u64_t th = {.f = rh}, tl = {.f = rl};
	POWL_DPRINTF("rh = %a\nrl = %a\n", rh, rl);
	long eh = th.u>>52, el = (tl.u>>52)&0x3ff, de = eh - el;
	// the high part is always positive, the low part can be positive or negative
	// represent the mantissa of the low part in two's complement format
	long ml = (tl.u&~(0xffful<<52))|1l<<52, sgnl = -(tl.u>>63);
	ml = (ml^sgnl) - sgnl;
	int64_t mlt;
	long sh = de-11;
	if(__builtin_expect(sh>63,0)){
		mlt = sgnl;
		if(__builtin_expect(sh-64>63,0))
			ml = sgnl;
    		else
		ml >>= sh-64;
	} else {
		mlt = ml>>sh;
		ml <<= 64-sh;
	}

	if(__builtin_expect((tl.u&~(1l<<63)) == 0, 0)) {
		// ml == 0. Note that ml cannot be subnormal
		ml = 0; mlt = 0;
	}
	// construct the mantissa of the long double number
	uint64_t mh = ((th.u<<11)|1l<<63);
	/* The evaluation of 2^(yln2(1+x)) was precise to 2^(-86.240) relative error.
	   Since yln2(1+x) was computed to 2^(-92.421) relative error, and must be
	   less than 2^14 lest the result be infinite, yln2(1+x) was known to
	   2^-78.421. This implies an additional relative error of 2^-78.949. The
	   total relative error computing x^y is therefore at most
	   2^-78.949 + 2^-92.421 + 2^-92.421*2^-78.949 = 2^-78.948
	*/
	int64_t eps = (mh >> (78 - 64));
	POWL_DPRINTF("tl_u = %016lx\n", tl.u);
	POWL_DPRINTF("ml = %016lx\n", ml);
	
	mh += mlt;
	if(__builtin_expect(!(mh>>63),0)){ // the low part is negative and
					     // can unset the msb so shift the
					     // number back
		mh = mh<<1 | (uint64_t)ml>>63;
		ml <<= 1;
		extra_exp--;
		eps <<= 1;
	}

	int wanted_exponent = extra_exp + 0x3c00 + eh;
	POWL_DPRINTF("wanted exponent : %x\n", wanted_exponent);
	POWL_DPRINTF("mh||ml = %lx%lx\n", mh, ml);

	if(__builtin_expect(wanted_exponent <= 0, 0)) {
		int shiftby = 1 - wanted_exponent;

		if(__builtin_expect(shiftby == 64, 0)) {
			return 0x1p-16445L * .75L;
		}

		if(__builtin_expect(shiftby > 64, 0)) {	
				return 0x1p-16445L * .25L;
		}
		ml = (uint64_t)ml >> shiftby;
		ml |= mh << (64 - shiftby);
		mh >>= shiftby;	
		eps >>= shiftby;
		wanted_exponent = 0;

		POWL_DPRINTF("Shifting by %u\n", shiftby);
		POWL_DPRINTF("mh||ml = %016lx %016lx\n", mh, ml);
	}

	if(rm==FE_TONEAREST){ // round to nearest
		mh += (uint64_t)ml>>63;
		ml ^= (1ul << 63);
	} else if((rm==FE_UPWARD && !invert) || (rm==FE_DOWNWARD && invert)) {
		mh += 1;
		// This is as if ml had a trailing 1.
		// We are not precise up to an LSB of ml anyway.
	}

	// This branch can only be taken if wanted_exponent != 0
  // Else we simply cannot have an overflow
	if(__builtin_expect(!mh, 0)) {
		ml = ml/2; // Signed semantics matter
	  eps >>= 1;
		mh = 1ull << 63;
	  wanted_exponent++;	
	}

	// We had a denormal but rounding made it into the smallest normal
	if(__builtin_expect((mh>>63) && !wanted_exponent, 0)) {
		wanted_exponent = 1;	
	}

	b80u80_t v;	
	v.m = mh; // mantissa
	v.e = wanted_exponent; // exponent
	if(__builtin_expect(invert, 0)) {v.e += (1<<16);}
	bool b1 = (uint64_t)(ml + eps) <= (uint64_t)(2*eps); 
	*fail = b1;

	// Denormals *inside* the computation don't seem to pause a problem
	// given the error analysis (we used absolute bounds mostly)
	// rounding test

	// Infinity output case
	if(__builtin_expect(wanted_exponent >= 32767, 0)) {
		return 0x1p16383L + 0x1p16383L;
	}
	return v.f;
}

inline static
bool is_integer(long double x) {
	const b80u80_t cvt = {.f = x};
	int e = (cvt.e & 0x7fff) - 16383; // 2^e <= |x| < 2^(e+1)
	if (e >= 63) return true; // ulp(x) >= ulp(2^63) = 1 thus x is integer
        else if (e <= -1) return false; // |x| < 1
        // bit 0 of cvt.m has weight 2^(e-63)
        // thus bit 62-e corresponds to weight 1/2
        // we need the low 63-e bits to equal 000...000
        // now 0 <= e <= 62
        else {
          uint64_t u = cvt.m << (e + 1);
          return u == 0;
        }
}

// return non-zero iff x is an odd integer
inline static
bool is_odd_integer(long double x) {
	const b80u80_t cvt = {.f = x};
        int e = (cvt.e&0x7fff) - 16383; // 2^e <= |x| < 2^(e+1)
	if (e >= 64) return false; // ulp(x) >= ulp(2^64) = 2 thus x is even
        else if (e <= -1) return false; // |x| < 1
        // bit 0 of cvt.m has weight 2^(e-63)
        // thus bit 63-e corresponds to weight 1
        // we need the low 64-e bits to equal 1000...000
        // now 0 <= e <= 63
        else {
          uint64_t u = cvt.m << e;
          return u == (uint64_t) 1 << 63;
        }
}

// return non-zero iff x is a NaN
inline static
int isnan(long double x) {
  const b80u80_t v = {.f = x};
  return ((v.e&0x7fff) == 0x7fff && (v.m != (1ul << 63)));
}

// return non-zero iff x is a signaling NaN
inline static
int issnan(long double x) {
	const b80u80_t v = {.f = x};
	return isnan(x) && (!((v.m>>62)&1));
}

long double cr_powl(long double x, long double y) {

	const b80u80_t cvt_x = {.f = x}, cvt_y = {.f = y};

	if(__builtin_expect(isnan(x) || isnan(y), 0)) {
          if (issnan (x) || issnan (y))
          {
            feraiseexcept(FE_INVALID);
            return __builtin_nanl(""); // Returns a quiet NaN, raises invalid operation
          }
          // qNaN^0 = 1^qNaN = 1
          if (__builtin_expect(cvt_y.m == 0 || x == 1.L, 0)) return 1.L;
          // Return a quiet NaN
          return __builtin_nanl("");
	}  

        // x^0 = 1^y = 1
	if(__builtin_expect(cvt_y.m == 0 || x == 1.L, 0)) return 1.L;

	const int x_exp = (cvt_x.e&0x7fff) - 16383;
        // 2^x_exp <= |x| < 2^(x_exp+1)
	const int y_exp = (cvt_y.e&0x7fff) - 16383;
        // 2^y_exp <= |y| < 2^(y_exp+1)

	bool invert = (cvt_x.e>>15) & is_odd_integer(y);
	const long double sign = ((cvt_x.e>>15) & is_odd_integer(y)) ? -1.L : 1.L;

	static const long double inf = __builtin_infl();	
	if(__builtin_expect(cvt_x.m == 0, 0)) { // x = +- 0
		if(cvt_y.e>>15) { // y < 0
			if(cvt_y.e != 0xffff) feraiseexcept(FE_DIVBYZERO); // If y != -inf
			return sign * inf;
		} else {
			return sign * 0L;
		}
	}

	// -inf < x < 0
	if(__builtin_expect(cvt_x.e >= 1<<15 && cvt_x.e != 0xffff, 0)) {
		if(!is_integer(y)) // Note that +-infty are (even) integers here
			{feraiseexcept(FE_INVALID); return __builtin_nanl("");}
		if(__builtin_expect(x == -1L, 0)) {
			return sign;
		}
	}

	// Now, the handling of forbidden values has been done
	// and the sign of the result computed in sign. We treat x as |x|.

	if(__builtin_expect((cvt_x.e&0x7fff) == 0x7fff, 0)) // x = +-inf
          return sign * ((cvt_y.e>>15) ? 0L : inf);

	bool lt1 = (x_exp < 0) ^ (cvt_y.e>>15);
        // x^y = 2^s with s < 0 (we have already dealt wuth x=1
        // thus x_exp >= 0 implies x > 1)

	if(__builtin_expect(y_exp >= 78, 0)) {
          /* For |y| >= 2^78, since |x| <> 1, the smallest value of
             |y * log2(x)| is attained for x = 1 - 2^-64, and is > 23637,
             thus |x^y| is smaller than the smallest positive subnormal
             2^-16445, or largest than MAX_LDBL = 2^16384*(1-2^-64). */
          // |y| >= 2^78 implies y is a (possibly infinite) even integer
		if(__builtin_expect(y_exp == 0x7fff - 16383, 0)) // y = +-infty
                  return (lt1) ? 0L : inf;
		else
                  return (lt1) ? 0x1p-16445L * .5L : 0x1p16383L + 0x1p16383L;
	} else if(__builtin_expect(y_exp <= -81, 0)) {
          /* For y_exp <= -81, we have |y| < 2^-80,
             thus since |log2(x)| <= 16445, we have |y*log2(x)| < 0x1.00f4p-66.
             Since for |t| <= 0x1.71547652b82fe176p-65, 2^t rounds to 1
             to nearest, we can deduce the correct rounding. */
          return (lt1) ? 1.L - 0x1p-65L : 1.L + 0x1p-65L;
	}

        // now -80 <= y_exp <= 77 thus 2^-80 <= |y| < 2^78

	if(__builtin_expect((int64_t)cvt_x.m < 0, 1)) { // x is normal
		double rh, rl;
		POWL_DPRINTF("x="SAGE_RE"\n",x);
		POWL_DPRINTF("y="SAGE_RE"\n",y);
		compute_log2pow(&rh, &rl, x, y);
    // rh + rl approximates y*log2(x)
		POWL_DPRINTF("get_hex(R(log2(x^y)-"SAGE_DD"))\n",rh,rl);
		long double r;
		if(__builtin_expect(rh <= -16445.2, 0)) {
			return (sign * 0x1p-16445L) * .5L;
		} else if(__builtin_expect(rh >= 16384.5, 0)) {
                  return sign * 0x1p16383L + sign * 0x1p16383L;
	  } else if(__builtin_expect(rh < 0x1p-66 && rh > -0x1p-66, 0)) {
		  return sign * 1.L + sign * rh;
	    // If |rh| is sufficiently small, even with the error margin we know
	    // how x^y rounds.
		} else {
			double resh, resl;
			int extra_exponent = exp2d(&resh, &resl, rh, rl);
			bool fail = false;
			r = fastpath_roundtest(resh, resl, extra_exponent, invert, &fail);
			if(__builtin_expect(fail, 0)) {return -1.;}
		}
		POWL_DPRINTF("get_hex(R(x^y-"SAGE_RE"))\n",r);
		return r;
	} else {return -1.;} // x subnormal
}
