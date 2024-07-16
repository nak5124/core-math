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

/* Split a number of exponent 0 into a high part on 34 bits an a low part on
31 bits exactly.
*/
static inline
void split(double* rh, double* rl, long double x) {
	static double C = 0x1.8p+31; // ulp(C)=2^-32 once cast to 80 bits
	long double y = (x + C) - C;
	/* Given the relative sizes of C and x, x + C has the binade of C. Therefore,
	   the difference is exact. Furthermore, ulp(x + C) = ulp(C) = 2^-63*2^31
	   = 2^-32. The rounding error in x + C is therefore less than 2^-32.
	   Thus, |x - y| < 2^-32. Note that since 2^31 <= x + C < 2^32 and the
	   difference is exact, y is a multiple of ulp(x + C) = 2^-32.
	   Since |y| < 4 given |x - y|, this ensures y fits in
	   32 - (-1) + 1 = 34 bits.
	*/
	*rh = y; // This conversion is exact by the argument above.
	*rl = x - y;
	/* 
		|x - y| < 2^-32. Note that x and y are both multiples of 2^-63; therefore
	  x - y is too. This implies that x - y fits in 63 - 33 + 1 = 31
	  mantissa bits and that the difference is exact.
	*/
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
	// FIXME directly analyze the absolute errors.
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

- 2^(-65-15) <= |y| < 2^(65 + 14)
- |x| is at least the smallest positive normal number
*/
static inline
void compute_log2pow(double* rh, double* rl, long double x, long double y) {
	b80u80_t cvt_x = {.f = x};
	int extra_int = (cvt_x.e&0x7fff) - 16383;
	cvt_x.e = 16383; // New wanted exponent
	x = cvt_x.f;

	double xh, xl; // a (resp b) bits
	split(&xh, &xl, x);

	POWL_DPRINTF("sx = " SAGE_RE "\nei = %d\n", x, extra_int);
	// Uses the high 7 bits of x's mantissa.
	lut_t l = coarse[cvt_x.m>>56 & 0x7f];
	POWL_DPRINTF("key=0x%lx\n", (cvt_x.m>>56 & 0x7f));

	/* If l.z is 1, then |x*r1 - 1| <= 0x1p-7.
	   If l.z is 0, then |(x/2)*r1 - 1| <= 0x1p-7.
	   In all cases, |mlogr1h + mlogr1l approximates log2(r1) with relative error
	   at most 2^-107. Note that |mlogr1h+mlogr1l| <= .505.
	   The tables are constructed in such a way that r fits in 9 mantissa bits.
	*/
	double r1      = l.r;
	double mlogr1h = l.mlogrh;
	double mlogr1l = l.mlogrl;
	extra_int     += l.z;

	POWL_DPRINTF("r1 = " SAGE_RR "\n", r1);
	
	if(l.z) {xh/=2; xl/=2; POWL_DPRINTF("sx = sx/2\nei+=1\n");}
	xh *= r1; xl *= r1; // xh fits in 43 bits at most, xl in 40.

	POWL_DPRINTF("get_hex(R(abs("SAGE_RR" - 1)))\n", xh);
	/* Note that now |xh - 1| <= 1p-7 (say)
	   Therefore, xh's mantissa is either
	   1.00000 00p or 1.11111 1q.
	   We skip the first 6 bits of the mantissa and use the 7 next to index
	   another lookup table. A quarter of the table is wasted!

		 We're looking at 1 + 2^-(5+7)k, 1 + 2^-(5+7)(k+1) when the high bit of k is
	   0. Else we're looking at 1 - 2^-7 + 2^(-6-7)*k', 1 - 2^-6 + 2^-13(k'+1)
	*/

	b64u64_t cvt_xh = {.f = xh};
	lut_t l2 = fine[(cvt_xh.u>>40) & 0x7f];
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
	        relative error at most 2^-107.
	   iv) |mlogr2h+mlogr2l| <= 2^-8.
	*/

	double mlogrh, mlogrl;
	high_sum(&mlogrh, &mlogrl, extra_int, mlogr1h, mlogr1l);
	//add22(&mlogrh, &mlogrl, mlogr1h, mlogr1l, mlogr2h, mlogr2l);
	POWL_DPRINTF("get_hex(R(-log2(r1)-" SAGE_DD"))\n", mlogr1h, mlogr1l);
	POWL_DPRINTF("get_hex(R(-log2(r2)-" SAGE_DD"))\n", mlogr2h, mlogr2l);

	//high_sum(&mlogrh, &mlogrl, scaled_eint, mlogrh, mlogrl);
	add22(&mlogrh, &mlogrl, mlogrh, mlogrl, mlogr2h, mlogr2l);
	/* If extra_int was not 0, then mlogrh + mlogrl >= 1/2 which implies that
	   |mlogrh| >= |mlogr2h|. If extra_int was 0, then either mlogr1 = 0 (in which
	   case everything is exact) or |mlogr1h| >= 0x1.6fe50b6ef0851p-7 >= |mlogr2h|
	   This ensures that the arguments are in the right order.

	   Sketch : If |mlogr| >= 2|mlogr2| the analysis goes well (this is always
	   true when extra_int != 0).  If this is not the case then we can only be
	   in a few different cases for mlogr1, mlogr2, which we need to manually
	   check (a few bits get cancelled).
	*/
	POWL_DPRINTF("get_hex(R(-log2(r1) - log2(r2)+ei- "SAGE_DD"))\n",
		mlogrh, mlogrl);

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
	   this may be enough ?
	*/
	POWL_DPRINTF("get_hex(R(r1*r2*sx - 1 - "SAGE_DD"))\n", xh, xl);
	POWL_DPRINTF("s = r1*r2*sx - 1\n");
	POWL_DPRINTF("get_hex(s)\n");

	polyeval(rh, rl, xh, xl);

	double yh = y; double yl = y - (long double)(yh);
	/* If mlogr != 0 then |mlogrh| >= 1.6p-12 >= 1.01p-12 >= |rh|
	*/
	add22(rh, rl, mlogrh, mlogrl, *rh, *rl);

	POWL_DPRINTF("get_hex(R(log2(x)) - "SAGE_DD")\n", *rh, *rl);
	d_mul(rh, rl, *rh, *rl, yh, yl);

}


/* computes 2^(xh + xl), assuming |xl| <= ulp(xh) and |xh| < 2^31 */
static inline
long double exp2d(double xh, double xl) {
	b64u64_t cvt = {.f = xh};
	bool do_red	= cvt.e >= -20 + 0x3ff;

	static const double C = 0x1.8p+32; // ulp is 2^-20
	b64u64_t y = {.f = xh + C};
	uint64_t fracpart = y.u;
	int16_t extra_exponent = y.u>>20;

	double rem = xh - (y.f - C);
	if(__builtin_expect(do_red, 1))
		fast_two_sum(&xh,&xl,rem,xl);

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

	double xsq = xh * xh;
	double orders23 = xsq * __builtin_fma(xh,0x1.c6b08d704a1cdp-5,
		0x1.ebfbdff82c696p-3);

	double order1h, order1l;
	static const double coeff1h = 0x1.62e42fefa39efp-1;
	static const double coeff1l = 0x1.abc9e3b369936p-56;
	d_mul(&order1h, &order1l, xh, xl, coeff1h, coeff1l);

	double finalh, finall;
	fast_two_sum(&finalh, &finall, 1, orders23);

	double tmp;
	fast_two_sum(&finalh, &tmp, finalh, order1h);

	finall = tmp + (finall + order1l);

	if(__builtin_expect(do_red,1)) {
		d_mul(&finalh, &finall, finalh, finall, xs_pow2_h, xs_pow2_l);
	} else {
		extra_exponent = 0;
	}

	const unsigned rm = get_rounding_mode();
	b64u64_t th = {.f = finalh}, tl = {.f = finall};
	long eh = th.u>>52, el = (tl.u>>52)&0x3ff, de = eh - el;
	// the high part is always positive, the low part can be positive or negative
	// represent the mantissa of the low part in two's complement format
	long ml = (tl.u&~(0xfffl<<52))|1l<<52, sgnl = -(tl.u>>63);
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
	// construct the mantissa of the long double number
	uint64_t mh = ((th.u<<11)|1l<<63);
	int64_t eps = (mh >> (87 - 64));
	
	mh += mlt;
	if(__builtin_expect(!(mh>>63),0)){ // the low part is negative and
					     // can unset the msb so shift the
					     // number back
		mh = mh<<1 | (uint64_t)ml>>63;
		ml <<= 1;
		extra_exponent--;
		eps <<= 1;
	}

	int wanted_exponent = extra_exponent + 0x3c00 + eh;
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
		POWL_DPRINTF("mh||ml = %lx%lx\n", mh, ml);
	}

	if(rm==FE_TONEAREST){ // round to nearest
		mh += (uint64_t)ml>>63;
		ml ^= (1ul << 63);
	} else if(rm==FE_UPWARD) { // round to +inf
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
	
	/*bool b1 = false;//(uint64_t)(ml + eps) <= (uint64_t)(2*eps); 

	// Denormals *inside* the computation don't seem to pause a problem
	// given the error analysis (we used absolute bounds mostly)
	// rounding test
  *is_accurate = b1;
	if(__builtin_expect(b1, 0)) {
		ri->fracpart = fracpart;
		ri->xs = y.f - C; // Being careful with the double rounding, we've done that
		ri->extra_exponent = extra_exponent;
	}*/

	// Infinity output case
	if(__builtin_expect(wanted_exponent >= 32767, 0)) {
		return 0x1p16383L + 0x1p16383L;
	}
	return v.f;
}

inline static
bool is_integer(long double x) {
	const b80u80_t cvt = {.f = x};
	int e = (cvt.e & 0x7fff) - 16383;
	if(e >= 63) { // Ulp is 2^(e - 63) >= 1
		return true;
	} else if(e >= -1) {
		return !(cvt.m & (-1ul >> (e + 1)));
	} else {return false;}
	// low bits must be 0
}

inline static
int is_odd_integer(long double x) {
	const b80u80_t cvt = {.f = x};
	if((cvt.e&0x7fff) - 16383 >= 64) return false;
	else return is_integer(x) && (cvt.m & (1ul << (63 - (cvt.e&0x7fff) + 16383)));
}

inline static
int isnan(long double x) {
  const b80u80_t v = {.f = x};
  return ((v.e&0x7fff) == 0x7fff && (v.m != (1ul << 63)));
}

inline static
int issnan(long double x) {
	const b80u80_t v = {.f = x};
	return isnan(x) && (!((v.m>>62)&1));
}

long double cr_powl(long double x, long double y) {

	const b80u80_t cvt_x = {.f = x}, cvt_y = {.f = y};
	if(__builtin_expect(issnan(x) || issnan(y), 0)) {
		feraiseexcept(FE_INVALID);
		return __builtin_nanl(""); // Returns a quiet NaN, raises invalid operation
	}  

	if(__builtin_expect(cvt_y.m == 0 || x == 1.L, 0)) return 1.L;

	const int x_exp = (cvt_x.e&0x7fff) - 16383;
	const int y_exp = (cvt_y.e&0x7fff) - 16383;
	const long double sign = ((cvt_x.e>>15) & is_odd_integer(y)) ? -1.L : 1.L;

	// Return a quiet NaN
	if(__builtin_expect(isnan(x) || isnan(y), 0)) return __builtin_nanl("");

	static const long double inf = __builtin_infl();	
	if(__builtin_expect(cvt_x.m == 0, 0)) { // x = +- 0
		if(cvt_y.e>>15) { // y < 0
			if(cvt_y.e != 0xffff) {feraiseexcept(FE_DIVBYZERO);} // If y != -inf
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

	if(__builtin_expect((cvt_x.e&0x7fff) == 0x7fff, 0)) { // x = +-inf
		if(cvt_y.e>>15) {return sign * 0L;} else {return sign * inf;}
	}

	bool lt1 = (x_exp < 0) ^ (cvt_y.e>>15); // x^y = 2^s with s < 0
	// If y is that big, necessarily |yln2(x)| >= 2^15
	// Note that sign == 1 here because y is a (possibly infinite) even integer
	if(__builtin_expect(y_exp >= 79, 0)) {
		if(__builtin_expect(y_exp == 0x7fff - 16383, 0)) { // y = +-infty
			if(lt1) {return 0L;}
			else {return inf;}
		} else {
			if(lt1) {return 0x1p-16445L * .5L;}
			else { return 0x1p16383L + 0x1p16383L;}
		}
	} else if(__builtin_expect(y_exp <= -81, 0)) {
		if(lt1) {return 1.L - 0x1p-65L;}
		else    {return 1.L + 0x1p-65L;}
	}

	/* Note that log2|x| < 2^15. Therefore, if |y| < 2^-80 we have
	   |ylog2|x|| <= 2^-65, which ensures that 2^(ylog2|x|) rounds the same way as
	   1 + sgn(ylog2|x|) * 1p-16445L.
	*/

	// Automatic giveup if x subnormal
	if(__builtin_expect((int64_t)cvt_x.m < 0, 1)) {
		double rh, rl;
		POWL_DPRINTF("x="SAGE_RE"\n",x);
		POWL_DPRINTF("y="SAGE_RE"\n",y);
		compute_log2pow(&rh, &rl, x, y);
		POWL_DPRINTF("get_hex(R(log2(x^y)-"SAGE_DD"))\n",rh,rl);
		long double r;
		if(__builtin_expect(rh <= -16446, 0)) {
			return sign * 0x1p-16445L * .5L;
		} else if(__builtin_expect(rh >= 16383.5, 0)) {
			return sign * 0x1p16383L + sign * 0x1p16383L;
		} else {
			// TODO Directed roundings change depending on sign in exp2d
			r = exp2d(rh, rl);
		}
		POWL_DPRINTF("get_hex(R(x^y-"SAGE_RE"))\n",r);
		return r;
	} else {return -1.;}
}
