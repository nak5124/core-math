#include <stdint.h>

#ifdef POWL_DEBUG
#include <stdio.h>
#define POWL_DPRINTF(...) printf(__VA_ARGS__)
#define SAGE_RR "R(\"%a\",16)"
#define SAGE_RE "R(\"%La\",16)"
#define SAGE_DD "(R(\"%a\",16)+R(\"%a\",16))"
#else
#define POWL_DPRINTF(...)
#endif

typedef union {long double f; struct {uint64_t m; uint16_t e;};} b80u80_t;
typedef union {double f; uint64_t u;} b64u64_t;

/* Split a number of exponent 0 into a w part on 11 bits and an a part
	 on 53 bits.
*/
static inline
void split(double* w, double* a, long double x) {
	static const double C = 0x1.8p+42; // ulp(C) = 2^-10
	*w = ((double)x + C) - C; // we get in w the bits of weight 1,...,2^-10
	*a = (double)x - *w; 
}

static inline
void add22(double* zh, double* zl, double xh, double xl, double yh, double yl) {
	double r,s;
	r = xh+yh;
	s = ((xh-r)+yh)+yl+xl;
	*zh = r+s;
	*zl = (r - (*zh)) + s;
}
/* Computes an approximation of a/b. Latency ~32 cycles. */
static inline
void s_div(double* rh, double* rl, double a, double b) {
	*rh = a/b;
	*rl = __builtin_fma(-b, *rh, a)/b;
}

static inline
void s_long_div(double* rh, double* rl, long double a, long double b) {
	long double k = a/b;

	// The result is correct to 2^-64 bits.
	double kh = k; double kl = k - (long double)kh;
	// This goes back to SSE I think
	double err = __builtin_fma(kl, (double)b, (double)(kh*b - a));
	kl -= err/b;

	*rh = kh;
	*rl = kl;
}

static inline
void a_mul(double* rh, double* rl, double a, double b) {
	*rh = a*b;
	*rl = __builtin_fma(a,b,-*rh);
}

static inline
void d_mul(double* rh, double* rl, double ah, double al,
                                   double bh, double bl)
{ double s;
	a_mul(rh, &s, ah, bh);
	double t = __builtin_fma(al, bh, s);
	*rl = __builtin_fma(ah, bl, t);
}

#include "powl_tables.h"

static inline
void compute_log2pow(double* rh, double* rl, long double x, long double y) {
	b80u80_t cvt_x = {.f = x};
	int extra_int = (cvt_x.e&0x7fff) - 16383;
	POWL_DPRINTF("extra_int=%d\n", extra_int);
	// This may not be efficient, but we avoid overflow/underflow problems.
	// Note that x starts in memory so this should not be too expensive.

	x = __builtin_ldexpl(1, -extra_int) * x; // Scale x
	double w,a; split(&w, &a, x);
	POWL_DPRINTF("w="SAGE_RR"\n",w);

	b64u64_t cvt_w = {.f = w};
	double logwh, logwl;
	logwh = t[(cvt_w.u>>42) & 0x3ff][0]; logwl = t[(cvt_w.u>>42) & 0x3ff][1];
	logwh += (double)extra_int; // exact

	/* To test : if we can shoulder an 80-bit division,
	   we can use the improved polynomial. This comes at a ~20 cycle penalty.
	   OTOH, directly using a/w means we have to sum the first 3 terms as dds
	   but makes a better use of pipelining.

	   Sollya gives a polynomial for log2(1+x)/x accurate to 2^-89.5.
	   This translates to a final relative error bound which gives 2^-10
	   probability to break out of the fastpath.
	   If we use the improved polynomial, we can aim for ~10 more bits of
	   precision OR the same precision but tables 2x smaller LUTs.

	   TODO Use Tang's algorithm ? May save on computation time.
	*/
	double th, tl; s_div(&th, &tl, a, w);
	double t2h, t2l; d_mul(&t2h, &t2l, th, tl, th, tl);
	double t4 = t2h * t2h;

	// Addition and fma throughput is .5 CPI. We can expect order01 and order23's
	// first two ops to run in parallel, with total latency 4*6 + 16 = 40
	double ord01h, ord01l; d_mul(&ord01h, &ord01l, th, tl,
	                             0x1.71547652b82fep+0, 0x1.777d0ffd4eca8p-56);
	add22(&ord01h, &ord01l, logwh, logwl, ord01h, ord01l);

	double ord23h, ord23l; d_mul(&ord23h, &ord23l, th, tl,
	                             -0x1.71547652b82fep-1, -0x1.777d11aa72848p-57);
	add22(&ord23h, &ord23l, 0x1.ec709dc3a03fdp-2, 0x1.d3231f5394554p-56,
	                        ord23h, ord23l);

	// We can expect the products here to execute concurrently as the high order
	// fmas.
	d_mul(&ord23h, &ord23l, t2h, t2l, ord23h, ord23l);
	add22(&ord01h, &ord01l, ord01h, ord01l, ord23h, ord23l);

	double acc = __builtin_fma(th,0x1.a6178bb07904cp-3,-0x1.ec709dc3f6604p-3);
	double bcc = __builtin_fma(th,0x1.2776c50ef8f2cp-2,-0x1.71547652b82fep-2);
	acc =  __builtin_fma(t2h, acc, bcc);
	acc *= t4;

	ord01l = __builtin_fma(t4, acc, ord01l);
	ord01l = __builtin_fma(t4, t4*-0x1.7151a1efb867p-3, ord01l);
	// Maybe we need a Fast2Sum here

	/*If y rounds to +-infty as a double, then, given that
		|log2(x)| >~ 2^-64 for x != 1, we had in fact |ylog2(x)| >> 16383.
	*/
	double yh = y; double yl = y - (long double)yh;
	d_mul(&ord01h, &ord01l, yh, yl, ord01h, ord01l);
	*rh = ord01h; *rl = ord01l;
}

//static inline
//void exp2(double* rh, double* rl, double xh, double xl) {
	/* More or less the same as exp's implementation */
//}

long double cr_powl(long double x, long double y) {
	double rh, rl;
	POWL_DPRINTF("x="SAGE_RE"\n",x);
	POWL_DPRINTF("y="SAGE_RE"\n",y);
	compute_log2pow(&rh, &rl, x, y);
	POWL_DPRINTF("get_hex(R(log2(x^y)-"SAGE_DD"))\n",rh,rl);
	return rl;
}
