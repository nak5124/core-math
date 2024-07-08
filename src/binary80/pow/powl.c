#include <stdint.h>
#include <fenv.h>
#include <stdbool.h>

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

/* Split a number of exponent 0 into a w part on 11 bits and an a part
	 on 53 bits.
*/
static inline
void split(double* w, double* a, long double x) {
	static const long double C = 0x1.8p+53L; // ulp(C) = 2^-10
	long double y = (x + C) - C; // we get in w the bits of weight 1,...,2^-10
	*w = y;
	*a = x - y; 
}

static inline
void add22(double* zh, double* zl, double xh, double xl, double yh, double yl) {
	double r,s;
	r = xh+yh;
	s = ((xh-r)+yh)+yl+xl;
	*zh = r+s;
	*zl = (r - (*zh)) + s;
}

/* Computes an approximation of a/(bh+bl). Latency ~32 cycles. */
static inline
void s_div(double* rh, double* rl, double a, double bh, double bl) {
	*rh = a/bh;
	double eh = __builtin_fma(-bh, *rh, a);
	*rl = __builtin_fma(-bl, *rh, eh)/bh;
}

static inline
void fast_two_sum(double* rh, double* rl, double a, double b) {
	*rh = a + b;
	double e = *rh - a;
	*rl = b - e;
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


/* Computes an approximation of ylog2|x| under the following conditions :

- 2^(-65-15) <= |y| < 2^(65 + 14)
- |x| is at least the smallest positive normal number

If |ylog2(x)| > 16383, one of the outputs may be +-inf.
*/
static inline
void compute_log2pow(double* rh, double* rl, long double x, long double y) {

	b80u80_t cvt_x = {.f = x};
	int extra_int = (cvt_x.e&0x7fff) - 16383; // We don't have to mask : x >= +0
	cvt_x.e = 16383; // New wanted exponent
	x = cvt_x.f;

	POWL_DPRINTF("extra_int=%d\n", extra_int);
	// This may not be efficient, but we avoid overflow/underflow problems.
	// Note that x starts in memory so this should not be too expensive.
	double w,a; split(&w, &a, x); // Note that w can be exactly 2

	// Note that on modern processors, an fma is as expensive as an add.
	double denomh, denoml;
	denomh = __builtin_fma(2,w,a);
	double e = __builtin_fma(2,w,-denomh); // -a - deltah
	denoml = e + a;

	//fast_two_sum(&denomh, &denoml, 2*w, a);
	POWL_DPRINTF("w="SAGE_RR"\n",w);
	POWL_DPRINTF("a="SAGE_RR"\n",a);
	POWL_DPRINTF("get_hex(w + a - "SAGE_RE")\n", x);

	b64u64_t cvt_w = {.f = w};
	double logwh, logwl;

	if(__builtin_expect(w == 2., 0)) {logwh = extra_int + 1; logwl = 0;}
	else {
		double ex;
		logwh = t[(cvt_w.u>>42) & 0x3ff][0]; logwl = t[(cvt_w.u>>42) & 0x3ff][1];
		fast_two_sum(&logwh, &ex, extra_int, logwh);
		logwl += ex;
	}
	/* We could construct the tables in such a way that logwh + extra_int is
	   exact. However, this would waste 14 bits of precision when extra_int is
	   0, and even more for small values of logw. This chain of operation is
	   not parallelisable and has latency ~16cycles. We can expect this is hidden
	   in the division latency.
	*/
	POWL_DPRINTF("get_hex(R(log2(w*2^%d) - "SAGE_DD"))\n",extra_int,logwh,logwl);
	
	double th, tl; s_div(&th, &tl, a, denomh, denoml);
	POWL_DPRINTF("t="SAGE_DD"\n", th, tl);
	POWL_DPRINTF("get_hex(R(t - a/(2*w + a)))\n");

	/* Evaluation strategy :
	   2/(ln 2) * t (1 + t^2/3 + t^4/5 + t^6/7)
	   = (2/3ln2) * t * (3 + t^2 + 3/5t^4 + 3/7t^6)
	*/
	double tsqh, tsql; d_mul(&tsqh, &tsql, th, tl, th, tl);

	double highorder = tsqh * tsqh;
	highorder *= __builtin_fma(tsqh, 0x1.b6db7860d5f78p-2, 0x1.333333333325ep-1);

	double lorderh, lorderl;
	fast_two_sum(&lorderh, &lorderl, 0x1.8p+1, tsqh);
	lorderl += tsql + (-0x1.00c47fd09p-68); // I think we still are normalized
	//add22(&lorderh, &lorderl, 0x1.8p+1, -0x1.dc61e618p-75, tsqh, tsql);

	double scaleh, scalel;
	d_mul(&scaleh, &scalel, 0x1.ec709dc3a03fdp-1, 0x1.d28197e2ad4ccp-55, th, tl);
	lorderl += highorder; /* Highorder <~= 2^-48 */
	d_mul(&lorderh, &lorderl, lorderh, lorderl, scaleh, scalel);

	POWL_DPRINTF("get_hex(R(2*atanh(t)/ln(2) - "SAGE_DD"))\n", lorderh, lorderl);

	double yh = y; double yl = y - (long double)(yh);
	fast_two_sum(rh, rl, logwh, logwl);
	*rl += lorderh + lorderl;
	//add22(rh, rl, logwh, logwl, lorderh, lorderl);
	d_mul(rh, rl, *rh, *rl, yh, yl);
}

static inline
long double exp2d(double xh, double xl) {
	b64u64_t cvt = {.f = xh};
	bool do_red	= cvt.e >= -20 + 0x3ff;

	static const double C = 0x1.8p+32;
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
	v.e = wanted_exponent; //extra_exponent + 0x3c00 + eh; // exponent
	
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

static
bool is_integer(long double x) {
	b80u80_t cvt = {.f = x};
	int e = (cvt.e & 0x7fff) - 16383;
	if(e >= 63) { // Ulp is 2^(e - 63) >= 1
		return true;
	} else if(e >= -1) {
		return !(cvt.m & (-1ul >> (e + 1)));
	} else {return false;}
	// low bits must be 0
}

static bool is_odd_integer(long double x) {
	b80u80_t cvt = {.f = x};
	if((cvt.e&0x7fff) - 16383 >= 64) return false;
	else return is_integer(x) && (cvt.m & (1ul << (63 - (cvt.e&0x7fff) + 16383)));
}

inline
static int isnan(long double x) {
  b80u80_t v = {.f = x};
  return ((v.e&0x7fff) == 0x7fff && (v.m != (1ul << 63)));
}

inline
static int issnan(long double x) {
	b80u80_t v = {.f = x};
	return isnan(x) && (!((v.m>>62)&1));
}

long double cr_powl(long double x, long double y) {

	const b80u80_t cvt_x = {.f = x}, cvt_y = {.f = y};
	if(__builtin_expect(issnan(x) || issnan(y), 0)) {
		return x + y; // Returns a quiet NaN, raises invalid operation
	}  

	if(__builtin_expect(cvt_y.m == 0, 0)) return 1.L;
	if(__builtin_expect(cvt_x.m == 0x8000000000000000ul
		&& cvt_x.e == 16383, 0)) return 1.L;

	int x_exp = (cvt_x.e&0x7fff) - 16383;
	long double sign = ((cvt_x.e>>15) & is_odd_integer(y)) ? -1.L : 1.L;

	if(__builtin_expect(isnan(x), 0)) return x + x; // Check for NaN, quiet it.
	if(__builtin_expect(isnan(y), 0)) return y + y;

	static long double inf = __builtin_infl();	
	if(__builtin_expect(cvt_x.m == 0, 0)) { // x = +- 0
		if(cvt_y.e>>15) { // Need to raise divide_by_zero if odd_integer here
			if(is_odd_integer(y)) {return sign * 1./0.;}
			else { return sign * inf; }
		} else {
			return sign * 0L;
		}
	}

	// -inf < x < 0
	if(__builtin_expect(cvt_x.e >= 1<<15 && cvt_x.e != 0xffff, 0)) {
		if(!is_integer(y)) {return 0.L/0.L;} // Raises invalid exception
		if(__builtin_expect(x == -1L, 0)) {
			return sign;
		}
	}

	// Now, the handling of forbidden values has been done
	// and the sign of the result computed in sign. We treat x as |x|.
	
	if((x_exp < 0) ^ (cvt_y.e >> 15)) { // 2^s with s < 0
		if(__builtin_expect((cvt_y.e&0x7fff) == 0x7fff, 0) ||
		   __builtin_expect((cvt_x.e&0x7fff) == 0x7fff, 0)) { // s == +-inf
			return sign * 0L;
		} else if(__builtin_expect((cvt_y.e&0x7fff) - 16383 >= 79, 0)) {
			return sign * 0x1p-16445L * .5L;
		}
	} else { // 2^s with s > 0
		if(__builtin_expect((cvt_y.e&0x7fff) == 0x7fff, 0) ||
	     __builtin_expect((cvt_x.e&0x7fff) == 0x7fff, 0)) { // s == +-inf
			return sign * inf;
		} else if(__builtin_expect((cvt_y.e&0x7fff) - 16383 >= 79, 0)) {
			return sign * (0x1p16383L + 0x1p16383L);
		}
	}

	// Automatic giveup if x subnormal
	if(__builtin_expect((int64_t)cvt_x.m < 0, 1)) {
		double rh, rl;
		POWL_DPRINTF("x="SAGE_RE"\n",x);
		POWL_DPRINTF("y="SAGE_RE"\n",y);
		compute_log2pow(&rh, &rl, x, y);
		POWL_DPRINTF("get_hex(R(log2(x^y)-"SAGE_DD"))\n",rh,rl);
		long double r = exp2d(rh, rl);
		POWL_DPRINTF("get_hex(R(x^y-"SAGE_RE"))\n",r);
		return sign * r;
	} else {return -1.;}
}
