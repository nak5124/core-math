#include <stdio.h>
#include <mpfr.h>
#include <fenv.h>
#include <math.h>
#include "powl_tables.h"

// define RATIO to analyze the maximal value of the ratios
// |mlogrh/mlogr12h|, |mlogr1h/mlogr12h|, |mlogr2h/mlogr12h|
// otherwise analyze the maximal relative error between the
// computed value mlogr12h + mlogr12l and the true value

static inline
void fast_two_sum(double* rh, double* rl, double a, double b) {
	*rh = a + b;
	double e = *rh - a;
	*rl = b - e;
}

/* Computes an approximation of a + bh + bl assuming a = 0 or |a| >= |bh|*/
static inline
void high_sum(double* rh, double* rl, double a, double bh, double bl) {
	double e;
	fast_two_sum(rh, &e, a, bh);
	*rl = bl + e;
}

// return |mlogr2h + mlogr2l - (extra_int + l1[i1] + l2[i2])| / |mlogr2l|
// where l1[i1] = coarse[i1].mlogrh + coarse[i1].mlogrl
// and   l2[i2] = fine[i2].mlogrh + fine[i2].mlogrl
static double
compute_error (int extra_int, int i1, int i2, double mlogr12h, double mlogr12l)
{
  mpfr_t x, y;
  mpfr_init2 (x, 127);
  mpfr_init2 (y, 127);
  // compute in x the sum with 127-bit precision
  mpfr_set_si (x, extra_int, MPFR_RNDN);
  mpfr_add_d (x, x, coarse[i1].mlogrh, MPFR_RNDN);
  mpfr_add_d (x, x, coarse[i1].mlogrl, MPFR_RNDN);
  mpfr_add_d (x, x, fine[i2].mlogrh, MPFR_RNDN);
  mpfr_add_d (x, x, fine[i2].mlogrl, MPFR_RNDN);
  // put in y the value of mlogr2h + mlogr2l
  mpfr_set_d (y, mlogr12h, MPFR_RNDN);
  mpfr_add_d (y, y, mlogr12h, MPFR_RNDN);
  mpfr_sub (x, x, y, MPFR_RNDN);
  mpfr_abs (x, x, MPFR_RNDN);
  double err = mpfr_get_d (x, MPFR_RNDN);
  mpfr_clear (x);
  mpfr_clear (y);
  return err / mlogr12h;
}

// find by brute-force the maximal relative error of the 2nd high_sum call
// in compute_log2pow()
static void
analyze_second_sum (void)
{
  double maxerr = 0, maxratio1 = 0, maxratio2 = 0, maxratio3 = 0;
  static int R[4] = {FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD};
  for (int r = 0; r < 4; r++)
  {
    fesetround (R[r]);
    for (int extra_int = -16382; extra_int < 16384; extra_int ++)
      for (int i1 = 0; i1 < 128; i1++)
      {
        double mlogr1h = coarse[i1].mlogrh;
        double mlogr1l = coarse[i1].mlogrl;
        double mlogrh, mlogrl;
        high_sum(&mlogrh, &mlogrl, (double) extra_int, mlogr1h, mlogr1l);
        for (int i2 = 0; i2 < 128; i2++)
        {
          if (32 <= i2 && i2 < 64)
            continue;
          double mlogr2h = fine[i2].mlogrh;
          double mlogr2l = fine[i2].mlogrl;
          // printf ("%d %d %d %d\n", r, extra_int, i1, i2);
          double mlogr12h, mlogr12l;
          high_sum(&mlogr12h, &mlogr12l, mlogrh, mlogr2h, mlogr2l);
          mlogr12l += mlogrl;
#ifdef RATIO
          double ratio1 = fabs (mlogrh / mlogr12h);
          double ratio2 = fabs (mlogr2h / mlogr12h);
          double ratio3 = fabs (mlogr1h / mlogr12h);
          if (ratio1 > maxratio1)
            printf ("r=%d extra_int=%d i1=%d i2=%d |mlogrh/mlogr12h|=%.16e\n",
                    r, extra_int, i1, i2, maxratio1 = ratio1);
          if (ratio2 > maxratio2)
            printf ("r=%d extra_int=%d i1=%d i2=%d |mlogr2h/mlogr12h|=%.16e\n",
                    r, extra_int, i1, i2, maxratio2 = ratio2);
          if (ratio3 > maxratio3)
            printf ("r=%d extra_int=%d i1=%d i2=%d |mlogr1h/mlogr12h|=%.16e\n",
                    r, extra_int, i1, i2, maxratio3 = ratio3);
          fflush (stdout);
#else
          double err = compute_error (extra_int, i1, i2, mlogr12h, mlogr12l);
          if (err > maxerr)
            printf ("r=%d extra_int=%d i1=%d i2=%d err=%f\n",
                    r, extra_int, i1, i2, maxerr = err);
#endif
        }
      }
  } 
}

int
main ()
{
  analyze_second_sum ();
  return 0;
}
