#include <mpfr.h>
#include "fenv_mpfr.h"

/* reference code using MPFR */
long double ref_powl(long double x, long double y) {
  mpfr_t z, _x, _y;
  mpfr_inits2(64, z, _x, _y, NULL);
  mpfr_set_d(_x, x, MPFR_RNDN);
  mpfr_set_d(_y, y, MPFR_RNDN);
  int inex = mpfr_pow(z, _x, _y, rnd2[rnd]);
  mpfr_subnormalize(z, inex, rnd2[rnd]);
  long double ret = mpfr_get_ld(z, rnd2[rnd]);
  mpfr_clears(z, _x, _y, NULL);
  return ret;
}
