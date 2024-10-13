#include <mpfr.h>
#include "fenv_mpfr.h"
#include <stdint.h>

#include "cm_types.h"

inline
static int is_nan(long double x) {
  b80u80_t v = {.f = x};
  return ((v.e&0x7fff) == 0x7fff && (v.m != (1ul << 63)));
}

inline
static int issnan(long double x) {
	b80u80_t v = {.f = x};
	return is_nan(x) && (!((v.m>>62)&1));
}

/* reference code using MPFR */
long double ref_powl(long double x, long double y) {
	if(issnan(x) || issnan(y)) {
		b80u80_t v;
		v.m = 0x8000000000000001ul;
		v.e = 0x7fff;
		return v.f;
	}
	mpfr_t z, _x, _y;
	mpfr_exp_t emin = mpfr_get_emin ();
  mpfr_set_emin (-16444);
  
	mpfr_inits2(64, z, _x, _y, NULL);
  mpfr_set_ld(_x, x, MPFR_RNDN);
  mpfr_set_ld(_y, y, MPFR_RNDN);
  int inex = mpfr_pow(z, _x, _y, rnd2[rnd]);
  mpfr_subnormalize(z, inex, rnd2[rnd]);
  long double ret = mpfr_get_ld(z, rnd2[rnd]);

  mpfr_clears(z, _x, _y, NULL);
  mpfr_set_emin (emin);
  return ret;
}
