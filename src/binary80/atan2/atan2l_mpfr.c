/* Correctly-rounded reference arc tangent function of two binary80 variables.

Copyright (c) 2026 Alexei Sibidanov <sibid@uvic.ca>.

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

#include <mpfr.h>
#include "fenv_mpfr.h"

/* reference code using MPFR */
long double ref_atan2l(long double y, long double x) {
  mpfr_t z, _x, _y;
  mpfr_exp_t emin = mpfr_get_emin ();
  mpfr_set_emin (-16444);
  mpfr_inits2(64, z, _x, _y, NULL);

  mpfr_set_ld(_x, x, MPFR_RNDN);
  mpfr_set_ld(_y, y, MPFR_RNDN);
  int inex = mpfr_atan2 (z, _y, _x, rnd2[rnd]);
  mpfr_subnormalize(z, inex, rnd2[rnd]);
  long double ret = mpfr_get_ld(z, rnd2[rnd]);

  mpfr_clears(z, _x, _y, NULL);
  mpfr_set_emin (emin);
  return ret;
}
