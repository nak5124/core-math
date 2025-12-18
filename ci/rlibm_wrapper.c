#include <fenv.h>
#include "rlibm.h"

float log10f (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_log10f (x);
  fesetround (rnd);
  return (float) temp;
}

float log2f (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_log2f (x);
  fesetround (rnd);
  return (float) temp;
}

float logf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_logf (x);
  fesetround (rnd);
  return (float) temp;
}

float exp10f (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_exp10f (x);
  fesetround (rnd);
  return (float) temp;
}

float exp2f (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_exp2f (x);
  fesetround (rnd);
  return (float) temp;
}

float expf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_expf (x);
  fesetround (rnd);
  return (float) temp;
}

float coshf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_coshf (x);
  fesetround (rnd);
  return (float) temp;
}

float sinhf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_sinhf (x);
  fesetround (rnd);
  return (float) temp;
}

float cospif (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_cospif (x);
  fesetround (rnd);
  return (float) temp;
}

float sinpif (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_sinpif (x);
  fesetround (rnd);
  return (float) temp;
}

float sinf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_sinf (x);
  fesetround (rnd);
  return (float) temp;
}

float cosf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_cosf (x);
  fesetround (rnd);
  return (float) temp;
}

float tanf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_tanf (x);
  fesetround (rnd);
  return (float) temp;
}

float atanf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_atanf (x);
  fesetround (rnd);
  return (float) temp;
}

float asinf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_asinf (x);
  fesetround (rnd);
  return (float) temp;
}

float acosf (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_acosf (x);
  fesetround (rnd);
  return (float) temp;
}
