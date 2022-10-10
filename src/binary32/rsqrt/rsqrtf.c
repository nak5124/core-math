#include <stdint.h>
#include <errno.h>
#include <fenv.h>

typedef union {float f; uint32_t u;} b32u32_u;

float cr_rsqrtf(float x){
  double xd = x;
  b32u32_u ix = {.f = x};
  if(__builtin_expect(ix.u >= 0xff<<23 || ix.u==0, 0)){
    if(!(ix.u << 1)) return 1.0f/x;
    if(ix.u >> 31){
      ix.u &= ~0u>>1;
      if(ix.u > 0xff<<23) return x;
      errno = EDOM;
      feraiseexcept (FE_INVALID);
      return __builtin_nanf("<0");
    }
    if(!(ix.u<<9)) return 0.0f;
    return x;
  }
  double r = (1.0/xd)*__builtin_sqrt(xd);
  double rx = r*xd, drx = __builtin_fma(r, xd, -rx);
  double h = __builtin_fma(r,rx,-1.0) + r*drx, dr = (r*0.5)*h;
  double rf = r - dr;
  return rf;
}

/* rsqrt function is not in glibc so define it here just to compile tests */
float rsqrtf(float x){
  return cr_rsqrtf(x);
}
