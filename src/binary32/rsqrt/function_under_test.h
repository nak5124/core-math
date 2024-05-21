#define cr_function_under_test cr_rsqrtf
#define ref_function_under_test ref_rsqrt

void doit (uint32_t n);
static inline uint32_t asuint (float f);
static inline float asfloat (uint32_t f);

static inline int doloop (void)
{
  uint32_t nmin = asuint (0x0p0f), nmax = asuint (0x1.fffffep127f);
#pragma omp parallel for schedule(dynamic,1024)
  for (uint32_t n = nmin; n <= nmax; n++)
  {
    doit (n);
    // also check -x
    float x = asfloat (n);
    doit (asuint (-x));
  }
  printf ("all ok\n");
  return 0;
}
