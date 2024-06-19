#define cr_function_under_test cr_rsqrtf
#define ref_function_under_test ref_rsqrt

void doit (uint32_t n);
static inline uint32_t asuint (float f);
static inline float asfloat (uint32_t f);

static inline int doloop (void)
{
  uint32_t nmin = asuint (0x0p0f), nmax = asuint (0x1.fffffep127f);
#pragma omp parallel for
  for (uint32_t n = nmin; n <= nmax; n++)
  {
    doit (n);
    doit (n | 0x80000000); // also test negative numbers
  }
  printf ("all ok\n");
  return 0;
}
