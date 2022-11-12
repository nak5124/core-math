#define cr_function_under_test cr_acospif
#define ref_function_under_test ref_acospi

void doit (uint32_t n);
static inline uint32_t asuint (float f);

static inline int doloop (void)
{
  uint32_t nmin = asuint (0x0p0f), nmax = asuint (0x1p0f);
#pragma omp parallel for
  for (uint32_t n = nmin; n <= nmax; n++)
  {
    doit (n);
    doit (n | 0x80000000);
  }
  printf ("all ok\n");
  return 0;
}
