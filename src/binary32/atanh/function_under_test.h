#define cr_function_under_test cr_atanhf
#define ref_function_under_test ref_atanh

void doit (uint32_t n);
static inline uint32_t asuint (float f);

static inline int doloop (void)
{
  /* atanh is defined over (-1,1) */
  uint32_t nmin = asuint (0x0p0), nmax = asuint (0x1p0);
#pragma omp parallel for schedule(dynamic,1024)
  for (uint32_t n = nmin; n < nmax; n++)
  {
    doit (n);
    doit (n | 0x80000000);
  }
  printf ("all ok\n");
  return 0;
}
