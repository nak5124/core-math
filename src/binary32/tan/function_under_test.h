#define cr_function_under_test cr_tanf
#define ref_function_under_test ref_tan

void doit (uint32_t n);
static inline uint32_t asuint (float f);

static inline int doloop (void)
{
  uint32_t nmin = asuint (0x0p0f), nmax = asuint (0x1.fffffep+127);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint32_t n = nmin; n <= nmax; n++)
  {
    doit (n);
    doit (n | 0x80000000);
  }
  printf ("all ok\n");
  return 0;
}
