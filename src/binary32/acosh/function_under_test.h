#define cr_function_under_test cr_acoshf
#define ref_function_under_test ref_acosh

void doit (uint32_t n);
static inline uint32_t asuint (float f);

static inline int doloop (void)
{
  /* acosh is only defined for x >= 1 */
  uint32_t nmin = asuint (0x1p0), nmax = asuint (0x1.fffffep127f);
#ifndef NO_OPENMP
#pragma omp parallel for
#endif
  for (uint32_t n = nmin; n <= nmax; n++)
    doit (n);
  printf ("all ok\n");
  return 0;
}
