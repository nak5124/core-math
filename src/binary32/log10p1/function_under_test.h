#define cr_function_under_test cr_log10p1f
#define ref_function_under_test ref_log10p1

void doit (uint32_t n);
static inline uint32_t asuint (float);
static inline float asfloat (uint32_t);

static inline int doloop (void)
{
  /* log1p is defined for x > -1 */
  uint32_t nmin = asuint (0x1p-149), nmax = asuint (0x1.fffffep127f);
#pragma omp parallel for schedule(dynamic,1024)
  for (uint32_t n = nmin; n <= nmax; n++)
  {
    doit (n);
    if (asfloat (n | 0x80000000) > -1)
      doit (n | 0x80000000);
  }
  printf ("all ok\n");
  return 0;
}
