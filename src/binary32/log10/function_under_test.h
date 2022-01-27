#define cr_function_under_test cr_log10f
#define ref_function_under_test ref_log10

void doit (uint32_t n);
static inline uint32_t asuint (float f);

static inline int doloop (void)
{
  /* log10 is defined for x > 0 */
  uint32_t nmin = asuint (0x1p-149f), nmax = asuint (0x1.fffffep127f);
#pragma omp parallel for schedule(dynamic,1024)
  for (uint32_t n = nmin; n <= nmax; n++)
    doit (n);
  printf ("all ok\n");
  return 0;
}
