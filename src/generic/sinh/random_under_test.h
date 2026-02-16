#ifndef XMIN
#define XMIN -10
#endif

#ifndef XMAX
#define XMAX 10
#endif

static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-10,10] by default */
  // return 20 * ((double) rand() / (double) RAND_MAX) - 10;
  return XMIN + ((double) rand() / (double) RAND_MAX) * (XMAX - XMIN);
}
