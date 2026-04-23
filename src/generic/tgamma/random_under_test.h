#ifndef XMAX
#define XMAX 20
#endif

#ifndef XMIN
#define XMIN -XMAX
#endif

static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-20,20] */
  return XMIN + (XMAX - XMIN) * ((double) rand() / (double) RAND_MAX);
}
