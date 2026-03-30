/* acos is defined over [-1,1] */

#ifndef XMIN
#define XMIN -1
#endif

#ifndef XMAX
#define XMAX 1
#endif

static inline TYPE_UNDER_TEST random_under_test (void)
{
  return XMIN + (XMAX - XMIN) * ((double) rand() / (double) RAND_MAX);
}
