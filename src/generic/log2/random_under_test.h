#include <math.h>
static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [exp(-1),exp(1)] */
  const double low = exp(-1), high = exp(1);
  return (high - low) * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) + low;
}
