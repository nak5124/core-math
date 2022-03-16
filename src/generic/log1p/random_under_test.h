#include <math.h>
static inline TYPE_UNDER_TEST random_under_test (void)
{
  const double low = exp(-0.5) - 1, high = exp(0.5) - 1;
  return (high - low) * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) + low;
}
