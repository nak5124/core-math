static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-5,5] */
#define XMIN -5.0
#define XMAX +5.0
  return (XMAX - XMIN) * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) + XMIN;
}
