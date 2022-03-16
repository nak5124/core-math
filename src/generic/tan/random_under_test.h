#define pi 0x1.921fb54442d18p+1

static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-pi,pi] */
  return 2 * pi * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - pi;
}
