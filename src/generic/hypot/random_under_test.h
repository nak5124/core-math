#define random_under_test_0 random_under_test
#define random_under_test_1 random_under_test

static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-10,10] */
  return 20 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - 10;
}
