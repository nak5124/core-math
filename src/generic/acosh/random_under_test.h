static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [1,21] */
  return 20 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) + 1;
}
