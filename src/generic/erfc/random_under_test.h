static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-5,5] */
  return 10.0 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - 5.0;
}
