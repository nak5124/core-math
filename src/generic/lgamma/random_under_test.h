static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-20,20] */
  return 40. * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - 20.0;
}
