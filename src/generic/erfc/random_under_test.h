static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-4.07,10.68] */
  return 14.75 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - 4.07;
}
