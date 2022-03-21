static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [-4.16,4.16] */
  return 8.32 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - 4.16;
}
