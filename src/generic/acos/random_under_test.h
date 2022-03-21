static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* acos is defined over [-1,1] */
  return 2 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) - 1;
}
