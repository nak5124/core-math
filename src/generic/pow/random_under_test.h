static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [0,10] since x^y gives NaN for x negative and y non-integer */
  return 10 * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX);
}
