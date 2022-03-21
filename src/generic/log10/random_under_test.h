static inline TYPE_UNDER_TEST random_under_test (void)
{
  /* sample in [exp(-1),exp(1)] */
  const double low = 0x1.78b56362cef38p-2, high = 0x1.5bf0a8b145769p+1;
  return (high - low) * ((TYPE_UNDER_TEST) rand() / (TYPE_UNDER_TEST) RAND_MAX) + low;
}
