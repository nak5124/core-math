#define cr_function_under_test cr_acosf

static inline float random_under_test (void)
{
  /* acos is defined over [-1,1] */
  return 2 * ((float) rand() / (float) RAND_MAX) - 1;
}

static const double pi = 0x1.921fb54442d18p+1;

static inline float randomize_under_test (float x)
{
  return (x / pi) * random_under_test();
}
