#define cr_function_under_test cr_acosf
#define function_under_test acosf

static inline float random_under_test (void)
{
  /* acos is defined over [-1,1] */
  return 2 * ((float) rand() / (float) RAND_MAX) - 1;
}
