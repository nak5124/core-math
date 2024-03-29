FUNCTION_UNDER_TEST := atan2f

include ../support/Makefile.bivariate

all:: check_special

check_special.o: check_special.c
	$(CC) $(CFLAGS) $(CORE_MATH_DEFINES) -c -fopenmp -o $@ $<

check_special: check_special.o $(CORE_MATH_OBJS) $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -lmpfr -lm

clean::
	rm -f check_special
