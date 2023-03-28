FUNCTION_UNDER_TEST := cospi

include ../support/Makefile.univariate

all:: check_special

check_special: check_special.o $(CORE_MATH_OBJS) $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) $(CORE_MATH_DEFINES) -fopenmp -o $@ $^ -c

clean::
	rm -f check_special
