FUNCTION_UNDER_TEST := cbrtl

include ../support/Makefile.univariate

all:: check_special check_exact

check_special: check_special.o $(CORE_MATH_OBJS) $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) $(CFLAGS) $(CORE_MATH_DEFINES) -fopenmp -o $@ $^ -c

clean::
	rm -f check_special check_exact

check_exact: check_exact.o $(CORE_MATH_OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

check_exact.o: check_exact.c
	$(CC) $(CFLAGS) $(CORE_MATH_DEFINES) -o $@ $^ -c
