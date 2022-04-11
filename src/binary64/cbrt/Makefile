FUNCTION_UNDER_TEST := cbrt

include ../support/Makefile.univariate

all:: check_special

check_special: check_special.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) -fopenmp -o $@ $^ -lmpfr -lm

clean::
	rm -f check_special
