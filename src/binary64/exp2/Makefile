FUNCTION_UNDER_TEST := exp2

include ../support/Makefile.univariate

all:: check_special

check_special: check_special.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) -fopenmp -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) -fopenmp -o $@ $^ -c

clean::
	rm -f check_special
