FUNCTION_UNDER_TEST := cosh

include ../support/Makefile.univariate

all:: check_special

check_special: check_special.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -c

clean::
	rm -f check_special
