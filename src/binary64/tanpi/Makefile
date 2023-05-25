FUNCTION_UNDER_TEST := tanpi

include ../support/Makefile.univariate

check_special: check_special.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -c
