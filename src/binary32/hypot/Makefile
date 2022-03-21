FUNCTION_UNDER_TEST := hypotf

include ../support/Makefile.bivariate

all:: check_special

check_special: check_special.o triples.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) -fopenmp -o $@ $^ -lmpfr -lm

triples.o: triples.c
	$(CC) $(CFLAGS) -c -fopenmp -o $@ $<

clean::
	rm -f check_special
