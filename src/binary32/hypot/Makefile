FUNCTION_UNDER_TEST := hypotf

include ../support/Makefile.bivariate

all:: check_special

check_special: check_special.o triples.o $(CORE_MATH_OBJS) $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(LDFLAGS) -fopenmp -o $@ $^ -lmpfr -lm

triples.o: triples.c
	$(CC) $(CFLAGS) $(CORE_MATH_DEFINES) -c -fopenmp -o $@ $<

clean::
	rm -f check_special
