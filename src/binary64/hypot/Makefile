FUNCTION_UNDER_TEST := hypot

include ../support/Makefile.bivariate

all:: check_special

check_special: check_special.o triples.o $(CORE_MATH_OBJS) $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(CFLAGS) $(LDFLAGS) -fopenmp -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) $(CFLAGS) $(CORE_MATH_DEFINES) -I../support -c -fopenmp -o $@ $<

triples.o: triples.c
	$(CC) $(CFLAGS) $(CORE_MATH_DEFINES) -I../support -c -fopenmp -o $@ $<

clean::
	rm -f check_special
