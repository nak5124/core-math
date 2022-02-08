FUNCTION_UNDER_TEST := hypotf

include ../support/Makefile.bivariate

all: check_special

check_special: check_special.o triples.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) -o $@ $^ -lmpfr -lm

clean::
	rm -f check_special
