FUNCTION_UNDER_TEST := pow

include ../support/Makefile.bivariate

random: random.c
	gcc -W -Wall random.c -O3 -march=native -o r.out -lmpfr -lm
	@./r.out

test: test.c
	gcc -W -Wall -lm -lmpfr test.c -march=native -O3 -fno-builtin -o t.out
	@./t.out

check_exact: check_exact.c
	gcc -W -Wall -Wextra $(OPENMP) -O3 -march=native pow.c check_exact.c -lm -lmpfr -o ce.out

check: check_pow.c
	gcc -W -Wall -lm -lmpfr check_pow.c -O3 -march=native -fno-builtin -o c.out

exact: exact.c
	gcc -W -Wall -O3 -march=native exact.c -lm -lmpfr -o e.out

check_special: check_special.o $(FUNCTION_UNDER_TEST).o $(FUNCTION_UNDER_TEST)_mpfr.o
	$(CC) $(OPENMP) -o $@ $^ -lmpfr -lm

check_special.o: check_special.c
	$(CC) $(OPENMP) -o $@ $^ -c

clean_all: clean
	rm -f *.out
