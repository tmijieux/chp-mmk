CC=gfortran
CFLAGS=-O3 -march=native -fimplicit-none -Wall -Wextra \
	-Wline-truncation -fwhole-file
TARGETS=run test test_grad


SRC_run=fonct.f90 grad.f90 main.f90
OBJ_run=$(SRC_run:.f90=.o)

SRC_test=fonct.f90 grad.f90 test.f90
OBJ_test=$(SRC_test:.f90=.o)

SRC_test_grad=fonct.f90 grad.f90 test_grad.f90
OBJ_test_grad=$(SRC_test_grad:.f90=.o)


all: $(TARGETS)

%.o: %.f90
	$(CC) $(CFLAGS) $^ -c -o $@

debug: CFLAGS=-g -ggdb -O0 -Wall -Wextra -fcheck=all -fimplicit-none \
	-fbounds-check -fbacktrace -Wline-truncation -Wcharacter-truncation \
	-Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter \
	-fwhole-file
debug: run

run: $(OBJ_run)
	$(CC) $(CFLAGS) $^ -o $@

test: $(OBJ_test)
	$(CC) $(CFLAGS) $^ -o $@

test_grad: $(OBJ_test_grad)
	$(CC) $(CFLAGS) $^ -o $@

clear:
	$(RM) -r sol
	mkdir sol

clean:
	$(RM) *.o run *.mod
