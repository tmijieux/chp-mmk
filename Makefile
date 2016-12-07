TARGET=projCHP test_driver
CFLAGS=-std=gnu99 -g -Wall -Wextra -fdiagnostics-color=auto
LDFLAGS=-lm
GENGETOPT=gengetopt
CC=mpicc

ifdef DEBUG
CFLAGS+=-ggdb -O0 -DDEBUG=1 -fsanitize=address -fsanitize=undefined
LDFLAGS+=-g-fsanitize=address -fsanitize=undefined
else
CFLAGS+=-O3 -march=native
endif

ifeq ($(strip $(BLASLIB)),)
else
LDFLAGS+= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \
	-Wl,--end-group ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a \
	-ldl -lpthread -lm -fopenmp
CFLAGS+=-DMKL
endif

SRC=    cmdline.c \
	main.c \
	error.c \
	perf/perf.c \
	util.c \
	grad.c \
	func.c \
	equation.c 

SRC_test_driver = \
	test_driver.c \
	util.c \
	grad.c \
	func.c

OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)

OBJ_test_driver=$(SRC_test_driver:.c=.o)

%.o: %.mod


all: $(TARGET)

-include $(DEP)

projCHP: $(OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

test_driver: $(OBJ_test_driver)
	$(CC) $^ -o $@ $(LDFLAGS) -lgfortran

%.o: %.c
	@$(CC) -MM $(CFLAGS) $*.c > $*.d
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	$(RM) $(OBJ_test_driver) $(OBJ) $(DEP) *.d *.o

mrproper: clean
	$(RM) $(TARGET)

genopt: projCHP.ggo
	$(GENGETOPT) < $^

