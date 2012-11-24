# Simulation parameters
N = 128
threads = 8

# Build parameters
CC = gcc 
CFLAGS += -Wall -pedantic
#CFLAGS += -march=native -mfpmath=sse -ffast-math
CFLAGS += -g
MPICC = mpicc 
MPICFLAGS += -Wall -pedantic -O3
MPICFLAGS += -std=c99
LDFLAGS += -lm

# File dependencies
yee_SRC = yee.c yee_common.c
yee_pthr_SRC = yee_pthr.c yee_common.c
yee_pthr_barrier_SRC = yee_pthr_barrier.c yee_common.c
yee_omp_SRC = yee_omp.c yee_common.c
yee_mpi_SRC = yee_mpi.c yee_common.c
yee_ref_SRC = yee_ref.c

yee_OBJ = $(yee_SRC:.c=.o)
yee_pthr_OBJ = $(yee_pthr_SRC:.c=.o)
yee_pthr_barrier_OBJ = $(yee_pthr_barrier_SRC:.c=.o)
yee_omp_OBJ = $(yee_omp_SRC:.c=.o)
yee_mpi_OBJ = $(yee_mpi_SRC:.c=.o)
yee_ref_OBJ = $(yee_ref_SRC:.c=.o)
OBJ = $(yee_OBJ) $(yee_pthr_OBJ) $(yee_pthr_barrier_OBJ) $(yee_omp_OBJ) $(yee_mpi_OBJ) $(yee_ref_OBJ)

yee_BIN = yee
yee_pthr_BIN = yee_pthr
yee_pthr_barrier_BIN = yee_pthr_barrier
yee_omp_BIN = yee_omp
yee_ref_BIN = yee_ref
yee_mpi_BIN = yee_mpi
BIN = $(yee_BIN) $(yee_pthr_BIN) $(yee_pthr_barrier_BIN) $(yee_omp_BIN) $(yee_mpi_BIN) $(yee_ref_BIN)

yee_DATA = $(yee_BIN:=_p.tsv) $(yee_BIN:=_u.tsv)
yee_pthr_DATA = $(yee_pthr_BIN:=_p.tsv) $(yee_pthr_BIN:=_u.tsv)
yee_pthr_barrier_DATA = $(yee_pthr_barrier_BIN:=_p.tsv) $(yee_pthr_barrier_BIN:=_u.tsv)
yee_omp_DATA = $(yee_omp_BIN:=_p.tsv) $(yee_omp_BIN:=_u.tsv)
yee_mpi_DATA = $(yee_mpi_BIN:=_p.tsv) $(yee_mpi_BIN:=_u.tsv)
yee_ref_DATA = $(yee_ref_BIN:=_p.tsv) $(yee_ref_BIN:=_u.tsv)
DATA = $(yee_DATA) $(yee_pthr_DATA) $(yee_pthr_barrier_DATA) $(yee_omp_DATA) $(yee_mpi_DATA) $(yee_ref_DATA)

GNUPLOT = $(DATA:.tsv=.plt)
PNG = $(GNUPLOT:.plt=.png)

.PHONY: clean

all: $(BIN) 

data: $(DATA)

plot: $(PNG)

$(yee_BIN): $(yee_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@ 

$(yee_pthr_BIN): $(yee_pthr_OBJ)
	$(CC) $^ $(LDFLAGS) -lpthread -o $@ 

$(yee_pthr_barrier_BIN): $(yee_pthr_barrier_OBJ)
	$(CC) $^ $(LDFLAGS) -lpthread -o $@ 

$(yee_omp_BIN): $(yee_omp_OBJ)
	$(CC) $^ $(LDFLAGS) -fopenmp -o $@ 

$(yee_ref_BIN): $(yee_ref_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@ 

yee_mpi: $(yee_mpi_OBJ)
	$(MPICC) $^ -o $@

yee_omp.o: yee_omp.c
	$(CC) $(CFLAGS) -c $< -fopenmp

yee_mpi.o: yee_mpi.c
	$(MPICC) $(MPICFLAGS) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

$(yee_DATA): $(yee_BIN)
	./yee -n $N

$(yee_ref_DATA): $(yee_ref_BIN)
	./yee_ref -n $N

$(yee_pthr_DATA): $(yee_pthr_BIN)
	./yee_pthr -n $N -t $(threads)

$(yee_pthr_barrier_DATA): $(yee_pthr_barrier_BIN)
	./yee_pthr_barrier -n $N -t $(threads)

$(yee_omp_DATA): $(yee_omp_BIN)
	./yee_omp -n $N -t $(threads)

$(yee_mpi_DATA): $(yee_mpi_BIN)
	mpirun -n $(shell dc -e '$(threads) 1 + p') yee_mpi -n $N

$(GNUPLOT): template.gnuplot
	sed 's/file/$(@:.plt=)/' $< > $@
	
%.png: %.plt %.tsv
	gnuplot $<

clean:
	$(RM) $(BIN) $(OBJ) $(DATA) $(PNG) $(GNUPLOT)
