# Simulation parameters
N = 1024
threads = 4

# Build parameters
CC = gcc 
CFLAGS += -Wall -pedantic -O3
#CFLAGS += -march=native -mfpmath=sse -ffast-math
MPICC = mpicc 
MPICFLAGS += -Wall -pedantic -O3 -std=c99
LDFLAGS += -lm

HOSTNAME = $(shell hostname)
SAVEDIR = tests_$(HOSTNAME)

# File dependencies
yee_SRC = yee.c yee_common.c
yee_pthr_SRC = yee_pthr.c yee_common.c
yee_omp_SRC = yee_omp.c yee_common.c
yee_mpi_SRC = yee_mpi.c yee_common.c
yee_ref_SRC = yee_ref.c

yee_OBJ = $(yee_SRC:.c=.o)
yee_pthr_OBJ = $(yee_pthr_SRC:.c=.o)
yee_omp_OBJ = $(yee_omp_SRC:.c=.o)
yee_mpi_OBJ = $(yee_mpi_SRC:.c=.o)
yee_ref_OBJ = $(yee_ref_SRC:.c=.o)
OBJ = $(yee_OBJ) $(yee_pthr_OBJ) $(yee_omp_OBJ) $(yee_mpi_OBJ) $(yee_ref_OBJ)

yee_BIN = yee
yee_pthr_BIN = yee_pthr
yee_omp_BIN = yee_omp
yee_ref_BIN = yee_ref
yee_mpi_BIN = yee_mpi
BIN = $(yee_BIN) $(yee_pthr_BIN) $(yee_omp_BIN) $(yee_mpi_BIN) $(yee_ref_BIN)

yee_DATA = $(yee_BIN:=_p.tsv) $(yee_BIN:=_u.tsv)
yee_pthr_DATA = $(yee_pthr_BIN:=_p.tsv) $(yee_pthr_BIN:=_u.tsv)
yee_omp_DATA = $(yee_omp_BIN:=_p.tsv) $(yee_omp_BIN:=_u.tsv)
yee_mpi_DATA = $(yee_mpi_BIN:=_p.tsv) $(yee_mpi_BIN:=_u.tsv)
yee_ref_DATA = $(yee_ref_BIN:=_p.tsv) $(yee_ref_BIN:=_u.tsv)
DATA = $(yee_DATA) $(yee_pthr_DATA) $(yee_omp_DATA) $(yee_mpi_DATA) $(yee_ref_DATA)

GNUPLOT = $(DATA:.tsv=.plt)
PNG = $(GNUPLOT:.plt=.png)

tests_SH = tests_scaling.sh tests_perf.sh
tests_DATA = $(tests_SH:.sh=.tsv)
tests_PNG = $(tests_DATA:.tsv=.png)

.PHONY: clean
.PRECIOUS: $(tests_DATA)

all: $(BIN) 

data: $(DATA)

plot: $(PNG)

tests: $(tests_PNG)

savetests: $(tests_PNG)
	[ ! -f $(SAVEDIR) ] && mkdir $(SAVEDIR) && echo "Create $(SAVEDIR)" || echo "Saving to $(SAVEDIR)"
	[ -d $(SAVEDIR) ] && cp $(tests_PNG) $(tests_DATA) $(SAVEDIR)

verify: $(BIN)
	./tests_integrity.sh

$(yee_BIN): $(yee_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@ 

$(yee_pthr_BIN): $(yee_pthr_OBJ)
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

$(yee_omp_DATA): $(yee_omp_BIN)
	./yee_omp -n $N -t $(threads)

$(yee_mpi_DATA): $(yee_mpi_BIN)
	mpirun -n $(shell dc -e '$(threads) 1 + p') yee_mpi -n $N

$(GNUPLOT): template.gnuplot
	sed 's/file/$(@:.plt=)/' $< > $@
	
%.png: %.plt %.tsv
	gnuplot $<

%.tsv: %.sh $(BIN)
	./$<

clean:
	$(RM) $(BIN) $(OBJ) $(DATA) $(PNG) $(GNUPLOT) $(tests_DATA) $(tests_PNG)
