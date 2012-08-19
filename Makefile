CC = gcc
CFLAGS += -Wall -pedantic -O3 
#CFLAGS += -march=native -mfpmath=sse -ffast-math
#CFLAGS += -g
LDFLAGS += -lm -lpthread
yee_SRC = yee.c yee_common.c
yee_pthr_SRC = yee_pthr.c yee_common.c
yee_pthr_barrier_SRC = yee_pthr_barrier.c yee_common.c
yee_ref_SRC = yee_ref.c
yee_OBJ = $(yee_SRC:.c=.o)
yee_pthr_OBJ = $(yee_pthr_SRC:.c=.o)
yee_pthr_barrier_OBJ = $(yee_pthr_barrier_SRC:.c=.o)
yee_ref_OBJ = $(yee_ref_SRC:.c=.o)
BIN = yee yee_pthr yee_pthr_barrier yee_ref
DATA = output_p.tsv output_u.tsv
GNUPLOT = $(DATA:.tsv=.plt)
PNG = $(DATA:.tsv=.png)

.PHONY: clean

all: $(BIN)

plot: $(PNG)

yee: $(yee_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@ 

yee_pthr: $(yee_pthr_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@ 

yee_pthr_barrier: $(yee_pthr_barrier_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@ 

yee_ref: $(yee_ref_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@ 

%.o:%.c
	$(CC) $(CFLAGS) -c $<

$(DATA): yee_ref
	./yee_ref

$(GNUPLOT): gnuplot.plt
	sed 's/file/output_p/' $< > output_p.plt
	sed 's/file/output_u/' $< > output_u.plt

$(PNG): $(GNUPLOT) $(DATA) 
	gnuplot $(GNUPLOT)

clean:
	$(RM) $(BIN) $(yee_OBJ) $(yee_pthr_OBJ) $(yee_pthr_barrier_OBJ) $(yee_ref_OBJ) $(DATA) $(PNG) $(GNUPLOT)
