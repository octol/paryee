CC = gcc
CFLAGS += -Wall -pedantic 
#CFLAGS += -g
LDFLAGS += -lm -lpthread
yee_SRC = yee.c
yee_OBJ = $(yee_SRC:.c=.o)
yee_pthr_SRC = yee_pthr.c
yee_pthr_OBJ = $(yee_pthr_SRC:.c=.o)
yee_pthr2_SRC = yee_pthr2.c
yee_pthr2_OBJ = $(yee_pthr2_SRC:.c=.o)
BIN = yee yee_pthr yee_pthr2
DATA = output_p.tsv output_u.tsv
GNUPLOT = $(DATA:.tsv=.plt)
PNG = $(DATA:.tsv=.png)

.PHONY: clean

all: $(BIN)

plot: $(PNG)

yee: $(yee_OBJ)
	$(CC) $< $(LDFLAGS) -o $@ 

yee_pthr: $(yee_pthr_OBJ)
	$(CC) $< $(LDFLAGS) -o $@ 

yee_pthr2: $(yee_pthr2_OBJ)
	$(CC) $< $(LDFLAGS) -o $@ 

%.o:%.c
	$(CC) $(CFLAGS) -c $<

$(DATA): yee_pthr
	./yee_pthr

$(GNUPLOT): gnuplot.plt
	sed 's/file/output_p/' $< > output_p.plt
	sed 's/file/output_u/' $< > output_u.plt

$(PNG): $(GNUPLOT) $(DATA) 
	gnuplot $(GNUPLOT)

clean:
	$(RM) $(BIN) $(yee_OBJ) $(yee_pthr_OBJ) $(yee_pthr2_OBJ) $(DATA) $(PNG) $(GNUPLOT)
