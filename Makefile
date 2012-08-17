CFLAGS += -Wall -pedantic -std=c99
LDFLAGS += -lm -lpthread
yee_SRC = yee.c
yee_OBJ = $(yee_SRC:.c=.o)
yee_pthr_SRC = yee_pthr.c
yee_pthr_OBJ = $(yee_pthr_SRC:.c=.o)
BIN = yee yee_pthr
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

%.o:%.c
	$(CC) $(CFLAGS) -c $<

$(DATA): yee
	./yee

$(GNUPLOT): gnuplot.plt
	sed 's/file/output_p/' $< > output_p.plt
	sed 's/file/output_u/' $< > output_u.plt

$(PNG): $(GNUPLOT) $(DATA) 
	gnuplot $(GNUPLOT)

clean:
	$(RM) $(BIN) $(yee_OBJ) $(yee_pthr_OBJ) $(DATA) $(PNG) $(GNUPLOT)
