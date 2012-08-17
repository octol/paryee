CFLAGS += -Wall -pedantic -std=c99
LDFLAGS += -lm -lpthread
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
BIN = yee
DATA = output_p.tsv output_u.tsv
GNUPLOT = $(DATA:.tsv=.plt)
PNG = $(DATA:.tsv=.png)

.PHONY: clean

all: $(BIN)

plot: $(PNG)

yee: $(OBJ)
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
	$(RM) $(BIN) $(OBJ) $(DATA) $(PNG) $(GNUPLOT)
