CC = gcc 
CFLAGS = -Wall -pedantic -Wextra -std=gnu99 -Ofast
#CFLAGS = -Wall -pedantic -Wextra -std=gnu99 -g -DDEBUG
MPICC = mpicc 
MPICFLAGS = -Wall -pedantic -Wextra -std=c99 -Ofast
#MPICFLAGS = -Wall -pedantic -Wextra -std=c99 -g -DDEBUG
LDFLAGS = -lm
#LDFLAGS = -lm -g
OPENMP = -fopenmp

include Makefile.common

