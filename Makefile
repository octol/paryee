CC = gcc 
CFLAGS = -Wall -pedantic -Wextra -std=gnu99 -Ofast
MPICC = mpicc 
MPICFLAGS = -Wall -pedantic -Wextra -std=c99 -Ofast
LDFLAGS = -lm
OPENMP = -fopenmp

include Makefile.common

