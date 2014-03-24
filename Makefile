CC = gcc 
MPICC = mpicc 

CFLAGS_BASIC = -Wall -Wextra -Wpedantic -std=gnu99 
LDFLAGS_BASIC = -lm

CFLAGS = $(CFLAGS_BASIC) -Ofast
LDFLAGS = $(LDFLAGS_BASIC)
#CFLAGS = $(CFLAGS_BASIC) -g -DDEBUG
#LDFLAGS = $(LDFLAGS_BASIC) -g

MPICFLAGS = $(CFLAGS_BASIC) -Ofast
MPILDFLAGS = $(LDFLAGS_BASIC)
#MPICFLAGS = $(CFLAGS_BASIC) -g -DDEBUG
#MPILDFLAGS = $(LDFLAGS_BASIC) -g

OPENMP_FLAG = -fopenmp

include Makefile.common

