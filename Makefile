ARCH ?= gcc

# -----------------------------------------------------------------------------
#  Environment
# -----------------------------------------------------------------------------

#
# Default GCC environment
# 
ifeq ($(ARCH), gcc)
CC = gcc 
MPICC = mpicc 
CFLAGS = -Wall -Wextra -Wpedantic -std=gnu99 
LDFLAGS = -lm
MPICFLAGS = $(CFLAGS)
MPILDFLAGS = $(LDFLAGS)
OPENMP_FLAG = -fopenmp

ifneq (,$(DEBUG))
CFLAGS += -g -DDEBUG
MPICFLAGS += -g -DDEBUG
LDFLAGS += -g
MPILDFLAGS += -g
else
CFLAGS += -Ofast
MPICFLAGS += -Ofast
endif
endif

#
# Solaris environment
#
ifeq ($(ARCH), solaris)
CC = cc
MPICC = mpicc 
CFLAGS = -m64 -fast
MPICFLAGS = -m64 -fast
LDFLAGS = -m64 -lm -fast
MPILDFLAGS = -m64 -fast
OPENMP = -xopenmp
endif

#
# SGI environment
#
ifeq ($(ARCH), irix)
CC = c99
MPICC = c99
CFLAGS = -64 -Ofast
MPICFLAGS = -64 -Ofast
LDFLAGS = -64 -lm -Ofast
MPILDFLAGS = -64 -lm -Ofast -lmpi
OPENMP = -mp
endif

#
# IA-64 (Itanium) environment
#
ifeq ($(ARCH), ia64)
CC = opencc
MPICC = mpicc
CFLAGS = -std=gnu99 -Ofast
MPICFLAGS = -std=gnu99 -Ofast
LDFLAGS = -lm -Ofast
MPILDFLAGS = -lm -Ofast
OPENMP = -mp
endif

# -----------------------------------------------------------------------------
# Build config
# -----------------------------------------------------------------------------

#
# Simulation parameters
#
N=128
threads=4

#
# Environment
#
HOSTNAME = $(shell hostname)
ARCHNAME = $(shell uname -s)_$(shell uname -m)
SAVEDIR = tests_$(HOSTNAME)
TESTDIR = test
SRCDIR = src
ifneq (,$(DESTDIR))
OBJDIR = $(DESTDIR)/$(ARCHNAME)/obj
BINDIR = $(DESTDIR)/$(ARCHNAME)/bin
OUTDIR = $(DESTDIR)/$(ARCHNAME)/output
else
OBJDIR = $(ARCHNAME)/obj
BINDIR = $(ARCHNAME)/bin
OUTDIR = $(ARCHNAME)/output
endif

OBJ =
BIN =
DATA =

# -----------------------------------------------------------------------------
#  Macros
# -----------------------------------------------------------------------------

# Macro to create linking target
# $(1) Target bin, $(2) Sources
define LINKER_TARGET
OBJ += $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(2))
BIN += $(1)

$(1): $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(2))
	$(CC) $$^ $(LDFLAGS) -o $$@ 

$(1)_DATA = $(patsubst $(BINDIR)/%,$(OUTDIR)/%.tsv,$(BINDIR)/$(1))
$(1)_DATA: $(1)
	$$< -n $N -o $$@
endef

# -----------------------------------------------------------------------------
#  Targets
# -----------------------------------------------------------------------------

#
# Meta
#

all: $(OBJDIR) $(BINDIR) bin

#
# Serial Yee
#
yee_SRC = $(SRCDIR)/yee.c \
	  $(SRCDIR)/yee_common.c
#yee_OBJ = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(yee_SRC))
#yee_BIN = $(BINDIR)/yee
#yee_DATA = $(patsubst $(BINDIR)/%,$(OUTDIR)/%.tsv,$(yee_BIN))
#OBJ += $(yee_OBJ)
#BIN += $(yee_BIN)
#DATA += $(yee_DATA)

$(eval $(call LINKER_TARGET, $(BINDIR)/yee, $(yee_SRC))) 

#
# Pthread
#
yee_pthr_SRC = $(SRCDIR)/yee_pthr.c \
	       $(SRCDIR)/yee_common.c 
yee_pthr_OBJ = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(yee_pthr_SRC))
yee_pthr_BIN = $(BINDIR)/yee_pthr
yee_pthr_DATA = $(patsubst $(BINDIR)/%,$(OUTDIR)/%.tsv,$(yee_pthr_BIN))
OBJ += $(yee_pthr_OBJ)
BIN += $(yee_pthr_BIN)
DATA += $(yee_pthr_DATA)

#
# OpenMP
#
yee_omp_SRC = $(SRCDIR)/yee_omp.c \
	      $(SRCDIR)/yee_common.c 
yee_omp_OBJ = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(yee_omp_SRC))
yee_omp_BIN = $(BINDIR)/yee_omp
yee_omp_DATA = $(patsubst $(BINDIR)/%,$(OUTDIR)/%.tsv,$(yee_omp_BIN))
OBJ += $(yee_omp_OBJ)
BIN += $(yee_omp_BIN)
DATA += $(yee_omp_DATA)

#
# MPI
#
yee_mpi_SRC = $(SRCDIR)/yee_mpi.c \
	      $(SRCDIR)/yee_common.c \
	      $(SRCDIR)/yee_common_mpi.c
yee_mpi_OBJ = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(yee_mpi_SRC))
yee_mpi_BIN = $(BINDIR)/yee_mpi
yee_mpi_DATA = $(patsubst $(BINDIR)/%,$(OUTDIR)/%.tsv,$(yee_mpi_BIN))
OBJ += $(yee_mpi_OBJ)
BIN += $(yee_mpi_BIN)
DATA += $(yee_mpi_DATA)

#
# MPI (Non-blocking)
#
yee_mpi2_SRC = $(SRCDIR)/yee_mpi2.c \
	       $(SRCDIR)/yee_common.c \
	       $(SRCDIR)/yee_common_mpi.c
yee_mpi2_OBJ = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(yee_mpi2_SRC))
yee_mpi2_BIN = $(BINDIR)/yee_mpi2
yee_mpi2_DATA = $(patsubst $(BINDIR)/%,$(OUTDIR)/%.tsv,$(yee_mpi2_BIN))
OBJ += $(yee_mpi2_OBJ)
BIN += $(yee_mpi2_BIN)
DATA += $(yee_mpi2_DATA)

#
# Tests
# 
unittests_SRC = $(SRCDIR)/yee_common_tests.c \
		$(SRCDIR)/yee_common.c
unittests_OBJ = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(unittests_SRC))
unittests_BIN = $(BINDIR)/unit_tests

#
# Output data
#
tests_DATA = $(OUTDIR)/tests_perf_4.tsv \
	     $(OUTDIR)/tests_perf_8.tsv
tests_PNG = $(tests_DATA:.tsv=.png)

GNUPLOT = $(DATA:.tsv=.plt) 
PNG = $(GNUPLOT:.plt=.png) 

.PHONY: clean
.PRECIOUS: $(tests_DATA)

# -----------------------------------------------------------------------------
#  Rules
# -----------------------------------------------------------------------------

#
# Meta targets 
#

bin: $(BIN)

data: all $(OUTDIR) $(DATA)

plot: data $(PNG)

integration-test: all
	$(TESTDIR)/integration-tests.sh

test: $(unittests_BIN)
	./$<

#
# General targets
#

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

#
# Specific targets
#

# Serial
#$(yee_BIN): $(yee_OBJ)
#	$(CC) $^ $(LDFLAGS) -o $@ 

#$(yee_DATA): $(yee_BIN)
#	$< -n $N -o $@

# Pthread
$(yee_pthr_BIN): $(yee_pthr_OBJ)
	$(CC) $^ $(LDFLAGS) -lpthread -o $@ 

$(yee_pthr_DATA): $(yee_pthr_BIN)
	$< -n $N -t $(threads) -o $@

# OpenMP
$(yee_omp_BIN): $(yee_omp_OBJ)
	$(CC) $^ $(LDFLAGS) $(OPENMP_FLAG) -o $@ 

$(OBJDIR)/yee_omp.o: $(SRCDIR)/yee_omp.c
	$(CC) $(CFLAGS) -c $< $(OPENMP_FLAG) -o $@

$(yee_omp_DATA): $(yee_omp_BIN)
	$< -n $N -t $(threads) -o $@

# MPI Common

$(OBJDIR)/yee_common_mpi.o: $(SRCDIR)/yee_common_mpi.c $(SRCDIR)/yee_common.c
	$(MPICC) $(MPICFLAGS) -c $< -o $@

# MPI
$(yee_mpi_BIN): $(yee_mpi_OBJ)
	$(MPICC) $^ $(MPILDFLAGS) -o $@

$(OBJDIR)/yee_mpi.o: $(SRCDIR)/yee_mpi.c
	$(MPICC) $(MPICFLAGS) -c $< -o $@

$(yee_mpi_DATA): $(yee_mpi_BIN)
	mpirun -n $(threads) $< -n $N -o $@

# MPI (Non-blocking)
$(yee_mpi2_BIN): $(yee_mpi2_OBJ)
	$(MPICC) $^ $(MPILDFLAGS) -o $@

$(OBJDIR)/yee_mpi2.o: $(SRCDIR)/yee_mpi2.c
	$(MPICC) $(MPICFLAGS) -c $< -o $@

$(yee_mpi2_DATA): $(yee_mpi_BIN)
	mpirun -n $(threads) $< -n $N -o $@

#
# Unit tests
#

$(unittests_BIN): $(unittests_OBJ)
	$(CC) $^ $(LDFLAGS) -lcunit -o $@ 

#
# Generate plots
#

$(GNUPLOT): $(SRCDIR)/template.gnuplot
	sed 's:file:$(@:.plt=):g' $< > $@

%.png: %.plt %.tsv
	if [ -f /usr/bin/gnuplot ]; then gnuplot $<; else touch $<; fi
 
%.tsv: %.sh $(BIN)
	./$<

$(BINDIR):
	mkdir -p $@

$(OBJDIR):
	mkdir -p $@

$(OUTDIR):
	mkdir -p $@

clean:
	$(RM) $(BIN) $(OBJ) $(DATA) $(PNG) $(GNUPLOT) $(unittests_BIN) $(unittests_OBJ)
	if [ -d $(OBJDIR) ]; then rmdir $(OBJDIR); fi
	if [ -d $(BINDIR) ]; then rmdir $(BINDIR); fi
	if [ -d $(OUTDIR) ]; then rmdir $(OUTDIR); fi
	if [ -d $(ARCHNAME) ]; then rmdir $(ARCHNAME); fi
