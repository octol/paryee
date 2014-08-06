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
gridlength=128
threads=4

#
# Environment
#
HOSTNAME = $(shell hostname)
ARCHNAME = $(shell uname -s)_$(shell uname -m)
TESTDIR = test
MPIRUN = mpirun
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
OBJ_TESTS =
BIN_TESTS =

# -----------------------------------------------------------------------------
#  Macros
# -----------------------------------------------------------------------------

# Macro to create linking target
# $(1) Target bin, $(2) Sources, $(3) Extra LDFLAGS, $(4) Output data
define DEF_BIN
OBJ += $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(2))
BIN += $(1)

$(1): $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(2))
	$(CC) $$^ $(LDFLAGS) $(3) -o $$@ 

DATA += $(4)
$(4): $(1)
	$$< -n $(gridlength) -o $$@
endef

# Macro to create linking target using MPI
# $(1) Target bin, $(2) Sources, $(3) Extra LDFLAGS, $(4) Output data
define DEF_BIN_MPI
OBJ += $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(2))
BIN += $(1)

$(1): $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(2))
	$(MPICC) $$^ $(MPILDFLAGS) $(3) -o $$@ 

DATA += $(4)
$(4): $(1)
	$(MPIRUN) -n $(threads) $$< -n $(gridlength) -o $$@
endef

# Macro to create linking target
# $(1) Target bin, $(2) Sources, $(3) Extra LDFLAGS
define DEF_TEST
OBJ_TESTS += $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(2))
BIN_TESTS += $(1)

$(1): $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(2))
	$(CC) $$^ $(LDFLAGS) -lcunit $(3) -o $$@ 
endef

# -----------------------------------------------------------------------------
#  Targets
# -----------------------------------------------------------------------------

#
# Meta
#
all: $(OBJDIR) $(BINDIR) bin

#
# Serial reference version
#
yee_SRC = $(SRCDIR)/yee.c \
	  $(SRCDIR)/yee_common.c
$(eval $(call DEF_BIN, $(BINDIR)/yee, $(yee_SRC),, $(OUTDIR)/yee.tsv)) 

#
# POSIX threads
#
yee_pthr_SRC = $(SRCDIR)/yee_pthr.c \
	       $(SRCDIR)/yee_common.c 
EXTRA_LDFLAGS = -pthread
$(eval $(call DEF_BIN, $(BINDIR)/yee_pthr, $(yee_pthr_SRC), $(EXTRA_LDFLAGS), $(OUTDIR)/yee_pthr.tsv)) 

#
# OpenMP
#
yee_omp_SRC = $(SRCDIR)/yee_omp.c \
	      $(SRCDIR)/yee_common.c 
EXTRA_LDFLAGS = $(OPENMP_FLAG)
$(eval $(call DEF_BIN, $(BINDIR)/yee_omp, $(yee_omp_SRC), $(EXTRA_LDFLAGS), $(OUTDIR)/yee_omp.tsv)) 

#
# MPI
#
yee_mpi_SRC = $(SRCDIR)/yee_mpi.c \
	      $(SRCDIR)/yee_common.c \
	      $(SRCDIR)/yee_common_mpi.c
$(eval $(call DEF_BIN_MPI, $(BINDIR)/yee_mpi, $(yee_mpi_SRC),, $(OUTDIR)/yee_mpi.tsv)) 

#
# MPI (Non-blocking)
#
yee_nonblock_mpi_SRC = $(SRCDIR)/yee_nonblock_mpi.c \
		       $(SRCDIR)/yee_common.c \
	               $(SRCDIR)/yee_common_mpi.c
$(eval $(call DEF_BIN_MPI, $(BINDIR)/yee_nonblock_mpi, $(yee_nonblock_mpi_SRC),, $(OUTDIR)/yee_nonblock_mpi.tsv)) 

#
# Tests
# 
unittests_SRC = $(SRCDIR)/yee_common_tests.c \
		$(SRCDIR)/yee_common.c
$(eval $(call DEF_TEST, $(BINDIR)/unit_tests, $(unittests_SRC)))

# #
# # Output data
# #
# tests_DATA = $(OUTDIR)/tests_perf_4.tsv \
# 	     $(OUTDIR)/tests_perf_8.tsv
# tests_PNG = $(tests_DATA:.tsv=.png)
 
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

test: $(BIN_TESTS)
	./$<

integration-test: all
	bindir=$(BINDIR) $(TESTDIR)/integration-tests.sh

data: all $(OUTDIR) $(DATA) 

plot: data $(PNG)

#
# Compilation
#

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%_omp.o: $(SRCDIR)/%_omp.c
	$(CC) $(CFLAGS) -c $< $(OPENMP_FLAG) -o $@

$(OBJDIR)/%_mpi.o: $(SRCDIR)/%_mpi.c
	$(MPICC) $(MPICFLAGS) -c $< -o $@

#
# Generate plots
#
 
$(GNUPLOT): $(SRCDIR)/template.gnuplot
	sed 's:file:$(@:.plt=):g' $< > $@

%.png: %.plt %.tsv
	gnuplot $<

#
# Environment
#

$(BINDIR) $(OBJDIR) $(OUTDIR):
	mkdir -p $@

clean:
	$(RM) $(BIN) $(BIN_TESTS) $(OBJ) $(OBJ_TESTS) \
	    $(DATA) $(PNG) $(GNUPLOT) $(unittests_BIN) $(unittests_OBJ)
	[ -d $(OBJDIR) ] && rmdir $(OBJDIR)
	[ -d $(BINDIR) ] && rmdir $(BINDIR)
	[ -d $(OUTDIR) ] && rmdir $(OUTDIR)
	[ -d $(ARCHNAME) ] && rmdir $(ARCHNAME)
