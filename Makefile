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
SHELL := /bin/bash
HOSTNAME = $(shell hostname)
ARCHNAME = $(shell uname -s)_$(shell uname -m)
TESTDIR = test
TESTDIR_PERF = perf-test
TESTDIR_PERF_OUT = perf-test-out/$(HOSTNAME)
SCRIPTDIR = helper-scripts
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

# Additional output for the reference version
#DATA += $(OUTDIR)/yee0.tsv
#$(OUTDIR)/yee0.tsv: $(BINDIR)/yee
#        $< -n $(gridlength) -o $@

#
# POSIX threads (blocked)
#
yee_block_pthr_SRC = $(SRCDIR)/yee_block_pthr.c \
	       $(SRCDIR)/yee_common.c 
EXTRA_LDFLAGS = -pthread
$(eval $(call DEF_BIN, $(BINDIR)/yee_block_pthr, $(yee_block_pthr_SRC), $(EXTRA_LDFLAGS), $(OUTDIR)/yee_block_pthr.tsv)) 

#
# POSIX threads (stride: 1)
#
yee_stride1_pthr_SRC = $(SRCDIR)/yee_stride1_pthr.c \
	       $(SRCDIR)/yee_common.c 
EXTRA_LDFLAGS = -pthread
$(eval $(call DEF_BIN, $(BINDIR)/yee_stride1_pthr, $(yee_stride1_pthr_SRC), $(EXTRA_LDFLAGS), $(OUTDIR)/yee_stride1_pthr.tsv)) 

#
# OpenMP (naive implementation, slow)
#
yee_naive_omp_SRC = $(SRCDIR)/yee_naive_omp.c \
	      $(SRCDIR)/yee_common.c 
EXTRA_LDFLAGS = $(OPENMP_FLAG)
$(eval $(call DEF_BIN, $(BINDIR)/yee_naive_omp, $(yee_naive_omp_SRC), $(EXTRA_LDFLAGS), $(OUTDIR)/yee_naive_omp.tsv)) 

#
# OpenMP (blocked version)
#
yee_block_omp_SRC = $(SRCDIR)/yee_block_omp.c \
	      $(SRCDIR)/yee_common.c 
EXTRA_LDFLAGS = $(OPENMP_FLAG)
$(eval $(call DEF_BIN, $(BINDIR)/yee_block_omp, $(yee_block_omp_SRC), $(EXTRA_LDFLAGS), $(OUTDIR)/yee_block_omp.tsv)) 

#
# OpenMP (stride: 1)
#
yee_stride1_omp_SRC = $(SRCDIR)/yee_stride1_omp.c \
	      $(SRCDIR)/yee_common.c 
EXTRA_LDFLAGS = $(OPENMP_FLAG)
$(eval $(call DEF_BIN, $(BINDIR)/yee_stride1_omp, $(yee_stride1_omp_SRC), $(EXTRA_LDFLAGS), $(OUTDIR)/yee_stride1_omp.tsv)) 

#
# MPI (blocking send/recieve)
#
yee_blocking_mpi_SRC = $(SRCDIR)/yee_blocking_mpi.c \
	      $(SRCDIR)/yee_common.c \
	      $(SRCDIR)/yee_common_mpi.c
$(eval $(call DEF_BIN_MPI, $(BINDIR)/yee_blocking_mpi, $(yee_blocking_mpi_SRC),, $(OUTDIR)/yee_blocking_mpi.tsv)) 

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

test-integration: all
	bindir=$(BINDIR) $(TESTDIR)/integration-tests.sh

data: all $(OUTDIR) $(DATA) 

plot: data $(PNG)

indent:
	indent -kr --no-tabs $(SRCDIR)/*.c
	indent -kr --no-tabs $(SRCDIR)/*.h

tags:
	ctags -R src

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
# Performance testing
#

test-scaling-compute: all
	./$(TESTDIR_PERF)/test_scaling.sh -n 8 -N 256 -s 4

#test-scaling: test-scaling-compute
test-scaling: 
	python3 $(SCRIPTDIR)/gather_data.py $(TESTDIR_PERF_OUT)/test_scaling_*.tsv \
	    > $(TESTDIR_PERF_OUT)/test_scaling.tsv
	cp $(TESTDIR_PERF)/test_scaling.plt $(TESTDIR_PERF_OUT)
	pushd $(TESTDIR_PERF_OUT) && \
	    gnuplot test_scaling.plt 

test-perf-compute: all
	#./$(TESTDIR_PERF)/test_perf.sh -n 4 -N 1000 -s 2 -t 50
	./$(TESTDIR_PERF)/test_perf.sh -n 4 -N 1000 -s 4 -t 1

#test-perf: test-perf-compute
test-perf: 
	python3 $(SCRIPTDIR)/gather_data.py $(TESTDIR_PERF_OUT)/test_perf_4_*.tsv \
	    > $(TESTDIR_PERF_OUT)/test_perf_4.tsv
	cp $(TESTDIR_PERF)/test_perf_4.plt $(TESTDIR_PERF_OUT)
	cp $(TESTDIR_PERF)/test_perf_8.plt $(TESTDIR_PERF_OUT)
	pushd $(TESTDIR_PERF_OUT) && \
	    gnuplot test_perf_4.plt 

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
