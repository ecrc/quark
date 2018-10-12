###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# This is the makefile for the dynamic scheduler component used by
# PLASMA
#
# @version 2.8.0
# @author Asim YarKhan
# @date 2010-11-15
#
###

CC       = gcc
CFLAGS   = -std=c99 -O2 -g -Wall
FC       = gfortran
FFLAGS   = -O2 -g -Wall
LOADER     = $(FC)
LDFLAGS  = -O2 -g 
LIBBLAS  = -L/home/yarkhan/local/brutus/lib -lreflapack -lgoto -lgfortran
LDLIBS   = ${LIBBLAS} -lpthread -lm
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = /usr/bin/ranlib

# Multi machine make file
ifeq ($(shell uname -n),Julien-Langous-Mac-Pro.local)
CC       = gcc
CFLAGS   = -std=c99 -O2 -g -Wall
FC       = gfortran
FFLAGS   = -O2 -g -Wall
LOADER     = $(FC)
LDFLAGS  = -O2 -g 
LDLIBS   = -framework vecLib -lgfortran -lpthread -lm
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = /usr/bin/ranlib
endif

ifeq ($(shell uname -n),scar)
CC       = gcc
CFLAGS   = -std=c99 -O0 -g -Wall
FC       = gfortran
FFLAGS   = -O0 -g -Wall  -lf2c
LOADER     = $(FC)
LDFLAGS  = -O2 -g 
LIBBLAS  = -L/home/yarkhan/local/brutus/lib -llapack -latlas -lgfortran
LDLIBS   = ${LIBBLAS} -lpthread -lm
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = /usr/bin/ranlib
endif

ifeq ($(shell uname -n),brutus)
CC       = gcc
CFLAGS   = -std=c99 -O2 -g -Wall
FC       = gfortran
FFLAGS   = -O2 -g -Wall
LOADER     = $(FC)
LDFLAGS  = -O2 -g 
LIBBLAS  = -L/home/yarkhan/local/brutus/lib -lreflapack -lgoto -lgfortran
LDLIBS   = ${LIBBLAS} -lpthread -lm
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = /usr/bin/ranlib
endif

ifeq ($(shell uname -n),zoot)
MKLPATH  = /mnt/scratch/sw/intel/C-11.0.083/mkl/lib/em64t
CC       = icc
CFLAGS   = -O2 -std=c99 -I${MKLPATH}/../../include
FC       = ifort
FFLAGS   = -O2
LOADER     = $(FC)
LDFLAGS  = -O2 -nofor_main 
LIBBLAS  = -L${MKLPATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LDLIBS   = ${LIBBLAS} -lpthread -lm
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = /usr/bin/ranlib
endif

ifeq ($(shell uname -n),ig.eecs.utk.edu)
MKLPATH  = /mnt/scratch/sw/intel/C-11.0.083/mkl/lib/em64t
# CC       = icc
# CFLAGS   = -O2 -g -std=c99 -fast -I${MKLPATH}/../../include
# FC       = ifort
# FFLAGS   = -O2 -g
# LOADER     = $(FC)
# LDFLAGS  = -O2 -g -nofor_main
# LIBBLAS  = -L${MKLPATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
CC       = gcc
CFLAGS   = -O2 -std=c99 -g -I${MKLPATH}/../../include -Wall
FC       = gfortran
FFLAGS   = -O2 -g -Wall 
LOADER     = $(FC)
LDFLAGS  = -O2 -g
LIBBLAS  = -L${MKLPATH} -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group
LDLIBS   = ${LIBBLAS} -lpthread -lm
ARCH     = xiar
ARCHFLAGS= cr
RANLIB   = /usr/bin/ranlib
endif

include ../../make.inc

####################
EXEC = tileLU tileChInv tileChol

all: tileLU tileChInv tileChol 

clean:	FORCE
	rm -f *.o $(EXEC) 

cleanall: clean
	rm -f *.a $(EXEC) profile.0.* tautrace.*.*.*.trc tau.trc.slog2 tau.trc tau.edf events.0.edf

tileLU: core_dgetrf.o core_dtstrf.o core_dgessm.o core_dssssm.o tileLU.o ../libscheduler.a
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

tileLU-run: tileLU
	./tileLU 28 252 7 2
	 #./tileLU 28 252 10 8
	 #./tileLU 20 200 40 8
	 #./tileLU 28 252 15 16

tileLU-run-numa: tileLU
	numactl --interleave=all ./tileLU 28 252 50 48

# Inversion
tileChInv-orig: tileChInv.o dpot03.o ../libscheduler.a
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

tileChInv: tileChInv.o dpot03.o ../libscheduler.a
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

tileChInv-run: tileChInv
	./tileChInv 30 10 1

# Nothing test
test: test.o core_dgetrf.o core_dtstrf.o core_dgessm.o core_dssssm.o ../libscheduler.a
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

test-run: tileLU
	./tileLU 28 252 7 2

# Cholesky Factorization
tileChol: tile_cholesky_3.o ../libscheduler.a
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

tileChol-run: tileChol
	./tileChol 30 8 
	 # ./tileChol BB THREADS  (NB=200 by default)


# Packaging
tar:
	tar cvzf scheduler-`date +'%Y%m%d'`.tgz Makefile icl_hash.h icl_hash.c icl_list.h icl_list.c scheduler.h scheduler.c unpack_args.h event_trace.h plasmaos.c core_*.f tileLU.c dpot03.f tileChInv.c


runall:	tileLU-run tileChInv-run tileChol-run

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

../libscheduler.a: FORCE
	( cd .. && $(MAKE) lib )

FORCE:
