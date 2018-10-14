###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.8.0
# @author Asim YarKhan
# @date 2010-11-15
#
###

PLASMA_DIR := $(PWD)/..
##include ../Makefile.internal
include ./make.inc

TRACE_SHARED_LIB_DESTINATION_DIRECTORY ?= $(PLASMA_DIR)/lib

all: libquark.a trace_shared_libs
lib: libquark.a trace_shared_libs

# Decide whether to enable EZTrace tracing, and to build the shared tracing library.
# Note, for binaries EZTrace should be in the include and library paths ( -leztrace -lfxt )
# EZTrace tracing is enabled by setting QUARK_TRACE=1 in a makefile (or make.inc). For example,
QUARK_TRACE ?= 0

ifeq ($(QUARK_TRACE), 1)

EZT_DIR ?= /usr
CFLAGS += -DTRACE_QUARK -I${EZT_DIR}/include
# LDFLAGS += -L${EZT_DIR}/lib -leztrace -lfxt
QUARK_TRACE_SO_LIB=$(TRACE_SHARED_LIB_DESTINATION_DIRECTORY)/libeztrace-convert-quark.so
QUARK_TRACE_OBJS=eztrace_convert_quark.o
$(QUARK_TRACE_SO_LIB): eztrace_convert_quark.o
	$(CC) -shared $(CFLAGS) $(INC) $^ -o $@
eztrace_convert_quark.o: eztrace_convert_quark.c quark_trace.h
	$(CC) -fPIC $(CFLAGS) $(INC) -c $< -o $@

QUARK_INSTALL_SO_LIB=install_eztrace_quark

install_eztrace_quark: $(QUARK_TRACE_SO_LIB)
	cp $(QUARK_TRACE_SO_LIB) $(prefix)/lib
endif

trace_shared_libs: $(QUARK_TRACE_SO_LIB)
install_shared_libs: $(QUARK_INSTALL_SO_LIB)

icl_hash.o: icl_hash.c icl_hash.h
icl_list.o: icl_list.c icl_list.h
quark.o: quark.c quark.h bsd_queue.h bsd_tree.h quark_trace.h
quarkos.o: quarkos.c quarkos-hwloc.c

libquark.a: icl_hash.o icl_list.o quarkos.o quark.o $(QUARK_TRACE_OBJS)
	$(ARCH) $(ARCHFLAGS) $@ $^
	$(RANLIB) $@

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

install: libquark.a install_shared_libs
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/lib
	cp quark.h             $(prefix)/include
	cp quark_unpack_args.h $(prefix)/include
	cp icl_hash.h          $(prefix)/include
	cp icl_list.h          $(prefix)/include
	cp libquark.a          $(prefix)/lib

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.PHONY: clean
