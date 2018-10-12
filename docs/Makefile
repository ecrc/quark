###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.8.0
# @author Jakub Kurzak
# @date 2010-11-15
#
###

.PHONY: asciidoc latex doxygen

all: asciidoc latex doxygen

clean: cleanall

asciidoc:
	( cd asciidoc && $(MAKE) )

latex:
	( cd latex && $(MAKE) )

doxygen:
	( cd doxygen && $(MAKE) )

cleanall:
	( cd asciidoc && $(MAKE) clean )
	( cd latex    && $(MAKE) clean )
	( cd doxygen  && $(MAKE) clean )
