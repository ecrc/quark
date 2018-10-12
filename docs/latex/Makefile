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

.PHONY: users_guide arch_guide

all: users_guide arch_guide

users_guide:
	( cd users_guide && $(MAKE) )
arch_guide:
	( cd arch_guide && $(MAKE) )


clean:
	( cd users_guide && $(MAKE) clean)
	( cd arch_guide && $(MAKE) clean)
	rm -f ../pdf/quark_users_guide.pdf
	rm -f ../pdf/quark_arch_guide.pdf
