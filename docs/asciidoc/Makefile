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

.PHONY: readme

all: readme

readme:
	asciidoc -n -a toc README.txt
	mv -f README.html ../html/README.html
	asciidoc -n -a toc -b html4 -o - README.txt | lynx -justify -dump -stdin > ../../README

clean:
	rm -f ../html/README.html
	rm -f ../../README
