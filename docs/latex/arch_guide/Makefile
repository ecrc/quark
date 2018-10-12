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

PDFLATEX = pdflatex -file-line-error --shell-escape
SHELL = /bin/sh
CC = gcc
LATEX = latex
PDFLATEX = pdflatex -file-line-error --shell-escape
VIEW = open
BIBTEX = bibtex

BASENAME = quark_arch_guide

.PHONY: $(BASENAME).pdf

all: $(BASENAME).pdf

view: $(BASENAME).pdf
	open $(BASENAME).pdf

$(BASENAME).pdf:
	$(PDFLATEX) $(BASENAME).tex
	cp $(BASENAME).pdf ../../pdf/$(BASENAME).pdf 

clean:
	rm -f *.aux *.log *.out *.toc *.bbl *.blg $(BASENAME).pdf
	rm -f ../pdf/$(BASENAME).pdf

