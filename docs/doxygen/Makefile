DOXYGEN ?= doxygen
VERSION ?= 2.3.1


all: install

install: doxygen install-html install-latex

install-html: doxygen
	( test ! -d ../html/doxygen && mkdir ../html/doxygen ) || true
	cp -rf out/html/* ../html/doxygen/.

install-latex: doxygen
	( test ! -d ../latex/doxygen && mkdir ../latex/doxygen ) || true
	( cd out/latex && $(MAKE) && acroread refman.pdf)
	cp out/latex/*.pdf ../latex/doxygen/.

doxygen:
	if	test -d ../../../quark; \
		then make doxygen-dev; \
		else make doxygen-usr; \
	fi

doxygen-usr: quark-usr.dox groups-usr.dox
	$(DOXYGEN) quark-usr.dox

doxygen-dev: quark-dev.dox groups-usr.dox groups-dev.dox
	$(DOXYGEN) quark-dev.dox

quark-usr.dox: conf-usr.dox
	sed s/!VERSION!/$(VERSION)/g conf-usr.dox > quark-usr.dox

quark-dev.dox: conf-dev.dox
	sed s/!VERSION!/$(VERSION)/g conf-dev.dox > quark-dev.dox

clean:
	rm -rf ../html/doxygen
	rm -rf ../latex/doxygen
	rm -rf html
	rm -rf latex
	rm -rf out
	rm -rf quark-usr.dox
	rm -rf quark-dev.dox

.PHONY: clean doxygen doxygen-usr doxygen-dev install-latex install-html install all
