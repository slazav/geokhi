all: theor.pdf

%.pdf: %.ps
	ps2pdf $<
%.ps: %.dvi
	dvips $<
%.dvi: %.tex
	latex $<

theor.dvi: theor.tex images/*.fig
	make -C images
	latex $<

clean:
	rm -f *.aux *.log *.bbl *.blg *.dvi *.ps *.pdf *.toc *.ind *.ilg
distclean: clean
	rm -f *.pdf

spell:
	ispell -d russian-lebedev ${parts_diss:=.tex} ${parts_auto:=.tex} svet.tex
