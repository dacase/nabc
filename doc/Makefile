#  set the lyx location here, if it is not in your path:
#LYX=/Applications/LyX.app/Contents/MacOS/lyx
LYX=lyx

SHELL = /bin/sh
.SUFFIXES: .txt .ps .html .lyx .pdf .tex .pic

.lyx.pdf:  $< amber.bib
	$(LYX) --export pdf2 $<
		
.lyx.tex:  $< amber.bib
	$(LYX) --export pdflatex $<
		
.lyx.txt:  $< amber.bib
	$(LYX) --export text $<
		
.lyx.ps:  $< amber.bib
	$(LYX) --export ps $<

.lyx.html: $< amber.bib
	elyxer.py --css=amber.css --title "$@" $<  $@
		
.pic.pdf: $<
	groff -p $< | ps2pdf - $@

Amber22.pdf: Amber.lyx amber.bib
	-$(LYX) --export pdf2 Amber.lyx
	/bin/mv Amber.pdf Amber22.pdf

clean::
	/bin/rm -f Amber21.pdf
