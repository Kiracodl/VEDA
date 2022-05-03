all: VEDA_manual.pdf

#bibtex requires multiple passes to resolve all references (kinda goofy but that's what you need to do!
FirstTimeUsers.pdf: FirstTimeUsersManual.tex Library.bib
	pdflatex FirstTimeUsersManual && bibtex FirstTimeUsersManual && pdflatex FirstTimeUsersManual && pdflatex FirstTimeUsersManual

VEDA_manual.pdf: VEDA_Manual.tex Library.bib
	pdflatex VEDA_Manual && bibtex VEDA_Manual && pdflatex VEDA_Manual && pdflatex VEDA_Manual

clean:
	rm -f VEDA_Manual.pdf *.toc *.aux *.out *.bbl *.end *.blg