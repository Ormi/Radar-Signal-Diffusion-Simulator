all: pdf

pdf: 2017-ExcelFIT-ormos.tex
	pdflatex 2017-ExcelFIT-ormos.tex
	bibtex 2017-ExcelFIT-ormos.aux
	pdflatex 2017-ExcelFIT-ormos.tex
	pdflatex 2017-ExcelFIT-ormos.tex
	#dvips 2017-ExcelFIT-ormos.dvi
	#ps2pdf 2017-ExcelFIT-ormos.ps

	rm 2017-ExcelFIT-ormos.aux
	#rm 2017-ExcelFIT-ormos.dvi
	#rm 2017-ExcelFIT-ormos.ps
	rm 2017-ExcelFIT-ormos.blg
	rm 2017-ExcelFIT-ormos.bbl
	rm 2017-ExcelFIT-ormos.log
	rm 2017-ExcelFIT-ormos.out
	#rm 2017-ExcelFIT-ormos.toc