.PHONY: all
all: workflow.svg workflow.pdf performance.svg performance.pdf ../motif-methods.svg ../dcode-method.svg

workflow.svg: workflow.dot ../theme_svg.py
	dot -Tsvg $< -o $@
	python3 ../theme_svg.py $@

workflow.pdf: workflow.dot
	dot -Tpdf $< -o $@

plot.tsv: metrics.tsv plot.R
	Rscript plot.R

performance.svg performance.pdf &: plot.tsv performance.gp performance-body.gp ../theme_svg.py
	gnuplot performance.gp
	python3 ../theme_svg.py performance.svg

../%.svg: ../%.dot ../theme_svg.py
	dot -Tsvg $< -o $@
	python3 ../theme_svg.py $@
