.PHONY: all
all: workflow.svg workflow.pdf performance.svg performance.pdf

workflow.svg: workflow.dot
	dot -Tsvg $< -o $@

workflow.pdf: workflow.dot
	dot -Tpdf $< -o $@

plot.tsv: metrics.tsv plot.R
	Rscript plot.R

performance.svg performance.pdf &: plot.tsv performance.gp performance-body.gp
	gnuplot performance.gp
