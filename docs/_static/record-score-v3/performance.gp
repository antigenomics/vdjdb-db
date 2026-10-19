# Shared 0-1 scale; points compare the saved grouped held-out forecasts.
set encoding utf8
set datafile separator "\t"
set terminal svg size 840,340 font "Arial,16" noenhanced
set output "performance.svg"
load "performance-body.gp"
set terminal pdfcairo size 7,2.83 font "Arial,9" noenhanced
set output "performance.pdf"
load "performance-body.gp"
unset output
