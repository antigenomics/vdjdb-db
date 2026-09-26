# The clustering bake-off, one figure per chain: six instruments against retention.
#
#   gnuplot -e "chain='TRB'" docs/tuning/tuning.gp
#
# Reads docs/tuning/scorecard_plot.tsv (committed) and writes out/reports/tuning/tuning_<chain>.svg.
# Retention is the x axis of every panel because it is the axis every instrument is confounded with:
# lift and F1 fall as it rises, Q and purity and coverage rise, and `trivial` -- one cluster per
# epitope, nothing excluded -- sits at retention 1.0 having clustered nothing.
#
# Colours are ColorBrewer Set1 (qualitative, 9-class), per CLAUDE.md section 7.
if (!exists("chain")) chain = "TRB"
F = "docs/tuning/scorecard_plot.tsv"
# 1 gene  2 algo  3 config  4 retention  5 lift  6 f1  7 q  8 purity  9 precision
# 10 epitopes  11 perc_med  12 cids  13 admissible
ALGOS = "dbscan hdbscan-eom hdbscan-leaf lumbermark hybrid hybrid-recruited legacy trivial"
NAMES = "DBSCAN-shipped HDBSCAN-eom HDBSCAN-leaf Lumbermark hybrid-enriched \
hybrid-recruited legacy-release trivial-do-nothing"

set terminal svg size 1500,880 font "Helvetica,12" background rgb "white"
set output "out/reports/tuning/tuning_".chain.".svg"
set datafile separator "\t"
set datafile missing "null"

set style line 1 lc rgb "#e41a1c" pt 7  ps 1.1          # DBSCAN, the shipped path
set style line 2 lc rgb "#377eb8" pt 9  ps 1.2          # HDBSCAN, excess of mass
set style line 3 lc rgb "#a65628" pt 11 ps 1.2          # HDBSCAN, every condensed-tree leaf
set style line 4 lc rgb "#4daf4a" pt 5  ps 1.0          # Lumbermark, ungated
set style line 5 lc rgb "#984ea3" pt 13 ps 1.3          # TCRNET enriched gate + Lumbermark
set style line 6 lc rgb "#ff7f00" pt 2  ps 1.3 lw 2     # TCRNET recruited gate + Lumbermark
set style line 7 lc rgb "black"   pt 12 ps 2.0 lw 3     # the legacy release
set style line 8 lc rgb "#999999" pt 6  ps 2.0 lw 3     # the instrument's blind spot

set border 3 lw 1
set tics nomirror out
set grid ytics lc rgb "#dddddd"
set xlabel "retention (fraction of cohort clonotypes clustered)"
set xrange [0:1.03]
unset key

set multiplot layout 2,3 title chain."  --  every instrument against retention" font ",15"

set logscale y
set ylabel "independent-study lift"
set title "lift falls monotonically with retention"
plot for [i=1:words(ALGOS)] F using \
  (strcol(1) eq chain && strcol(2) eq word(ALGOS,i) ? $4 : NaN):(column(5)) ls i
unset logscale y

set ylabel "F1, clustered predicts replicated"
set title "F1: the statistic with no retention bias"
plot for [i=1:words(ALGOS)] F using \
  (strcol(1) eq chain && strcol(2) eq word(ALGOS,i) ? $4 : NaN):(column(6)) ls i

set ylabel "Q = 2hp/(h+p)"
set title "Q rises toward the do-nothing corner"
plot for [i=1:words(ALGOS)] F using \
  (strcol(1) eq chain && strcol(2) eq word(ALGOS,i) ? $4 : NaN):(column(7)) ls i

set ylabel "purity (epitope)"
set title "purity: the only axis trivial fails"
plot for [i=1:words(ALGOS)] F using \
  (strcol(1) eq chain && strcol(2) eq word(ALGOS,i) ? $4 : NaN):(column(8)) ls i

set ylabel "epitopes with at least one cluster"
set title "coverage is also maximal for doing nothing"
plot for [i=1:words(ALGOS)] F using \
  (strcol(1) eq chain && strcol(2) eq word(ALGOS,i) ? $4 : NaN):(column(10)) ls i

set key outside center bottom horizontal maxrows 2 samplen 1 spacing 1.1 font ",11"
set ylabel "median per-epitope percolation"
set title "percolation: largest cluster's share, lower is better"
plot for [i=1:words(ALGOS)] F using \
  (strcol(1) eq chain && strcol(2) eq word(ALGOS,i) ? $4 : NaN):(column(11)) \
  ls i title word(NAMES,i)

unset multiplot
