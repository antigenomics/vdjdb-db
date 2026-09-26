# The per-epitope distributions behind the pooled scorecard.
#
#   gnuplot docs/tuning/per_epitope.gp     ->  out/reports/tuning/per_epitope.svg
#
# Reads docs/tuning/per_epitope_ecdf.tsv. The top row is why a pooled retention figure is not a
# description of the database: on TRB the median epitope has under 1 % of its clonotypes clustered
# while the pooled figure is 32 %. The bottom row is percolation -- the largest cluster's share of
# each epitope's clustered clonotypes -- which is the failure a lift figure cannot see.
# Colours: ColorBrewer Set1.
F = "docs/tuning/per_epitope_ecdf.tsv"
# 1 gene  2 method  3 stat  4 value  5 frac
set terminal svg size 1180,780 font "Helvetica,12" background rgb "white"
set output "out/reports/tuning/per_epitope.svg"
set datafile separator "\t"

set style line 1 lc rgb "#e41a1c" lw 3      # TCRNET
set style line 2 lc rgb "#377eb8" lw 3 dt 2 # TCREMP
set border 3 lw 1
set tics nomirror out
set grid xtics ytics lc rgb "#dddddd"
set ylabel "fraction of epitopes at or below"
set yrange [0:1]
load "docs/tuning/pooled.gp"
unset key

set multiplot layout 2,2 title "per-epitope distributions, human, shipped files" font ",14"
set xrange [0:1]
do for [s in "retention percolation"] {
  do for [g in "TRA TRB"] {
    set title g."  --  ".s
    set xlabel s." within one epitope"
    unset arrow
    if (s eq "retention") {
      # The pooled retention, which is what a scorecard reports, against the distribution it
      # averages over. The gap between this rule and the median is the whole point of the panel.
      set arrow from (g eq "TRA" ? pooled_TRA_tcrnet : pooled_TRB_tcrnet), graph 0 \
                  to (g eq "TRA" ? pooled_TRA_tcrnet : pooled_TRB_tcrnet), graph 1 \
                  nohead lc rgb "#e41a1c" dt 3 lw 2
      set arrow from (g eq "TRA" ? pooled_TRA_tcremp : pooled_TRB_tcremp), graph 0 \
                  to (g eq "TRA" ? pooled_TRA_tcremp : pooled_TRB_tcremp), graph 1 \
                  nohead lc rgb "#377eb8" dt 3 lw 2
      set label 1 "dotted: pooled" at graph 0.55, 0.30 font ",10" tc rgb "#555555"
    } else {
      unset label 1
    }
    if (s eq "percolation" && g eq "TRB") {
      set key at graph 0.62, 0.22 samplen 3 spacing 1.2 font ",11"
    } else { unset key }
    plot for [i=1:2] F using \
      (strcol(1) eq g && strcol(3) eq s && strcol(2) eq word("tcrnet tcremp",i) ? $4 : NaN):5 \
      with steps ls i title word("TCRNET TCREMP",i)
  }
}
unset multiplot
