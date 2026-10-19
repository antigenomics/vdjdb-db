set border 3 lw 0.7
set tics nomirror out
set yrange [0:1]
set ytics 0,0.25,1
set xrange [0.5:3.5]
set xtics ("AP" 1, "Accuracy" 2, "Balanced\naccuracy" 3)
set ylabel "Held-out performance" offset 1,0
set grid ytics lc rgb "#DDDDDD" dt 3
set key inside left bottom vertical samplen 1 spacing 1.1
set multiplot layout 1,2 margins 0.13,0.98,0.24,0.83 spacing 0.10,0.05
do for [ep in "GLCTLVAML YLQPRTFLL"] {
  set title (ep eq "GLCTLVAML" ? "a  GLC (193 interactions)" : "b  YLQ (421 interactions)") offset 0,0.2
  plot "plot.tsv" using (strcol(1) eq ep ? (($0-1)%3+1)-0.13 : 1/0):3 with points pt 6 ps 1.1 lc rgb "#E69F00" title "ipTM", \
       "" using (strcol(1) eq ep ? (($0-1)%3+1) : 1/0):4 with points pt 4 ps 1.1 lc rgb "#009E73" title "TRA TCREmp", \
       "" using (strcol(1) eq ep ? (($0-1)%3+1)+0.13 : 1/0):5 with points pt 7 ps 1.1 lc rgb "#0072B2" title "Record score v3"
  unset key
  unset ylabel
}
unset multiplot
