#!/usr/bin/gnuplot

#set size ratio -1
set term png
set output "speedup2.png"
set xlabel "recouvrement"
set ylabel "speedup"

plot "speedup.data" using 1:2 title "speedup", \
     1 with lines title "y=1"
