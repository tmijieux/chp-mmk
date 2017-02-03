#!/usr/bin/gnuplot

set size ratio -1
set xlabel "nombre de processeurs"
set ylabel "speedup"

plot "speedup.data" using 1:2 title "speedup", \
    x with lines title "x=y"

pause -1
