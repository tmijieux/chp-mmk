#!/usr/bin/gnuplot

splot "numerique.dat" using 1:2:3, \
      "exacte.dat" using 1:2:3
pause -1
