#!/usr/bin/gnuplot

splot "numeric.dat" using 1:2:3, \
	sin(x)+cos(y)
#      "exact.dat" using 1:2:3
pause -1
