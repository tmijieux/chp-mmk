#!/usr/bin/gnuplot

set xrange[0:1]
set yrange[0:1]

splot "numeric.dat.1.0" using 1:2:3,  \
	"numeric.dat.1.1" using 1:2:3, \
	"numeric.dat.1.2" using 1:2:3, \
	"numeric.dat.1.3" using 1:2:3, \
	x*(1-x)*y*(1-y) w p

pause -1
