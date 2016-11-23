#!/usr/bin/gnuplot

set xrange[0:10]
set yrange[0:10]

splot "numeric.dat.1.0" using 1:2:3,  \
	"numeric.dat.1.1" using 1:2:3, \
	"numeric.dat.1.2" using 1:2:3, \
	"numeric.dat.1.3" using 1:2:3, \
      	sin(x)+cos(y)

pause -1
