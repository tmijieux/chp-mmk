#!/usr/bin/gnuplot

splot "numeric.dat.1.0" using 1:2:3,  \
	"numeric.dat.1.1" using 1:2:3,  \
	"numeric.dat.1.2" using 1:2:3,  \
	"numeric.dat.1.3" using 1:2:3
pause -1
