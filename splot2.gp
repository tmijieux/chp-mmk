#!/usr/bin/gnuplot

do for [i=1:2000]{ 
	splot "sol/sol".i.".dat" w pm3d ; pause 0.05
}
