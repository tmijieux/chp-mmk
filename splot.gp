#!/usr/bin/gnuplot

do for [i=1:2000]{ 
	splot "sol/sol".i.".dat" ; pause 0.05
}
