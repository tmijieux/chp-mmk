#!/usr/bin/gnuplot

set zrange [-0.09:0.09]
do for [i=1:2000]{ 
	splot "sol/sol".i.".dat" w pm3d ; pause 0.05
}
