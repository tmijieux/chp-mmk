#!/usr/bin/gnuplot

set zrange [ -0.2 : 1.0 ]
do for [i=1:2000:3]{ 
	splot "sol/sol".i.".dat.0" w p title "sol".i.".0", \
              "sol/sol".i.".dat.1" w p title "sol".i.".1", \
              "sol/sol".i.".dat.2" w p title "sol".i.".2", \
              "sol/sol".i.".dat.3" w p title "sol".i.".3"
        pause 0.005 
}
pause -1
