set key default
set xtics offset 0,graph 0.02
set tics font "name[,7]"
set yrange [0:1]
unset xlabel
unset ylabel
set size nosquare 0.48,0.36
set origin 0.49,0.06
set key samplen 2 bottom left
set object 1 rectangle from graph -0.1,-0.1 to graph 1.08,1.08 behind fc rgb "#FFFFFF" fs noborder
plot filename using 1:5 title "H/H0" with lines lt 2
unset multiplot
