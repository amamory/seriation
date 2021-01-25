set terminal pngcairo size 800, 600
set datafile separator "\t"
set multiplot
set tics font "name[,10]"
set yrange [0:1]
set xlabel "Tempo"
set ylabel "Custo"
set size 1,1
set origin 0,0
set key samplen 2 bottom left
plot filename using 1:5 title "H/H0" with lines lt 2
set key default
set xtics offset 0,graph 0.05
set tics font "name[,7]"
set yrange [0:*]
unset xlabel
unset ylabel
set format y "%.e"
set size 0.41,0.35
set origin 0.54,0.3
set key samplen 2 font ",8" top right
plot filename using 1:6 title "Temperatura/H0" with lines lt 1
set key default
set yrange [0:*]
unset xlabel
unset ylabel
set format y "%.e"
set size 0.41,0.35
set origin 0.54,0.6
set key samplen 2 font ",8" bottom right
plot filename using 1:4 title "Trocas" with lines lt 3
unset multiplot	 
