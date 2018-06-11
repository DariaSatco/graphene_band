set terminal postscript landscape color enhanced "Times-Roman" 16
#set terminal png
set output 'bands_cnt.eps'

set autoscale
set ytic auto
set xtic auto
#set xrange [-3.14:3.14]
#set yrange [-5.5:5.5]


# set logscale y
# set logscale x

set xlabel "k,  {/Symbol p}/T" font "Times-Roman, 20"
set ylabel "energy, eV" font "Times-Roman, 20"

unset key

set style line 1 lt 1 lc rgb "#00008b" lw 2 pt 13 ps 0.5  #dark-blue
set style line 2 lt 1 lc rgb "#ff1493" lw 2 pt 13 ps 0.5  #dark-pink
set style line 3 lt 1 lc rgb "#228b22" lw 2 pt 13 ps 0.5  #forest-green
set style line 4 lt 1 lc rgb "#9400d3" lw 2 pt 13 ps 0.5  #dark-violet
set style line 5 lt 1 lc rgb "#00ced1" lw 2 pt 13 ps 0.5  #dark-turquoise
set style line 6 lt 1 lc rgb "#306080" lw 2 pt 13 ps 0.5  #steelblue
set style line 7 lt 1 lc rgb "#00ff00" lw 2 pt 13 ps 0.5  #green
set style line 8 lt 1 lc rgb "#00ff00" lw 2 pt 13 ps 0.5

plot 'cnt_bands_etb.dat' u 1:4 w l ls 4, \
'cnt_bands_etb.dat' u 1:5 w l ls 2, \
'cnt_bands_etb.dat' u 1:6 w l ls 3, \
'cnt_bands_etb.dat' u 1:7 w l ls 1, \
'cnt_bands_etb.dat' u 1:8 w l ls 1, \
'cnt_bands_etb.dat' u 1:9 w l ls 3, \
'cnt_bands_etb.dat' u 1:10 w l ls 2 #, \
#'cnt_bands_etb.dat' u 1:($11 > 1. ? $11 : 1/0) w l ls 4
#'cnt_bands_etb.dat' u 1:11 w p ls 4 #, \
