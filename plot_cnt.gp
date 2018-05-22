set terminal postscript landscape color enhanced "Times-Roman" 16
#set terminal png
set output 'bands_cnt.eps'

set autoscale
set ytic auto
set xtic auto
#set xrange [-3.5:3.5]


# set logscale y
# set logscale x

set xlabel "k" font "Times-Roman, 20"
set ylabel "energy, eV" font "Times-Roman, 20"

unset key

set style line 1 lt 1 lc rgb "#00008b" lw 2 pt 13 ps 0.5

plot 'cnt_bands_etb.dat' u 1:4 w l ls 1, \
'cnt_bands_etb.dat' u 1:5 w l ls 1, \
'cnt_bands_etb.dat' u 1:6 w l ls 1, \
'cnt_bands_etb.dat' u 1:7 w l ls 1, \
'cnt_bands_etb.dat' u 1:8 w l ls 1, \
'cnt_bands_etb.dat' u 1:9 w l ls 1, \
'cnt_bands_etb.dat' u 1:10 w l ls 1, \
'cnt_bands_etb.dat' u 1:11 w l ls 1
