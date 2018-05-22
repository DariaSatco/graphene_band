set terminal postscript landscape color enhanced "Times-Roman" 16

set output 'bands_etb_linear1.eps'

set autoscale
set ytic auto
set xtics ("{/Symbol G}" 0.0, "K" 0.770, "M" 1.155, "{/Symbol G}" 1.821)

# set xlabel  "k, {/Symbol p}/a" font "Times-Roman, 20"

set ylabel "energy, eV" font "Times-Roman, 20"

# set key right top	#set the legend

unset key
set style line 1 lt 1 lc rgb "#00008b" lw 2 pt 7 ps 1
set style line 2 lt 0 lc rgb "#a0a0a0" lw 1 pt 7 ps 1
set style line 3 lt 1 lc rgb "#00ff00" lw 1 pt 1 ps 1

f(x)=0.0

set grid xtics ls 2

plot 'graphene_bands_etb_plot.dat' u 1:4 with line ls 1, \
'graphene_bands_etb_plot.dat' u 1:5 with line ls 1, \
'graphene_bands_etb_plot.dat' u 1:6 with line ls 1, \
'graphene_bands_etb_plot.dat' u 1:7 with line ls 1, \
'graphene_bands_etb_plot.dat' u 1:8 with line ls 1, \
'graphene_bands_etb_plot.dat' u 1:9 with line ls 1, \
'graphene_bands_etb_plot.dat' u 1:10 with line ls 1, \
'graphene_bands_etb_plot.dat' u 1:11 with line ls 1, \
'graphene_bands_tb_plot.dat' u 1:4 w l ls 3, \
'graphene_bands_tb_plot.dat' u 1:5 w l ls 3, \
f(x) with line ls 2
