set terminal postscript landscape color enhanced "Times-Roman" 16

set output 'bands_etb_linear.eps'

set autoscale

set ytic auto

set xtic auto

unset xlabel  #"k, {/Symbol p}/a" font "Times-Roman, 20"

set ylabel "energy, eV" font "Times-Roman, 20"

# set key right top	#set the legend

unset key

plot 'graphene_bands_etb_plot.dat' u 1:4 with line linewidth 2, \
'graphene_bands_etb_plot.dat' u 1:5 with line linewidth 2, \
'graphene_bands_etb_plot.dat' u 1:6 with line linewidth 2, \
'graphene_bands_etb_plot.dat' u 1:7 with line linewidth 2, \
'graphene_bands_etb_plot.dat' u 1:8 with line linewidth 2, \
'graphene_bands_etb_plot.dat' u 1:9 with line linewidth 2, \
'graphene_bands_etb_plot.dat' u 1:10 with line linewidth 2, \
'graphene_bands_etb_plot.dat' u 1:11 with line linewidth 2
