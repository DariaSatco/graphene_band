set terminal postscript landscape color enhanced "Times-Roman" 16

set output 'bands.eps'

set autoscale

set ytic auto

set xtic auto

set pm3d

# set logscale y

# set logscale x

set xlabel "k_x, {/Symbol p}/a" font "Times-Roman, 20"

set ylabel "k_y, {/Symbol p}/a" font "Times-Roman, 20"

set zlabel "Energy, eV" font "Times-Roman, 20"

# set key right top	#set the legend

unset key

splot 'graphene_bands.dat' u 1:2:3 with pm3d linewidth 2, \
'graphene_bands.dat' u 1:2:4 with pm3d linewidth 2
