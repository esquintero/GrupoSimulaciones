set terminal png size 500,800
set grid
# Punto a
set output 'SJ_a.png'
set multiplot
set title 'Sol-Jupiter'
plot "S-J_original_axis.dat" using 3:4 title "Jupiter" w l, "S-J_original_axis.dat" using 1:2 title "Sol" w l
set origin .25,.25
set size square .5,.5
clear
unset key
unset object
plot "S-J_original_axis.dat" using 1:2 title "Sol" w l linecolor "dark-green"
unset multiplot
unset output
unset key
unset object
unset origin
unset size
unset title
# Punto b
set output "SJ_b.png"
set multiplot
set title 'Sol-Jupiter_rotado.png'
set xtics -100,200,1100
plot "S-J_rotated_axis.dat" using 3:4 title "Jupiter" w l, "S-J_rotated_axis" using 1:2 title "Sol" w l
set origin .2,.2
set size .2,.5
clear
unset key
unset object
unset xtics
set xtics -2,2
plot "S-J_rotated_axis.dat" using 1:2 w l linecolor "dark-green"
unset multiplot
unset output
unset object
unset origin
unset size
unset title

