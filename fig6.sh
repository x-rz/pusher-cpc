#!/bin/bash
size="12"
font="times"
f_name="fig6"
for str in a b c d
do
	echo "#define ${str}_txt" > mydef.h
	g++ ${f_name}.cc -O3 -mavx -mfma
	./a.out > ./${str}.txt
done
gnuplot << EOF
set term pdf font "${font},${size}" lw 2 size 6.4,6.4/3
set output "./${f_name}.pdf"
if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .95
if (!exists("MP_BOTTOM")) MP_BOTTOM = .20
if (!exists("MP_TOP"))    MP_TOP = .9
if (!exists("MP_xGAP"))   MP_xGAP = 0.08
if (!exists("MP_yGAP"))   MP_yGAP = 0.0
set multiplot layout 1,2 columnsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_xGAP, MP_yGAP
set xrange [0:6.e-4]
set yrange [0e2:4e3]
set xtics (0.,2.e-4,4.e-4,6.e-4)
set xtics format "%1.0tx10^{%S}"
set ytics ("0" 0,"1000" 1000.,"2000" 2000,"3000" 3000,"4000" 4000,"5000" 5000)
set xlabel 'time[second]'
set ylabel 'Lorentz factor'
unset key
set label 1 "(a)" at graph 0.06,0.9
plot "./a.txt" u 1:2 w l lt 1 lc rgb "gray", "./b.txt" u 1:2 w p lt 6 ps 1. lc rgb "black"
set yrange [0e2:1e4]
set ytics ("0" 0,"3000" 3000.,"6000" 6000,"9000" 9000)
set label 1 "(b)" at graph 0.06,0.9
plot "./c.txt" u 1:2 w l lt 1 lc rgb "gray", "./d.txt" u 1:2 w p lt 6 ps 1. lc rgb "black"
unset multiplot
EOF
rm a.out mydef.h a.txt b.txt c.txt d.txt
