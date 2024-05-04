#!/bin/bash
size="12"
font="times"
f_name="fig2"
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
set multiplot layout 1,3 columnsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_xGAP, MP_yGAP
set xrange [0:1.e-4]
set yrange [8e2:2.2e3]
set xtics (0.,4.e-5,8.e-5,12.e-5)
set ytics ("1000" 1000.,"1500" 1500,"2000" 2000)
set xlabel 'time[second]'
set ylabel 'Lorentz factor'
unset key
set label 1 "(a)" at graph 0.1,0.9 font "${font},${size}"
plot "./a.txt" u 1:2 w l lt 1 lc rgb "gray", "./b.txt" u 1:2 w p lt 6 ps 1. lc rgb "black"
set label 1 "(b)" at graph 0.1,0.9 font "${font},${size}"
plot "./a.txt" u 1:2 w l lt 1 lc rgb "gray", "./c.txt" u 1:2 w p lt 6 ps 1. lc rgb "black"
set label 1 "(c)" at graph 0.1,0.9 font "${font},${size}"
plot "./a.txt" u 1:2 w l lt 1 lc rgb "gray", "./d.txt" u 1:2 w p lt 6 ps 1. lc rgb "black"
unset multiplot
EOF
rm a.out mydef.h a.txt b.txt c.txt d.txt
