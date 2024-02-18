#!/bin/bash
size="12"
font="times"
f_name="fig5"
gyr_r=5.689452e+05
for str in a b d
do
	echo "#define ${str}_txt" > mydef.h
	g++ ${f_name}.cc -O3  -mavx -mfma
	./a.out > ./${str}.txt
done
gnuplot << EOF
set term pdf font "${font},${size}" lw 2 size 6.4,6.4/3
set output "./${f_name}.pdf"
set size square
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
unset border
unset tics
unset key
set multiplot layout 1,2
set xrang[-0.2*${gyr_r}:2.2*${gyr_r}]
set yrang[-1.2*${gyr_r}:1.2*${gyr_r}]
set label 1 "(a)" at graph 0.1,0.1
plot "./a.txt" u 6:7 w l lt 1 lc rgb "gray", "./b.txt" u 6:7 w lp lt 6 ps 1. lc rgb "black", '-' u 1:2 w p lt 2 lc rgb "black"
0 0
e
set label 1 "(b)" at graph 0.1,0.1
plot "./a.txt" u 6:7 w l lt 1 lc rgb "gray", "./d.txt" u 6:7 w lp lt 6 ps 1. lc rgb "black", '-' u 1:2 w p lt 2 lc rgb "black"
0 0
e
EOF
rm a.out mydef.h a.txt b.txt d.txt
