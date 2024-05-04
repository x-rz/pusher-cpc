#!/bin/bash
flag="-mavx -mfma"
echo -e "################ performance test ################" |tee perf.dta
echo -e "implementation\tO2\t\tSD\t\tO3\t\tSD\t\tOfast\t\tSD" |tee -a perf.dta
echo -ne "\"imp. 1\"\t" |tee -a perf.dta 
echo "#define imp1" > mydef.h
for opt in O2 O3 Ofast
do
	g++ fig9.cc -${opt} ${flag} && ./a.out |tee -a perf.dta
done
echo -ne "\n\"imp. 2\"\t" |tee -a perf.dta 
echo "#define imp2" > mydef.h
for opt in O2 O3 Ofast
do
	g++ fig9.cc -${opt} ${flag} && ./a.out |tee -a perf.dta
done
echo -ne "\noptimized\t" |tee -a perf.dta
echo "#define opti" >mydef.h
for opt in O2 O3 Ofast
do
	g++ fig9.cc -${opt} ${flag} && ./a.out |tee -a perf.dta
done
echo -ne "\n"
size="12"
font="times"
gnuplot << EOF
set term pdf font "${font},${size}" lw 2 size 6.4/2,6.4/3
set output "fig9.pdf"
if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .9
if (!exists("MP_BOTTOM")) MP_BOTTOM = .20
if (!exists("MP_TOP"))    MP_TOP = .9
if (!exists("MP_xGAP"))   MP_xGAP = 0.0
if (!exists("MP_yGAP"))   MP_yGAP = 0.0
set style data histogram
set xrange [-0.8:2.8]
set yrange [0:350] noextend
set y2range [0:350] noextend
set boxwidth 1.
set style fill solid 1.00 border lt -1
set style histogram errorbars lw 1.0
set datafile columnhead separator tab
set format y "%.1f"
set key box opaque
set ylabel 'time/part/step[ns]'
set key nobox
#set xlabel 'analytical algorithm'
set ytics ("0" 0,"100" 100.,"200" 200,"300" 300)
set y2label 'processor cycles/part/step'
set y2tics ("0" 0,"300" 100.,"600" 200,"900" 300)
plot 'perf.dta' using 2:3 fs solid 0.5 lw 1.0 ti 'O2', \
     '' using 4:5:xticlabel(1) fs solid 0.5 lw 1.0 ti 'O3',\
     '' using 6:7:xticlabel(1) fs solid 0.5 lw 1.0 ti 'Ofast'
EOF
rm a.out mydef.h perf.dta
