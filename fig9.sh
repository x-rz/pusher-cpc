#!/bin/sh
echo -e "#========analytical algorithm========" |tee perf_anal.dta
echo -e "optim\ttau\t\tSD\t\tupdate\t\tSD" |tee -a perf_anal.dta
for opt in O2 O3 Ofast
do
	echo "#define anal_t" > mydef.h
	echo -ne "${opt}\t" |tee -a perf_anal.dta
	g++ fig9.cc -${opt} && ./a.out |tee -a perf_anal.dta
	echo -ne "\t" |tee -a perf_anal.dta
	echo "#define anal_p" > mydef.h
	g++ fig9.cc -${opt} && ./a.out |tee -a perf_anal.dta
	echo -ne "\n" |tee -a perf_anal.dta
done
echo -e "#==========this work(Imp. 1)==========" |tee perf_imp1.dta
echo -e "optim\ttau\t\tSD\t\tupdate\t\tSD" |tee -a perf_imp1.dta
for opt in O2 O3 Ofast
do
	echo "#define imp1_t" > mydef.h
	echo -ne "${opt}\t" |tee -a perf_imp1.dta
	g++ fig9.cc -${opt} && ./a.out |tee -a perf_imp1.dta
	echo -ne "\t" |tee -a perf_imp1.dta
	echo "#define imp1_p" > mydef.h
	g++ fig9.cc -${opt} && ./a.out |tee -a perf_imp1.dta
	echo -ne "\n" |tee -a perf_imp1.dta
done
echo -e "#==========this work(Imp. 2)==========" |tee perf_imp2.dta
echo -e "optim\ttau\t\tSD\t\tupdate\t\tSD" |tee -a perf_imp2.dta
for opt in O2 O3 Ofast
do
	echo "#define imp2_t" > mydef.h
	echo -ne "${opt}\t" |tee -a perf_imp2.dta
	g++ fig9.cc -${opt} && ./a.out |tee -a perf_imp2.dta
	echo -ne "\t" |tee -a perf_imp2.dta
	echo "#define imp2_p" > mydef.h
	g++ fig9.cc -${opt} && ./a.out |tee -a perf_imp2.dta
	echo -ne "\n" |tee -a perf_imp2.dta
done
size="12"
font="times"
gnuplot << EOF
set term pdf font "${font},${size}" lw 2 size 6.4,6.4/3
set output "fig9.pdf"
if (!exists("MP_LEFT"))   MP_LEFT = .1
if (!exists("MP_RIGHT"))  MP_RIGHT = .9
if (!exists("MP_BOTTOM")) MP_BOTTOM = .20
if (!exists("MP_TOP"))    MP_TOP = .9
if (!exists("MP_xGAP"))   MP_xGAP = 0.0
if (!exists("MP_yGAP"))   MP_yGAP = 0.0
set multiplot layout 1,3 columnsfirst margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_xGAP, MP_yGAP
set style data histogram
set yrange [0:500] noextend
set y2range [0:500] noextend
set boxwidth 1.
set style fill solid 1.00 border lt -1
set style histogram errorbars lw 1
set datafile columnhead separator tab
set format y "%.1f"
set key box opaque
set ylabel 'time/part/step[ns]'
set key nobox
set xlabel 'analytical algorithm'
set ytics ("0" 0,"200" 200.,"400" 400,"600" 600)
set label 1 "(a)" at graph 0.1,0.9
plot 'perf_anal.dta' using 2:3 fs solid 0.5 lw 1 ti 'computing {/Symbol t}', \
     '' using 4:5:xticlabel(1) fs solid 0.0 lw 1 ti 'updating u^{/Symbol m},x^{/Symbol m}'
unset ylabel
unset ytics
set xlabel 'this work(Imp. 1)'
set label 1 "(b)" at graph 0.1,0.9
plot 'perf_imp1.dta' using 2:3 fs solid 0.5 lw 1 ti 'computing {/Symbol t}', \
     '' using 4:5:xticlabel(1) fs solid 0.0 lw 1 ti 'updating u^{/Symbol m},x^{/Symbol m}'
set xlabel 'this work(Imp. 2)'
set label 1 "(c)" at graph 0.1,0.9
set y2label 'processor cycles/part/step'
set y2tics ("0" 0,"600" 200.,"1200" 400,"1800" 600)
plot 'perf_imp2.dta' using 2:3 fs solid 0.5 lw 1 ti 'computing {/Symbol t}', \
     '' using 4:5:xticlabel(1) fs solid 0.0 lw 1 ti 'updating u^{/Symbol m},x^{/Symbol m}'
unset multiplot
EOF
rm a.out mydef.h perf_anal.dta perf_imp1.dta perf_imp2.dta
