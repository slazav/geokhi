#!/usr/bin/gnuplot

#set terminal "fig" metric
#set output "plot_t2.fig"

set nokey
plot [:4] \
 "fit02.dat" using ($1*1000):2 with lines,\
 "fit04.dat" using ($1*1000):2 with lines,\
 "fit06.dat" using ($1*1000):2 with lines,\
 "fit08.dat" using ($1*1000):2 with lines,\
0
pause -1

#set output "plot_t3.fig"

plot [:4] \
 "fit02.dat" using ($1*1000):3 with lines,\
 "fit04.dat" using ($1*1000):3 with lines,\
 "fit06.dat" using ($1*1000):3 with lines,\
 "fit08.dat" using ($1*1000):3 with lines,\
0
pause -1

#set output "plot_t4.fig"

plot [:4] \
 "fit02.dat" using ($1*1000):4 with lines,\
 "fit04.dat" using ($1*1000):4 with lines,\
 "fit06.dat" using ($1*1000):4 with lines,\
 "fit08.dat" using ($1*1000):4 with lines,\
0

pause -1
