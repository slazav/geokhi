#!/usr/bin/gnuplot

#set terminal "fig" metric
#set output "plot_t1.fig"

set nokey
plot [:4] \
 "sol_s.dat" using ($2*1000):3 with lines,\
 "sol_m.dat" using ($2*1000):3 with lines,\
 "sol_l.dat" using ($2*1000):3 with lines,\
0
pause -1