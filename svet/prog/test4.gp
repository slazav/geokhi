#!/usr/bin/gnuplot

#set terminal fig metric
#set output "a.fig"

pi=3.1415926

unset key
plot [-2:6] [-3:3]\
 "test4.txt" using 1:4 with lines,\
  exp(-x*x) lc 3,0 lc 3

pause -1
