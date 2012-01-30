#!/usr/bin/gnuplot

unset key

#set terminal fig metric
#set output "u1.fig"

plot [0:0.003]\
  'sol_u1.dat' using 2:3 with lines lc 3,\
  'sol_u2.dat' using 2:3 with lines lc 3,\
  'sol_u3.dat' using 2:3 with lines lc 3,\
  'sol_u4.dat' using 2:3 with lines lc 3,\
  'sol_u5.dat' using 2:3 with lines lc 3,\
0
pause -1

#set output "u2.fig"

plot [0:0.003] []\
  'sol_q1.dat' using 2:3 with lines lc 3,\
  'sol_q2.dat' using 2:3 with lines lc 3,\
  'sol_q3.dat' using 2:3 with lines lc 3,\
  'sol_q4.dat' using 2:3 with lines lc 3,\
  'sol_q5.dat' using 2:3 with lines lc 3,\
0
pause -1

plot [0:0.003] []\
  'sol_t1.dat' using 2:3 with lines lc 3,\
  'sol_t2.dat' using 2:3 with lines lc 3,\
  'sol_t3.dat' using 2:3 with lines lc 3,\
  'sol_t4.dat' using 2:3 with lines lc 3,\
  'sol_t5.dat' using 2:3 with lines lc 3,\
0
pause -1

