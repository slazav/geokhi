#!/usr/bin/gnuplot

unset key

#set terminal fig metric
#set output "u1.fig"

plot [0:0.003]\
  'sol_u_1.dat' using 2:3 with lines lc 3,\
  'sol_u_2.dat' using 2:3 with lines lc 3,\
  'sol_u_3.dat' using 2:3 with lines lc 3,\
  'sol_u_4.dat' using 2:3 with lines lc 3,\
  'sol_u_5.dat' using 2:3 with lines lc 3,\
0
pause -1

#set output "u2.fig"

plot [0:0.003] []\
  'sol_q_1.dat' using 2:3 with lines lc 3,\
  'sol_q_2.dat' using 2:3 with lines lc 3,\
  'sol_q_3.dat' using 2:3 with lines lc 3,\
  'sol_q_4.dat' using 2:3 with lines lc 3,\
  'sol_q_5.dat' using 2:3 with lines lc 3,\
0
pause -1

plot [0:0.003] []\
  'sol_t_1.dat' using 2:3 with lines lc 3,\
  'sol_t_2.dat' using 2:3 with lines lc 3,\
  'sol_t_3.dat' using 2:3 with lines lc 3,\
  'sol_t_4.dat' using 2:3 with lines lc 3,\
  'sol_t_5.dat' using 2:3 with lines lc 3,\
0
pause -1

