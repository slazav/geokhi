#!/usr/bin/gnuplot

unset key

#set terminal fig metric
#set output "u1.fig"

plot [0:0.003]\
  'solz_u1.dat' using 2:3 with lines lc 3,\
  'solz_u2.dat' using 2:3 with lines lc 3,\
  'solz_u3.dat' using 2:3 with lines lc 3,\
  'solz_u4.dat' using 2:3 with lines lc 3,\
  'solz_u5.dat' using 2:3 with lines lc 3,\
0
pause -1

#set output "u2.fig"

plot [0:0.003] []\
  'solz_q1.dat' using 2:3 with lines lc 3,\
  'solz_q2.dat' using 2:3 with lines lc 3,\
  'solz_q3.dat' using 2:3 with lines lc 3,\
  'solz_q4.dat' using 2:3 with lines lc 3,\
  'solz_q5.dat' using 2:3 with lines lc 3,\
0
pause -1

plot [0:0.003] []\
  'solz_t1.dat' using 2:3 with lines lc 3,\
  'solz_t2.dat' using 2:3 with lines lc 3,\
  'solz_t3.dat' using 2:3 with lines lc 3,\
  'solz_t4.dat' using 2:3 with lines lc 3,\
  'solz_t5.dat' using 2:3 with lines lc 3,\
0
pause -1

plot \
  'solr_t1.dat' using 1:3 with lines lc 3,\
  'solr_t2.dat' using 1:3 with lines lc 3,\
  'solr_t3.dat' using 1:3 with lines lc 3,\
  'solr_t4.dat' using 1:3 with lines lc 3,\
  'solr_t5.dat' using 1:3 with lines lc 3,\
0
pause -1

