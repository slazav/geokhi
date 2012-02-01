#!/usr/bin/gnuplot

unset key

f1(x) = a1 - b1 * x**2
f2(x) = a2 - b2 * x**2
f3(x) = a3 - b3 * x**2
f4(x) = a4 - b4 * x**2
f5(x) = a5 - b5 * x**2

xr=0.000333*1e3

fit [0:xr] f1(x) 'solr_t1.dat' using ($1*1e3):($3*1e6) via a1,b1
b2=b1
fit [0:xr] f2(x) 'solr_t2.dat' using ($1*1e3):($3*1e6) via a2,b2
b3=b2
fit [0:xr] f3(x) 'solr_t3.dat' using ($1*1e3):($3*1e6) via a3,b3
b4=b3
fit [0:xr] f4(x) 'solr_t4.dat' using ($1*1e3):($3*1e6) via a4,b4
b5=b4
fit [0:xr] f5(x) 'solr_t5.dat' using ($1*1e3):($3*1e6) via a5,b5

plot [0:0.0011][0:]\
  'solr_t1.dat' using 1:3 with lines lc 3, f1(x*1e3)/1e6 lc 4,\
  'solr_t2.dat' using 1:3 with lines lc 3, f2(x*1e3)/1e6 lc 4,\
  'solr_t3.dat' using 1:3 with lines lc 3, f3(x*1e3)/1e6 lc 4,\
  'solr_t4.dat' using 1:3 with lines lc 3, f4(x*1e3)/1e6 lc 4,\
  'solr_t5.dat' using 1:3 with lines lc 3, f5(x*1e3)/1e6 lc 4,\
0

pause -1

plot [0:0.0008]\
  'solr_t1.dat' using 1:($3-f1($1*1e3)/1e6) with lines lc 3,\
  'solr_t2.dat' using 1:($3-f2($1*1e3)/1e6) with lines lc 3,\
  'solr_t3.dat' using 1:($3-f3($1*1e3)/1e6) with lines lc 3,\
  'solr_t4.dat' using 1:($3-f4($1*1e3)/1e6) with lines lc 3,\
  'solr_t5.dat' using 1:($3-f5($1*1e3)/1e6) with lines lc 3,\
0

print 0.0, a1, b1
print 0.5, a2, b2
print 1.0, a3, b3
print 2.0, a4, b4
print 4.0, a5, b5

pause -1

