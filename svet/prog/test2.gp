#!/usr/bin/gnuplot

pi=3.1415926

n1=1.0
n2=1.40
a1=45.0 * pi/180.0
a2=asin(sin(a1)*n1/n2)

print a2*180/pi
plot [-3:5] [-2.5:2]\
 "test2.txt" using 1:4 with lines,\
 -x, tan(pi/4.0-a2)*x
pause -1
