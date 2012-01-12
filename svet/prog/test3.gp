#!/usr/bin/gnuplot

#set terminal fig metric
#set output "a.fig"

pi=3.1415926

n1=1.0
n2=1.5
a1=45.0 * pi/180.0
a2=asin(sin(a1)*n1/n2)
b=6

f(x) = -x + b

print a2*180/pi
plot [-6:5] [-1:]\
 "test3.txt" using 1:4 with lines,\
 0, x+1.03, -x-0.5
pause -1
