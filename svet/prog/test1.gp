#!/usr/bin/gnuplot

#set terminal fig metric
#set output "a.fig"

pi=3.1415926

n1=1.0
n2=1.5
a1=45.0 * pi/180.0
a2=asin(sin(a1)*n1/n2)
b=0.119

f(x) = tan(pi/4.0-a2)*x + b

unset key
print a2*180/pi
plot [-3:5] [-2.5:3]\
 "test1.txt" using 1:4 with lines,\
 -x, f(x), 0 lc 3
pause -1
