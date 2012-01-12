#include <stdio.h>
#include <math.h>

// рассеянье света на клине n = n0 + dn * exp(-x*x -y*y)

double n0=1;
double dn=0.4;
double k=2;

double n(double x, double y){
  return n0 + dn * ( 1+ tanh(k*(x+y)))/2.0;
}

double nx(double x, double y){
  double c = cosh(k*(x+y));
  return k*dn/2/c/c;
}

double ny(double x, double y){
  double c = cosh(k*(x+y));
  return k*dn/2/c/c;
}


int
main(){
  double x1 = -5;
  double x2 = 5;
  double dx = 0.01;

  double y1 = -2;
  double y2 = 2;
  double dy = 0.25;

  double x,y;
  for (y=y1; y<y2; y+=dy){
    double f, fp=y, fpp=y;
    for (x=x1; x<x2; x+=dx){
      double fx=(fp-fpp)/dx;
      double fxx=(ny(x,y)-fx*nx(x,y))*(1+fx*fx)/n(x,y);
      f=2*fp-fpp+fxx*dx*dx;

      printf("%10f %10f %10f %10f\n", x, y, n(x,y), f);
      fpp=fp; fp=f;
    }
    printf("\n");
  }
}