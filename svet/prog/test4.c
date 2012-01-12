#include <stdio.h>
#include <math.h>

// рассеивание света на пятне  n0 + dn * exp(-x*x-y*y);

double n0=1.0;
double dn=0.1;

double n(double x, double y){
  return n0 + dn * exp(-x*x-y*y);
}

double nx(double x, double y){
  return - 2 * x * dn * exp(-x*x-y*y);
}

double ny(double x, double y){
  return - 2 * y * dn * exp(-x*x-y*y);
}


int
main(){
  double x1 = -5;
  double x2 = 10;
  double dx = 0.01;

  double y1 = -5;
  double y2 = 5;
  double dy = 0.2;

  double x,y;
  for (y=y1; y<y2; y+=dy){
    double f, fp=y, fpp=y;
    for (x=x1; x<x2; x+=dx){
      double fx=(fp-fpp)/dx;
      double fxx=(ny(x,fp)-fx*nx(x,fp))*(1+fx*fx)/n(x,fp);

      f=2*fp-fpp+fxx*dx*dx;

      printf("%10f %10f %10f %10f\n", x, y, n(x,y), f);

      fpp=fp; fp=f;
    }
    printf("\n");
  }
}
