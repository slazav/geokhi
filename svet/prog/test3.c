#include <stdio.h>
#include <math.h>

// ����������� ����� �� ����� n0 + dn * (1.0 + tanh(k*(x+y)))/2.0

double n0=1.0;
double dn=0.8;
double k=1.0;

double n(double x, double y){
  return n0 + dn * (1.0 + tanh(k*y))/2.0;
}

double nx(double x, double y){
  return 0;
}

double ny(double x, double y){
  double c = cosh(k*y);
  return k*dn/2.0/c/c;
}


int
main(){
  double x1 = -10;
  double x2 = 10;
  double dx = 0.01;

  double y1 = 5;
  double y2 = 10;
  double dy = 0.5;

  double x,y;
  for (y=y1; y<y2; y+=dy){
    double f, fp=y, fpp=y+dx;
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
