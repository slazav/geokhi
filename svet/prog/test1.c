#include <stdio.h>
#include <math.h>

// рассеянье света на пузыре n = n0 + dn * exp(-x*x -y*y)

double n0=1;
double dn=0.2;

double n(double x, double y){
  return n0 + dn * exp(-x*x -y*y);
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
  double x2 = 5;
  double dx = 0.01;

  double y1 = -2;
  double y2 = 2;
  double dy = 0.25;

  double x,y;
  for (y=y1; y<y2; y+=dy){
    double f, fp=y, fpp=y;
    for (x=x1; x<x2; x+=dx){
      f=(fpp*(dx*nx(x,y) + 2*n(x,y)) - 4*fp*n(x,y) + 2*dx*dx*ny(x,y))/(dx*nx(x,y) - 2*n(x,y));
      printf("%10f %10f %10f %10f\n", x, y, n(x,y), f);
      fpp=fp; fp=f;
    }
    printf("\n");
  }
}