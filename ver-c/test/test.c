#include <stdio.h>
#include "test.h"

void make_tree(int *nvor_b, 
               double poi_x[*nvor_b],
               double poi_y[*nvor_b], 
               double cir[*nvor_b], 
               double cor[*nvor_b])
{
  int i = 0;
  printf("make tree\n");
  qtree.max_level = 22;
  for (i=0;i<*nvor_b;i++){
      printf("%f\n",poi_x[i]);
  }
}

void biot_tree(double *x,double *y, double *u, double *v){
  printf("max level:%i\n",qtree.max_level);
  printf("x,y: %f %f\n",*x,*y);
  *u = *x + *y;
  *v = *x - *y;
  qtree.max_level = 12;
}