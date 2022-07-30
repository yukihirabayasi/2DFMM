#include <stdio.h>
#include <stdlib.h>
#include "morton_order.h"
#include "quadtree.h"
#include "fmm.h"

void read_s_file(Vortex *vor);

int main(){
  int max_level = 10;
  int nterms = 12;
  double theta = 0.5;
  double r_min;
  QuadTree qtree;
  Vortex *vor;

  read_s_file(vor);
  qtree = make_tree(vor, nterms, theta, r_min);
  S2M();
  M2M();
  

  sleep(30);

  biot_tree();
}

// void read_s_file(Vortex *vor){
//   FILE *s_file;
//   char fname[] = "../s_file/s_00344_step.dat";
//   double dummy1,dummy2;
//   double nvor_b = 3;
//   s_file = fopen(fname,"r");
//   if (s_file == NULL) {
//     printf("%s not open.\n",fname);
//     return;
//   }

//   fscanf(s_file,"%f,%f",&dummy1,&dummy2);
//   printf("%f,%f\n",dummy1,dummy2);

//   fclose(s_file);
//   vor = (Vortex*)malloc(sizeof(Vortex) * nvor_b);
// }

