#include <stdio.h>
#include <complex.h>
#include "morton_order.h"
#include "quadtree.h"



int main(){
  int i;
  int p_index;
  int c_index[4];
  int level;

  unsigned long pos[2];

  unsigned long x;

  double array1[2],array2[2];
  
  double a,b;
  double complex c;

  a = 2.0;
  b = 4.0;
  c = a+4*I;
  


  p_index = 2;
  getChild(p_index,c_index);
  for (i=0; i<4; i++){
    printf("get child(p=2): %i \n",c_index[i]);
  }

  i = getParent(19);
  printf("get parent(c=21): %i \n",i);

                     
  level = getLevel(13);
  printf("getLevel(l=13): %i \n",level);
  i = global2local(13,level);
  printf("global2local(13):%i \n",i); 
  p_index = local2global(i,level);
  printf("global2local(%i):%i \n",i,p_index); 
  
  
  index2pos(45,pos);
  printf("index2pos(45): %i,%i \n",pos[0],pos[1]);


  x = pos2index(3,6);
  printf("pos2index(3,6): %lu\n",x);

  
  array1[0] = 0;
  array1[1] = 0;
  array2[0] = 1;
  array2[1] = -1;
  printf("quadrant_pos:%i\n",quadrant_pos(array2,array1));



  return 0;
}
