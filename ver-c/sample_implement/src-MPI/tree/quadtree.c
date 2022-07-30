#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h> 
#include "quadtree.h"
#include "morton_order.h"
#include "fmm.h"

QuadTree qtree;
Vortex *vors;

void read_fmm_param(){
  FILE *file;
  file = fopen("../cal_cond/fmm_param.dat","r");
  if(file==NULL){
    printf("Error: Could not read FMM parameters.");
    return;
  }

  fscanf(file,"%i    nterms",&qtree.nterms);
  fscanf(file,"%i   num_p",&qtree.num_p);
  fscanf(file,"%lf theta",&qtree.theta);  
  fscanf(file,"%lf  r_min ! for Tree",&qtree.r_min);
  
  fclose(file);
  return;
}

int quadrant_pos(double pos[2], double c_pos[2]){
  int quad_index;
  quad_index = 0;
  if (pos[0] >= c_pos[0]){ quad_index++; }
  if (pos[1] <= c_pos[1]){ quad_index += 2; }
  return quad_index;
}

double maxval(double a[],size_t n){
  int i;
  double  max;
  max = a[0];
  for(i=1; i<n; i++){
    if(a[i] > max){
      max = a[i];
    }
  }
  // for (i=0;i<10;i++){printf("aa:%f",a[i]);}
  return max;
}

double minval(double a[],size_t n){
  int i;
  double min;
  min = a[0];
  for(i=1; i<n; i++){
    if(a[i] < min){
      min = a[i];
    }        
  }
  return min;
}

void calc_CenterPos(int quad_i, double *p_centerPos, double region, double *c_centerPos){
  // printf("\n\n%i",quad_i);
  if (quad_i%2 != 0) {
    // printf("x+\n");
    c_centerPos[0] = p_centerPos[0] + 0.25*region;
  } else {
    // printf("x-\n");
    c_centerPos[0] = p_centerPos[0] - 0.25*region;
  }

  if (quad_i / 2 == 0 ){
    // printf("y+\n");
    c_centerPos[1] = p_centerPos[1] + 0.25*region;
  } else {
    // printf("y-\n");
    c_centerPos[1] = p_centerPos[1] - 0.25*region;
  }

}

void calc_region(double vor_b[][2], int nvor_b){
  int i;
  double values[nvor_b];
  double region[2];

  for (i=0; i<nvor_b; i++){ values[i] = vor_b[i][0]; }
  region[0] = fabs(maxval(values,nvor_b) - minval(values,nvor_b)) + 0.0001;
  for (i=0; i<nvor_b; i++){ values[i] = vor_b[i][1]; }
  region[1] = fabs(maxval(values,nvor_b) - minval(values,nvor_b)) + 0.0001;
  if (region[0] > region[1]) {
    qtree.node[0].region = region[0];
  } else {
    qtree.node[0].region = region[1];
  }
  for (i=0; i<nvor_b; i++){ values[i] = vor_b[i][0]; }
  qtree.node[0].center_pos[0] = 0.5*(maxval(values,nvor_b) + minval(values,nvor_b));
  for (i=0; i<nvor_b; i++){ values[i] = vor_b[i][1]; }
  qtree.node[0].center_pos[1] = 0.5*(maxval(values,nvor_b) + minval(values,nvor_b));
}


void make_tree(int nvor_b, 
               double vor_b[nvor_b][2],
               double cir[nvor_b], 
               double cor[nvor_b])
{
  int i;
  int num_nodes; // number of nodes
  static int max_level = 1;
  double region[2];

  // printf("make tree\n");
  vors= (Vortex *)malloc(sizeof(Vortex) * nvor_b);
  for (i=0; i<nvor_b; i++){
    vors[i].pos[0] = vor_b[i][0];
    vors[i].pos[1] = vor_b[i][1];
    vors[i].cir = cir[i];
    vors[i].cor = cor[i]; 
  }
  // qtree.nterms = 0;
  // qtree.theta = 0;
  // qtree.r_min = 0; 
  // qtree.num_p = 250;

  read_fmm_param();

  // qtree.nterms = nterms;
  // qtree.theta = theta;
  // qtree.r_min = r_min;

  // searching max_level and make tree
  // stat:  
  //   max_level is sufficient:     qtree.stat transitions from 1 to 0 and end the loop.
  //   max_level is not sufficient: After qtree.stat goes from 1 to 0, it goes to -1 in split_cell, detects it and goes back to 1 again.
  qtree.stat = 1;
  while (qtree.stat == 1){ 
    // printf("max_level: %i\n",max_level);
    for (i=0; i<15; i++) { qtree.level_hist[i] = 0; }
    qtree.max_level = max_level;

    // initialize node
    for (i=0; i<sizeof(qtree.node)/sizeof(qtree.node[0]); i++){
      if(qtree.node[i].outer_coeffs){ free(qtree.node[i].outer_coeffs); }
    }
    if (qtree.node) { free(qtree.node); }
    num_nodes = (pow(4,qtree.max_level+1)-1)/3;
    qtree.num_nodes = num_nodes;
    qtree.node = (Node*)malloc(sizeof(Node) * num_nodes);
    for (i=0; i<num_nodes; i++){
      qtree.node[i].is_leaf = 0;
      qtree.node[i].num_vor = 0;
      qtree.node[i].level = getLevel(i);
    }

    // qtree.node[0].interactionset = 0
    // initialize root cell
    qtree.node[0].level = 0;
    qtree.node[0].num_vor = nvor_b;
    qtree.node[0].outer_coeffs = (double complex *)malloc(sizeof(double complex) * qtree.nterms);
    // qtree.node[0].inner_coeffs = (double complex *)malloc(sizeof(double complex) * qtree.nterms);
    for (i=0; i<qtree.nterms; i++){ qtree.node[0].outer_coeffs[i] = 0; }

    // set vortex on root cell
    qtree.node[0].vor_id = (int*)malloc(sizeof(int) * nvor_b);
    for (i=0; i<nvor_b; i++){ qtree.node[0].vor_id[i] = i; }

    calc_region(vor_b, nvor_b);


    // split cell (recurcive function)
    qtree.stat = 0;
    split_cell(0,1);
    // printf("-----------\n");
    // for (i=0;i<sizeof(qtree.level_hist)/sizeof(qtree.level_hist[0]);i++){
    //   printf("level:%i\n",qtree.level_hist[i]);
    // }

    // Judgment of End
    if (qtree.stat == -1) {
      qtree.stat = 1;
      max_level++;
    }
  }

  S2M();
  M2M();




  // quadrant_pos(pos,c_pos);

}

// void calc_centerPos(){
// }

void split_cell(int p_index, int flag){
  int p_num_vor;
  int c_index[4];
  int i,j;
  int quad_index,node_index;
  Vortex tmp_vor;

  if (qtree.node[p_index].level + 1 > qtree.max_level){
    qtree.stat = -1;
    return;
  }

  p_num_vor = qtree.node[p_index].num_vor;
  
  // initialize child cell 
  getChild(p_index,c_index);
  for (i=0;i<4;i++){ 
    qtree.node[c_index[i]].vor_id = (int *)malloc(sizeof(int) * p_num_vor);
    for (j=0; j<p_num_vor; j++){
      qtree.node[c_index[i]].vor_id[j] = -1;
    }
  }

  // put vortex to child cell
  // printf("---------\n");
  for (i=0; i<p_num_vor; i++){
    tmp_vor = vors[qtree.node[p_index].vor_id[i]];
    quad_index = quadrant_pos(tmp_vor.pos,qtree.node[p_index].center_pos);
    node_index = c_index[quad_index];
    qtree.node[node_index].num_vor = qtree.node[node_index].num_vor + 1;
    qtree.node[node_index].vor_id[qtree.node[node_index].num_vor-1] = qtree.node[p_index].vor_id[i];
  }

  for (i=0;i<4;i++){

    calc_CenterPos(i, 
                      qtree.node[p_index].center_pos,
                      qtree.node[p_index].region, 
                      qtree.node[c_index[i]].center_pos);
    // printf("%i:%i\n",c_index[i],qtree.node[c_index[i]].num_vor);
    // printf("center:%f %f\n",qtree.node[c_index[i]].center_pos[0],qtree.node[c_index[i]].center_pos[1]);
    qtree.node[c_index[i]].id = c_index[i];
    qtree.node[c_index[i]].level = qtree.node[p_index].level + 1;
    qtree.node[c_index[i]].region = 0.5*qtree.node[p_index].region;
    qtree.node[c_index[i]].outer_coeffs = (double complex *)malloc(sizeof(double complex)*qtree.nterms);
    for (j=0; j<qtree.nterms; j++) { qtree.node[c_index[i]].outer_coeffs[j] = 0.0; }

    if (qtree.node[c_index[i]].num_vor < qtree.num_p) {
      if (flag) {
        qtree.node[c_index[i]].is_leaf = 1;
        qtree.level_hist[qtree.node[c_index[i]].level] += 
              qtree.node[c_index[i]].num_vor;
      }
      if (qtree.node[c_index[i]].level < qtree.max_level){
        if (qtree.stat != -1) { split_cell(c_index[i], 0); }
      }
    } else {
      if (qtree.stat != -1) { split_cell(c_index[i], 1); }
    }
  }

}
