#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "quadtree.h"
#include "fmm.h"
#include "morton_order.h"

int combination(int n, int r)
{
    if (r == 0 || r == n)
        return (1);
    else if (r == 1)
        return (n);
    return (combination(n - 1, r - 1) + combination(n - 1, r));
}


void S2M(){
  int i,j,k;
  double complex z_j,z_star;
  for (i=0; i<qtree.num_nodes; i++){
    if (qtree.node[i].is_leaf) {
      z_star = qtree.node[i].center_pos[0]+qtree.node[i].center_pos[1] * I;
      for (k=0; k<qtree.nterms; k++){
        for (j=0; j<qtree.node[i].num_vor; j++){
          z_j = vors[qtree.node[i].vor_id[j]].pos[0]
                + vors[qtree.node[i].vor_id[j]].pos[1] * I;
          qtree.node[i].outer_coeffs[k] += vors[qtree.node[i].vor_id[j]].cir*cpow((z_j-z_star),k); 
        }
      }
    }
  }
}

void M2M(){
  int i,j,k,l,n;
  int p_index;
  double complex coeffs[qtree.nterms], p_pos,c_pos;
  int level;

  for (i=qtree.num_nodes-1; i>4; i-=4){
    p_index = getParent(i);
    p_pos = qtree.node[p_index].center_pos[0] + qtree.node[p_index].center_pos[1] * I;
    for (j=0;j>-4;j--){
      if (qtree.node[i+j].num_vor<1) {continue;}
      c_pos = qtree.node[i+j].center_pos[0] + qtree.node[i+j].center_pos[1] * I;
      for (n=0;n<qtree.nterms;n++){ coeffs[n] = 0.0; }
      for (n=0; n<qtree.nterms; n++){
        for (k=0; k<n+1; k++){
          coeffs[n] += qtree.node[i+j].outer_coeffs[k]*combination(n,k)*cpow(-(p_pos-c_pos),n-k);
        }
      }
      for (n=0;n<qtree.nterms;n++){
        qtree.node[p_index].outer_coeffs[n] += coeffs[n];
      }
    }
  }
}

void M2S(double complex z_i, double complex *phi, int calculation_node){
  double complex z_star;
  int k;
  z_star = qtree.node[calculation_node].center_pos[0]+qtree.node[calculation_node].center_pos[1]*I;
  for (k=0;k<qtree.nterms-1;k++){
    *phi += cpow((z_i - z_star),(-k-1)) * qtree.node[calculation_node].outer_coeffs[k];
  }
}

void biot_tree(double x,double y, double *u,double *v){
  double complex phi;
  double complex z_i,z_star;
  double u_direct[2];
  int i,k;
  int c_index[4];
  double r;
  int *fifo_index;
  int index;
  int t_index;

  phi = 0.0;
  u_direct[0] = 0.0; u_direct[1] = 0.0;

  fifo_index = malloc(sizeof(int)*qtree.num_nodes);
  fifo_index[0] = 0;
  for (i=1; i<qtree.num_nodes; i++){
    fifo_index[i] = -1;
  } 

  index = 0;
  while(1){

    t_index = fifo_index[index];
    fifo_index[index] = -1;
    index--;

    z_i = x + y*I;
    if (qtree.node[t_index].is_leaf){
      direct_evaluate(x,y,t_index,u_direct);
    } else {
      getChild(t_index,c_index);
      for (i=0;i<4;i++){
        r = sqrt(pow(x - qtree.node[c_index[i]].center_pos[0],2) + pow(y - qtree.node[c_index[i]].center_pos[1],2));
        if (qtree.node[c_index[i]].region > qtree.theta * r){           
          index++;
          fifo_index[index] = c_index[i];
        } else {
          M2S(z_i, &phi, c_index[i]);
        }
      }
    }  
    if (index == -1){ break; }
  }
  free(fifo_index);

  phi = 0.5*(-cimag(phi)+creal(phi)*I)/M_PI;
  *u =-creal(phi) + u_direct[0];
  *v = cimag(phi) + u_direct[1];
}


void direct_evaluate(double x, double y,  int t_index, double u_direct[2]){
  int i;
  double rx,ry,r,rv2,xai,rv,ax,ay;
  int vor_id;

  ax = 0.0; ay = 0.0;
  for (i=0;i<qtree.node[t_index].num_vor;i++){
    vor_id = qtree.node[t_index].vor_id[i];
    rx = x - vors[vor_id].pos[0];
    ry = y - vors[vor_id].pos[1];
    r = sqrt(pow(rx,2)+pow(ry,2));
    if ( r > 1.0e-10) {
      rv2 = 1.0 / pow(r,2);
      xai = r / vors[vor_id].cor;
      rv = vors[vor_id].cir * (1.0 - exp(-xai*xai));
      ax -= rv*ry*rv2;
      ay += rv*rx*rv2;
    }
  }
  u_direct[0] += 0.5*ax/M_PI;
  u_direct[1] += 0.5*ay/M_PI;
}

