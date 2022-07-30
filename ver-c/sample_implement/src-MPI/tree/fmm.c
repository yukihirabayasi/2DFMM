#include <stdio.h>
// #include <complex.h>
#include <math.h>
#include "quadtree.h"
#include "fmm.h"
#include "morton_order.h"

#define pi atan(0)

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
  double  complex z_j,z_star;
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

void biot_tree(double x,double y, double *u,double *v){
  double complex phi;
  double u_direct[2];
  phi = 0.0;
  u_direct[0] = 0.0; u_direct[1] = 0.0;
  tree_evaluate(x,y,0,&phi,u_direct);
  phi = 0.5*(-cimag(phi)+creal(phi)*I)/pi;
  *u =-creal(phi) + u_direct[0];
  *v = cimag(phi) + u_direct[1];
  // printf("%i",qtree.num_nodes);
}

void tree_evaluate(double x,double y, int p_index, double complex *phi, double u_direct[2] ){
  double complex z_i,z_star;
  double r;
  int c_index[4];
  int i,k;
  z_star = 0.0;


  z_i = x + y*I;
  if (qtree.node[p_index].is_leaf) {
    direct_evaluate(x,y,p_index,u_direct);
  } else {
    getChild(p_index,c_index);
    for (i=0;i<4;i++){
      r = sqrt(pow(x - qtree.node[c_index[i]].center_pos[0],2) + pow(y - qtree.node[c_index[i]].center_pos[1],2));
      if (qtree.r_min > r){
        direct_evaluate(x,y,c_index[i],u_direct);
      } else {
        if (qtree.node[c_index[i]].region > qtree.theta*r) {
          tree_evaluate(x,y,c_index[i],phi,u_direct);
        } else {
          z_star = qtree.node[c_index[i]].center_pos[0]+qtree.node[c_index[i]].center_pos[1]*I;
          for (k=0;k<qtree.nterms-1;k++){
            *phi += cpow((z_i - z_star),(-k-1)) * qtree.node[c_index[i]].outer_coeffs[k];
          }
        }
      }
    }
  }
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
  u_direct[0] += 0.5*ax/pi;
  u_direct[1] += 0.5*ay/pi;
}
