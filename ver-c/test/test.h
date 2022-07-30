#include <complex.h>

typedef struct {
  unsigned long interaction_set[27];
  int level;
  int id;
  int *vor_id;
  int num_vor;
  double center_pos[2];
  double region;
  double complex *outer_coeffs;
  double complex *inner_coeffs;
  int is_leaf;
} Node;

typedef struct {
  Node *node;
  int max_level;
  int nterms;
  double theta;
  int level_hist[15];
  double r_min;
  int stat; // for maxLevel loop
} QuadTree;


void make_tree(int *nvor_b, double poi_x[*nvor_b], double poi_y[*nvor_b], double cir[*nvor_b], double cor[*nvor_b]);

void biot_tree(double *x,double *y,double *u,double *v);

QuadTree qtree;

