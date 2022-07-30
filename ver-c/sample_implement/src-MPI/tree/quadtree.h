#include <complex.h>
#include <stdlib.h>

typedef struct {
  // unsigned long interaction_set[27];
  int level;
  int id;
  int *vor_id;
  int num_vor;
  double center_pos[2];
  double region;
  double complex *outer_coeffs;
  // double complex *inner_coeffs;
  int is_leaf;
} Node;

typedef struct {
  Node *node;
  int max_level;
  int nterms;
  double theta;
  int level_hist[15];
  double r_min;
  int num_p;
  int stat; // for maxLevel loop
  int num_nodes;
} QuadTree;

typedef struct {
  double pos[2];
  double cir;
  double cor;
} Vortex;


int quadrant_pos(double pos[2], double c_pos[2]);

void make_tree(int nvor_b, 
               double vor_b[2][nvor_b], 
               double cir_b[nvor_b], 
               double cor_b[nvor_b]);

void add_vortex();

void calc_CenterPos(int quad_i, double *p_centerPos, double region, double *c_centerPos);

void split_cell(int p_index, int flag);

extern QuadTree qtree;
extern Vortex *vors;