#include <complex.h>
#include <stdlib.h>

typedef struct {
  // unsigned long interaction_set[27];
  int level;
  // int id;
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

int main(){
    QuadTree qtree;
    qtree.node = (Node*)malloc(sizeof(Node) * 4);
    qtree.node[0].outer_coeffs = (double complex*)malloc(sizeof(double complex) * qtree.nterms);
    qtree.node[0].inner_coeffs = (double complex*)malloc(sizeof(double complex) * qtree.nterms);
}