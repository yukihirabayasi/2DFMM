#include <complex.h>

void S2M();

void M2M();

void biot_tree(double x, double y, double *u, double *v);

void tree_evaluate(double x, double y, int p_index, double complex *phi, double u_direct[2]);

void direct_evaluate(double x, double y,  int t_index, double u_direct[2]);
