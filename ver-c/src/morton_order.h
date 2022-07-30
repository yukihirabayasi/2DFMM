#include <math.h>

extern int global2local(int globalIndex, int level);

extern int local2global(int localIndex, int level);

extern int getLevel(int index);

extern void getChild(int p_index, int c_index[4]);

extern int getParent(int c_index);

