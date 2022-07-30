#include <math.h>

extern int global2local(int globalIndex, int level);

extern int local2global(int localIndex, int level);

extern int getLevel(int index);

extern void getChild(int p_index, int c_index[4]);

extern int getParent(int c_index);

extern unsigned long BitSeparate32( unsigned long n );

extern unsigned long CompactBy1(unsigned long x);

extern unsigned long pos2index( unsigned long x, unsigned long y );

extern void index2pos( unsigned long index, unsigned long pos[2] );

