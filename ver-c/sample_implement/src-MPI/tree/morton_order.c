#include <math.h>
#include "morton_order.h"

int global2local(int globalIndex, int level){
  int localIndex;
  localIndex = globalIndex - (pow(4,level)-1)/3;
  return localIndex;
}

int local2global(int localIndex, int level){
  int globalIndex;
  globalIndex = localIndex + (pow(4,level)-1)/3;
  return globalIndex;
}

int getLevel(int index){
  if (index == 0) {
    return 0;
  } else {
    return (int)log(3.0*(double)index)/log(4.0);
  }
}

void getChild(int p_index, int c_index[4]){
  int level;
  int index;
  int i;
  level = getLevel(p_index);
  index = global2local(p_index, level);
  index = index << 2;
  for (i=0; i<4; i++){
    c_index[i] = local2global(index + i,level+1);
  }
}

int getParent(int c_index){
  int p_index;
  int level;
  level = getLevel(c_index);
  c_index = global2local(c_index,level);
  p_index = c_index >> 2;
  p_index = local2global(p_index,level-1);
  return p_index;
}


unsigned long BitSeparate32( unsigned long n )
{
   n = (n|(n<<8)) & 0x00ff00ff;
   n = (n|(n<<4)) & 0x0f0f0f0f;
   n = (n|(n<<2)) & 0x33333333;
   return (n|(n<<1)) & 0x55555555;
}

unsigned long CompactBy1(unsigned long x)
{
    x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

unsigned long pos2index( unsigned long x, unsigned long y )
{
   return (BitSeparate32(x) | (BitSeparate32(y)<<1));
}

void index2pos(unsigned long index, unsigned long pos[2])
{            
  pos[0] = CompactBy1(index & 0x0f0f0f0f);
  pos[1] = CompactBy1((index & 0xf0f0f0f0) >> 1);
}

