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
    int level, value, i;
    level = 0;
    value = index;
    while(1){
      value -= (int)pow(4,level);
      if(value>=0){
        level++;
      } else {
        break;
      }
    }
    return level;
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



