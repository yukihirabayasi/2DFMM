import matplotlib.pyplot as plt

def local2global():
  pass

def global2local(g_index,level):
  l_index = g_index - (4**level-1)/3
  return int(l_index)

def num_particles_cell(index,level,vorPos,region):
  # How many particles the cell has 
  # index start from 1, not 0 and linear 

  cell_size = region/2**level
  index = global2local(index,level) - 1
  index = format(index, 'b')
  ix = int(index[::-2][::-1],2)
  iy = int(index[-2::-2][::-1],2)
  #ix = 1
  #iy = 0
  x = (ix+0.5)*cell_size - 0.5*region 
  y = -(iy+0.5)*cell_size + 0.5*region 
  print("pos",x,y)

  counter = 0
  for i in range(0,len(vorPos)):
    if (x-0.5*cell_size <= vorPos[i][0] and vorPos[i][0] <= x+0.5*cell_size and
        y-0.5*cell_size <= vorPos[i][1] and vorPos[i][1] <= y+0.5*cell_size):
      counter += 1

  print("particles in the cell:",counter)

def lattice_plot(pos,minPos,maxPos,level):
  plt.xlim([-0.5*region,0.5*region])
  plt.ylim([-0.5*region,0.5*region])
  pos = list(zip(*pos))
  plt.scatter(pos[0],pos[1],marker=".",s=0.03,color='black',alpha=1.0)
  plt.gca().set_aspect('equal', adjustable='box')
  
  plt.grid()
  grid = []
  for i in range(-2**(level-1),2**(level-1)):
    grid.append(i*region/(2**level)) 
  plt.xticks(grid )
  plt.yticks(grid )
  plt.tick_params(labelbottom=False,
                  labelleft=False,
                  labelright=False,
                  labeltop=False)
  
  plt.show()


if __name__ == '__main__':
  pos = []
  f = open('./vor.dat', 'r')
  for line in f:
      pos.append(line.replace("\n","").replace(" ","").split(","))
  f.close()
  
  maxPos = [0,0]
  minPos = [10,10]
  
  for i in range(0,len(pos)):
    pos[i][0] = float(pos[i][0])
    pos[i][1] = float(pos[i][1])
  
  for i in range(0,len(pos)):
    if (maxPos[0] < pos[i][0]):
      maxPos[0] = pos[i][0]
    if (minPos[0] > pos[i][0]):
      minPos[0] = pos[i][0]
  
    if (maxPos[1] < pos[i][1]):
      maxPos[1] = pos[i][1]
    if (minPos[1] > pos[i][1]):
      minPos[1] = pos[i][1]

  centerPos = [0,0]
  #centerPos[0] = 0.5*(maxPos[0]+minPos[0])
  #centerPos[1] = 0.5*(maxPos[1]+minPos[1])
  
  if (abs(maxPos[0]-minPos[0]) < abs(maxPos[1]-minPos[1])):
    region = abs(maxPos[1]-minPos[1])
  else:
    region = abs(maxPos[0]-minPos[0])
  
  for i in range(0,len(pos)):
    pos[i][0] = pos[i][0] - centerPos[0]
    pos[i][1] = pos[i][1] - centerPos[1]
  
  index = 18 
  level = 10
  num_particles_cell(index,level,pos,region)
  lattice_plot(pos,minPos,maxPos,level)

