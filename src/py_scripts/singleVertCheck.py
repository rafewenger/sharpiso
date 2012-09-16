#!/usr/bin/python
import sys 
import numpy



class gridVertType:
  ind = 0
  x=0
  y=0
  z=0
  isIntersect = False
  numIntersection = 0
  gridLoc=[]
  ptsLocs=[]

def readOffFile ():
  print 'reading the off file'
  fileName = sys.argv[1]
  print 'fname', fileName
  f = open (fileName, 'r')
  return f.readlines()

  
def computeVertLoc (data):
  
  numv=int(data[1].split(' ')[0])
  locList =[]
  for i in range (2,2+numv):
    loc = data[i].split(' ')
    locList.append(float(loc[0]))
    locList.append(float(loc[1]))
    locList.append(float(loc[2]))
  return [locList, numv]
  

def buildGrid(vertLoc, height):
  '''
  build the grid, which should be {2^height}^3. DEBUG 
  each point has a isIntersect which is gonna be fixed by the nrrd data
  the numPoints and there locations are set by the ...
  '''
  print 'height', height
  grid=[]
  numVertPerDirection = 2**height+1 ## DEBUG 
  for i in range (numVertPerDirection):
    for j in range (numVertPerDirection):
      for k in range (numVertPerDirection):
        n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
        gv = gridVertType()
        gv.ind = n
        gv.isIntersect = False
        gv.numIntersection = 0 
        gv.ptsLocs=[]
        gv.gridLoc=[]
        gv.x=i
        gv.y=j
        gv.z=k
        grid.append(gv)
  return grid  



def setNrrdData(grid, height):
  fi = open (sys.argv[4],'r')
  nrrddata=fi.read().split(' ')
  numVertPerDirection = 2**height+1
  '''
  tempgrid store the scalar values from the nrrd
  must be of size , 2^d+1
  '''
  tempGrid=numpy.zeros((numVertPerDirection, numVertPerDirection, numVertPerDirection))
  indx = 0
  for i in range (numVertPerDirection):
    for j in range (numVertPerDirection):
      for k in range (numVertPerDirection):        
        tempGrid[i][j][k]=float(nrrddata[indx])
        indx=indx+1
  '''
  print 'tempGrid'
  for i in range (numVertPerDirection):
    for j in range (numVertPerDirection):
      for k in range (numVertPerDirection):
        print tempGrid[i][j][k],
      print ' '
    print ' '
  return [0,0,0]
  '''        
  for i in range (numVertPerDirection):
    for j in range (numVertPerDirection):
      for k in range (numVertPerDirection):
        sum=0
        if i < (numVertPerDirection -1) and j < (numVertPerDirection -1) and k < (numVertPerDirection -1):
					sum = tempGrid[i][j][k]+tempGrid[i][j][k+1]+\
					tempGrid[i][j+1][k]+tempGrid[i][j+1][k+1]+tempGrid[i+1][j][k]+\
					tempGrid[i+1][j][k+1]+tempGrid[i+1][j+1][k]+tempGrid[i+1][j+1][k+1]
					if sum >=1 and sum <=7:
						ind = i*numVertPerDirection*numVertPerDirection+j*numVertPerDirection+k
						tempgv = gridVertType()
						tempgv = grid[ind]
						tempgv.isIntersect = True
          
  return grid
        
  
def main():
  if  len(sys.argv) != 8:
    print ' parameters 1> .OFF file \n 2> height\n 3>scale 4> text file with nrrd data \n 5>origin'
    print sys.argv, 'length',len(sys.argv)
  else :
    data=[]
    #read the off file 
    data  = readOffFile ()
    #compute the vertLoc
    vertLocAndNumv = computeVertLoc(data)
    
    vertLoc = vertLocAndNumv[0]
    numv = vertLocAndNumv[1]
    #build grid
    height = int (sys.argv[2])
    
    grid = buildGrid(vertLoc, height)
    
    #unit length
    scale = float(sys.argv[3])
    
    print 'scale', scale
    origin = [int(sys.argv[5+x]) for x in range (3)]
    print 'origin',origin
    ul = (1.0/scale)/(2**height)
    print ' unit length ' , ul
    numVertPerDirection = 2**height+1
    for i in range (numv):
      baseLoc = [ int ( (vertLoc[3*i+x] - origin[x]+0.5) /ul ) for x in  range (3)]
      ind = baseLoc[0]*numVertPerDirection*numVertPerDirection + baseLoc[1]*numVertPerDirection+ baseLoc[2]
      tempgv = gridVertType()
      tempgv = grid[ind]
      tempgv.gridLoc=baseLoc
      tempgv.numIntersection = tempgv.numIntersection + 1 
      tempgv.ptsLocs.append([vertLoc[3*i+0] - origin[0]+0.5, vertLoc[3*i+1] - origin[1]+0.5, vertLoc[3*i+2] - origin[2]+0.5])
    
    #read Nrrd data and set isIntersect
    NrrdsetGrid = setNrrdData(grid, height)
    '''    
    for i in  range (len(grid)):
      tempgv = gridVertType()
      tempgv = grid[i]
      if tempgv.numIntersection > 1:
        print  tempgv.x,',',tempgv.y,',',tempgv.z,' base loc  [',tempgv.gridLoc,'] -- numvert ', tempgv.numIntersection
        print 'loc ',tempgv.ptsLocs
        
    '''    
    '''
    for i in range (numVertPerDirection):
			for j in range (numVertPerDirection):
				for k in range (numVertPerDirection):
					n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
					tempgv = gridVertType()
					tempgv = grid[n]
					if tempgv.isIntersect == True:
						print  tempgv.x,tempgv.y,tempgv.z   
    '''
    print 'OUT'
    for i in range (numVertPerDirection):
			for j in range (numVertPerDirection):
				for k in range (numVertPerDirection):
					n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
					tempgv = gridVertType()
					tempgv = grid[n]
					if tempgv.numIntersection  < 1 and tempgv.isIntersect == True  :
						print  tempgv.x,tempgv.y,tempgv.z     
if __name__ == "__main__":
    main()
    
