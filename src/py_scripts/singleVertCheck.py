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
  sm=0
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
        #n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
        n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
        #print '(',i,j,k,')--', n
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
  
  fw = open ('asciiFile.txt', 'w')
  print 'tempGrid', numVertPerDirection
  for i in range (numVertPerDirection):
    for j in range (numVertPerDirection):
      for k in range (numVertPerDirection):
        #print >>fw, '(',i,j,k,')',tempGrid[i][j][k],
        print >>fw,tempGrid[i][j][k],
      print >>fw,' '
    print >>fw,' '
  
      
  for i in range (numVertPerDirection):
    for j in range (numVertPerDirection):
      for k in range (numVertPerDirection):
        sm=0
        if i < (numVertPerDirection -1) and j < (numVertPerDirection -1) and k < (numVertPerDirection -1):
					for m in range (2):
					  for n in range (2):
					    for p in range (2):
					      sm = sm + tempGrid[i+m][j+n][k+p]
					      #print >>fw1,'(',i+m,j+n,k+p,')=',tempGrid[i+m][j+n][k+p],
					      
					ind = i*numVertPerDirection*numVertPerDirection+j*numVertPerDirection+k
					#print ' --',i,j,k,'ind',ind,'--',sm
					tempgv = gridVertType()
					tempgv = grid[ind]
					tempgv.sm = sm
					if sm >=1 and sm <=7:
						tempgv.isIntersect = True
          
  return grid
    
'''
compute the centroid of the object supplied as input 
parameters
input the vertex locations and the number of the vertices
'''        
def computeCentroid(vertLoc, numv):
  sm = [0, 0, 0]
  for i in range (numv):
    sm[0] = sm[0] + vertLoc[3*i]
    sm[1] = sm[1] + vertLoc[3*i+1]
    sm[2] = sm[2] + vertLoc[3*i+2]    
  centroid =[sm[x] / (numv) for x in range(3)]
  return centroid
  
def main():
  if  len(sys.argv) != 9:
    print ' parameters 1> .OFF file \n 2> height\n 3>scale 4> text file with nrrd data \n 5>origin 6 > longestDim'
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
    print '***********************'
    print 'scale', scale
    originCompare = [int(sys.argv[5+x]) for x in range (3)]
    print 'originCompare', originCompare
    
    L = float (sys.argv[8])
    print 'longest dim ', L
    print  'Dimesion (', L/scale ,')'
    print 'new origin', [originCompare[x]-( (L*0.5/scale))  for x in range (3)]
    NewOrigin = [originCompare[x]-( (L*0.5/scale))  for x in range (3)]
    centroid = computeCentroid(vertLoc, numv)
    halfl = (L/scale)/2.0
    ul = (L /scale)/(2**height)
    print 'unit length ' , ul
    numVertPerDirection = 2**height+1
    print "numVertPerDirection ", numVertPerDirection 
    
    print 'computed centroid ', centroid
    print '***********************'
    for i in range (numv):
      baseLoc = [ int ( (vertLoc[3*i+x] - ( centroid[x] - halfl)) /ul ) for x in  range (3)]
      #debug
      print baseLoc, 
      ind = baseLoc[0]*numVertPerDirection*numVertPerDirection + baseLoc[1]*numVertPerDirection+ baseLoc[2]
      print ind, 
      tempgv = gridVertType()
      tempgv = grid[ind]
      tempgv.gridLoc=baseLoc
      print tempgv.gridLoc
      tempgv.numIntersection = tempgv.numIntersection + 1 
      tempgv.ptsLocs.append([vertLoc[3*i+0] , vertLoc[3*i+1], vertLoc[3*i+2]])
    
    #read Nrrd data and set isIntersect
    NrrdsetGrid = setNrrdData(grid, height)
    
    print 'Should not have  a  point but does '  
    for i in range (numVertPerDirection):
	  for j in range (numVertPerDirection):
		  for k in range (numVertPerDirection):
			  n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
			  tempgv = gridVertType()
			  tempgv = grid[n]
			  if  tempgv.numIntersection > 0  and tempgv.isIntersect == False:
			    print    'case 1',tempgv.x*ul ,tempgv.y*ul,tempgv.z*ul
    
    print 'Should have  a  point but does not '	     
    for i in range (numVertPerDirection):
	  for j in range (numVertPerDirection):
		  for k in range (numVertPerDirection):
			  n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
			  tempgv = gridVertType()
			  tempgv = grid[n]
			  if  tempgv.numIntersection == 0  and tempgv.isIntersect == True:
			    print    'case 2', tempgv.x*ul , tempgv.y*ul, tempgv.z*ul 
    
    print 'OUT'
    fw2 = open ('out1.off','w')
    for i in range (numVertPerDirection):
			for j in range (numVertPerDirection):
				for k in range (numVertPerDirection):
					n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
					tempgv = gridVertType()
					tempgv = grid[n]
					if  tempgv.numIntersection == 1 : 
						#print  tempgv.gridLoc
						#print   tempgv.gridLoc[0],' ',tempgv.gridLoc[1],' ',tempgv.gridLoc[2]
						#print >>fw2, " ".join(map(str, tempgv.ptsLocs[0])), ' grid location ',
						#print >>fw2, tempgv.gridLoc
						print >> fw2 , ' '.join(map(str,[tempgv.ptsLocs[0][x] - ( centroid[x] - halfl) for x in range (3)]))
    
    print 'OUT2'
    fw3 = open ('out2.off','w')
    for i in range (numVertPerDirection):
			for j in range (numVertPerDirection):
				for k in range (numVertPerDirection):
					n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
					tempgv = gridVertType()
					tempgv = grid[n]
					if  tempgv.numIntersection > 0:
					  print >>fw3,'NUM POINTS ', tempgv.numIntersection
					  print >>fw3, 'grid loc ',n,' -- ',[tempgv.gridLoc[x]  for x in range (3)],'--',[tempgv.gridLoc[x]*ul  for x in range (3)],' \n '
					  for m in  range (len(tempgv.ptsLocs)): 
						  #print >>fw3," ".join(map(str, tempgv.ptsLocs[m])),' grid location ',
						   
						  print >> fw3 , ' '.join(map(str,[tempgv.ptsLocs[m][x] - ( centroid[x] - halfl) for x in range (3)])),
						  print >> fw3 , '[',' '.join(map(str,[tempgv.ptsLocs[m][x] for x in range (3)])),']',
						  print >> fw3 ,'diff', [(tempgv.ptsLocs[m][x] - ( centroid[x] - halfl))-tempgv.gridLoc[x] *ul  for x in range (3)]
						  diff = [(tempgv.ptsLocs[m][x] - ( centroid[x] - halfl))-tempgv.gridLoc[x] *ul  for x in range (3)]
						  if abs(diff[0] - 0.0)<0.001 or  abs(diff[1] - 0.0)<0.001 or abs(diff[2] - 0.0)<0.001:
						    print >>fw3 , 'FALSE POSITIVE'
					  print >>fw3, ' '

    print 'OUT3'#####  all points with more than 1 points
    fw4 = open ('outMore.off','w')
    for i in range (numVertPerDirection):
			for j in range (numVertPerDirection):
				for k in range (numVertPerDirection):
					n = (numVertPerDirection)*(numVertPerDirection)*i + (numVertPerDirection)*j + k
					tempgv = gridVertType()
					tempgv = grid[n]
					if  tempgv.numIntersection > 1 :
					  for m in  range (len(tempgv.ptsLocs)): 
						  print >> fw4 , ' '.join(map(str,[tempgv.ptsLocs[m][x]- ( centroid[x] - halfl) for x in range (3)]))						  
						  #print >> fw3, tempgv.gridLoc
if __name__ == "__main__":
    main()
    
