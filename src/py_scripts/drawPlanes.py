#!/usr/bin/python 
import sys
import numpy
base=[]
fw = open('test.quad', 'w')

def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],\
         a[2]*b[0] - a[0]*b[2],\
         a[0]*b[1] - a[1]*b[0]]

    return c
    
def fileWrite (endPt):
  print 'endPt ', endPt
  for i in endPt:
    print >>fw,i,
  print >>fw,'\n',
  
  
  
        
def computePlane(vertexPoint, vertexGrad, planePoint, Mag):
  print '**compute Plane**'
  #print vertexPoint,
  #print vertexGrad,
  #print planePoint,
  Mag = float(Mag)
  '''
  P1=[]
  p = (float(vertexGrad.split(',')[0]))*float(vertexPoint.split(',')[0])*Mag + \
  float(vertexGrad.split(',')[1])*float(vertexPoint.split(',')[1])*Mag + \
  float(vertexGrad.split(',')[2])*float(vertexPoint.split(',')[2])*Mag 
  
  
  P1.append(float(vertexPoint.split(',')[0])) # p1x
  P1.append(float(vertexPoint.split(',')[1])) # p1y
  P1.append( (p-float(vertexGrad.split(',')[0])*Mag*float(vertexPoint.split(',')[0]) - \
  float(vertexGrad.split(',')[1])*Mag*float(vertexPoint.split(',')[1]))/(float(vertexGrad.split(',')[2])*Mag))
  
  P2=[]
  P2.append(float(vertexPoint.split(',')[0])) # p2x
  P2.append( (p-float(vertexGrad.split(',')[0])*Mag*float(vertexPoint.split(',')[0]) - \
  float(vertexGrad.split(',')[2])*Mag*float(vertexPoint.split(',')[2]))/(float(vertexGrad.split(',')[1])*Mag))
  P2.append(float(vertexPoint.split(',')[2])) # p2z

  P3=[]
  P3.append(float(planePoint.split(',')[0]))
  P3.append(float(planePoint.split(',')[1]))
  P3.append(float(planePoint.split(',')[2]))  
  '''
  print 'magnitude', Mag
  vertexPt = vertexPoint.split(',')
  vertexPt = [float(x) for x in vertexPt ]
  
  vertexGrd = vertexGrad.split(',')
  vertexGrd = [float(x) for x in vertexGrd]
  
  planePt = planePoint.split(',')
  planePt = [float(x) for x in planePt ]
  
  print vertexPt,
  print vertexGrd,
  print planePt
  index_min=vertexGrd.index(min([abs(x) for x in vertexGrd]))
  print 'index min ',index_min
  
  
  baseDirc=[0.0,0.0,0.0]
  
  for x in range (3):
    if x == index_min:
      baseDirc[x]=1.0
  
  print 'updated base direction ' , baseDirc
  u=cross (baseDirc,vertexGrd)
  v= cross (vertexGrd, u)
  print 'U ', u, ' V ',v
  maxDistU=1
  maxDistV=1
  endPt1=[planePt[0]+maxDistU*u[0]+maxDistV*v[0],planePt[1]+maxDistU*u[1]+maxDistV*v[1], planePt[2]+maxDistU*u[2]+maxDistV*v[2]]
  maxDistU=1
  maxDistV=-1
  endPt2=[planePt[0]+maxDistU*u[0]+maxDistV*v[0],planePt[1]+maxDistU*u[1]+maxDistV*v[1], planePt[2]+maxDistU*u[2]+maxDistV*v[2]]

  maxDistU=-1
  maxDistV=-1
  endPt3=[planePt[0]+maxDistU*u[0]+maxDistV*v[0],planePt[1]+maxDistU*u[1]+maxDistV*v[1], planePt[2]+maxDistU*u[2]+maxDistV*v[2]]

  maxDistU=-1
  maxDistV=1
  endPt4=[planePt[0]+maxDistU*u[0]+maxDistV*v[0],planePt[1]+maxDistU*u[1]+maxDistV*v[1], planePt[2]+maxDistU*u[2]+maxDistV*v[2]]
  
  fileWrite (endPt1)
  fileWrite (endPt2)
  fileWrite (endPt3)
  fileWrite (endPt4)
 

  
  
  '''
  p = Mag*vertexGrd[0]*planePt[0] + Mag*vertexGrd[1]*planePt[1] + Mag*vertexGrd[2]*planePt[2]
  print 'p',p
  
  endPt1=[base[0]-1, base[1]-1, base[2]-1]
  endPt3=[base[0]+2, base[1]+2, base[2]+2]
  
  endPt2=[base[0]-1, base[1]-1, base[2]+2]
  endPt4=[base[0]+2, base[1]+2, base[2]-1]
  print endPt1
  print endPt2
  print endPt3
  10 14 print endPt4
  
  print >> fw, 'test'
  print >> fw , endPt1[0], endPt1[1],
  pt = (p-(Mag*vertexGrd[0]*endPt1[0]+Mag*vertexGrd[1]*endPt1[1]))/(Mag*vertexGrd[2])
  print >>fw , pt

  print >> fw , endPt2[0], endPt2[1],
  pt = (p-(Mag*vertexGrd[0]*endPt2[0]+Mag*vertexGrd[1]*endPt2[1]))/(Mag*vertexGrd[2])
  print >>fw , pt
  
  print >> fw , endPt3[0], endPt3[1],
  pt = (p-(Mag*vertexGrd[0]*endPt3[0]+Mag*vertexGrd[1]*endPt3[1]))/(Mag*vertexGrd[2])
  print >>fw , pt  
     
  print >> fw , endPt4[0], endPt4[1],
  pt = (p-(Mag*vertexGrd[0]*endPt4[0]+Mag*vertexGrd[1]*endPt4[1]))/(Mag*vertexGrd[2])
  print >>fw , pt
  '''

  
  
  
def readFile():
  print '**read the file info**'
  f = open('test.txt', 'r')
  for line in f:
        lineTrns = line.translate(None, "()?!/;:")
        largeInfo=lineTrns.split(' ')
        #print largeInfo
        flag_call=False
        for idx,data in enumerate(largeInfo):
          if data=='Point':
            vertexPoint = largeInfo[idx+3]
            flag_call=True
          if data=='grad':
            vertexGrad = largeInfo[idx+1]
          if data=='point':
            planePoint = largeInfo[idx+1]
          if data=='Mag':
            Mag = largeInfo[idx+1]
        if flag_call:
          print "[",vertexPoint,
          print vertexGrad,
          print planePoint,"]"
          computePlane(vertexPoint, vertexGrad, planePoint, Mag)
        
def main():
  print 'read a set of points and draw the planes'
  readFile()
if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.' 
    print 'Argument List:', str(sys.argv)
    if len(sys.argv)==4:
      print 'base cube position is (', float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),")" 
      base.append(float(sys.argv[1]))
      base.append(float(sys.argv[2]))
      base.append(float(sys.argv[3]))
      print >> fw, "QUAD"
      main()
    else:
      print 'Not enough arguments'
    
