#!/usr/bin/python 
import sys
import numpy
base=[]
fw = open('test.quad', 'w')
fwList = open ('geom-test-2.list','w')
fwLine = open ('cube.vect','w')
FirstPoint = [1]

def drawCube(vertexPt):
  print >>fwLine,'VECT\n3 6 3'
  for i in range (3):
    print >>fwLine,'2',
  print >>fwLine,''
  for i in range (3):
    print >>fwLine,'1',
  print >> fwLine,'\n',vertexPt[0],vertexPt[1],vertexPt[2], vertexPt[0]+1,vertexPt[1],vertexPt[2]
  print >> fwLine,vertexPt[0],vertexPt[1],vertexPt[2], vertexPt[0],vertexPt[1]+1,vertexPt[2]
  print >> fwLine,vertexPt[0],vertexPt[1],vertexPt[2], vertexPt[0],vertexPt[1],vertexPt[2]+1
  print >> fwLine,'1.0 0.0 0.0 1.0'
  print >> fwLine,'0.0 1.0 0.0 1.0'
  print >> fwLine,'0.0 0.0 1.0 1.0'
  
    
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],\
         a[2]*b[0] - a[0]*b[2],\
         a[0]*b[1] - a[1]*b[0]]

    return c

# write a endPoint to the quad file    
def fileWrite (endPt):
  print 'endPt ', endPt
  for i in endPt:
    print >>fw,i,
  print >>fw,'\n',
  
# write a shere to the list file 
def writeSphere (r,c,color):
  print >>fwList,"{appearance {  material {diffuse ",
  print >>fwList, color[0]," ",color[1]," ",color[2],"}}"
  print >>fwList,"SPHERE"
  print >>fwList,r
  print >>fwList,c[0],c[1],c[2]," }"
  
  
  
        
def computePlane(vertexPoint, vertexGrad, planePoint, Mag):
  print '**compute Plane**'
  #print vertexPoint,
  #print vertexGrad,
  #print planePoint,
  Mag = float(Mag)
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
  absVertexGrd = [abs(x) for x in vertexGrd]
  index_min=absVertexGrd.index(min(absVertexGrd))
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
  writeSphere(0.1,planePt,[0,1,0])
  writeSphere(0.1,vertexPt,[0,0,1])
  print vertexPt
  fileWrite (endPt1)
  fileWrite (endPt2)
  fileWrite (endPt3)
  fileWrite (endPt4)
  
  
def readFile():
  print '**read the file info**'
  f = open('test.txt', 'r')
  for line in f:
        lineTrns = line.translate(None, "()?!/;:")
        largeInfo=lineTrns.split(" ")
        print largeInfo
        flag_call=False
        for idx,data in enumerate(largeInfo):
          if data=='Point' and largeInfo[idx+1]=='':
            vertexPoint = largeInfo[idx+3]
            flag_call=True
          if data=='Point' and largeInfo[idx+1]!='':
            vertexPoint = largeInfo[idx+2]
            flag_call=True
          if data=='grad':
            vertexGrad = largeInfo[idx+1]
          if data=='point':
            planePoint = largeInfo[idx+1]
          if data=='Mag':
            Mag = largeInfo[idx+1]
        if flag_call:
          print "[vertexPoint ",vertexPoint,
          print "vertexGrad ",vertexGrad,
          print "planePoint ",planePoint,"]"
          computePlane(vertexPoint, vertexGrad, planePoint, Mag)
        
def main():
  print 'read a set of points and draw the planes'
  readFile()
if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.' 
    print 'Argument List:', str(sys.argv)
    if len(sys.argv)==7:
      print 'base cube position is (', float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),")" 
      base.append(float(sys.argv[1]))
      base.append(float(sys.argv[2]))
      base.append(float(sys.argv[3]))
      intersectPt=[]      
      intersectPt.append(float(sys.argv[4]))
      intersectPt.append(float(sys.argv[5]))
      intersectPt.append(float(sys.argv[6]))
      print intersectPt

      print >>fw,"appearance { +transparent material {alpha 0.4}}"
      print >> fw, "QUAD"
      print >> fwList, "LIST \n< test.quad"
      print >> fwList, "LIST \n< cube.vect"
      writeSphere(0.1,base,[1,0,0])
      drawCube(base)
      writeSphere(0.1,intersectPt,[1,1,0])
      main()
    else:
      print 'Not enough arguments, first give BASE then the intersection'
    
