import sys
import struct
'''
the input name should not have any '.'s
'''

class Vertex:
  x=0
  y=0 
  z=0
class Face:
  numSides = 0
  indices=[]
    
def main():
  fi = open (sys.argv[1],'r')
  fileinput = sys.argv[1]
  fileoutput = sys.argv[2]
  
  fw = open (fileoutput, 'w')
  
  allLines = fi.read();
  
  f=open(sys.argv[1],'rb')
  lines = f.readlines()
  numv = int(lines[2].split(' ')[2])
  
  numf = int(lines[6].split(' ')[2])
  
  
  print >>fw,'OFF'
  print >>fw,numv, numf , 0
  l=0
  for i in range (9):
    l = l + len (lines[i])
  
  vertexLocation=[]  
  indicesList=[]
  
  for i in range (numv):
    v = Vertex()
    v.x = struct.unpack('>f', allLines[l:l+4])[0]
    l = l+4
    v.y = struct.unpack('>f', allLines[l:l+4])[0]
    l = l+4
    v.z = struct.unpack('>f', allLines[l:l+4])[0]
    l = l+4
    vertexLocation.append(v)

  for i in range(numv):
    v = Vertex()
    v = vertexLocation[i]
    print >>fw,v.x,v.y,v.z


  for i in range (numf):
    f = Face()
    n = int(struct.unpack('>B', allLines[l:l+1])[0])
    f.numSides= n
    l=l+1
    temp =[]
    for j in range (f.numSides):
      temp.append(struct.unpack('>i', allLines[l:l+4])[0])
      f.indices= temp
      l = l+4
    indicesList.append(f)
   
   
  for i in range (numf):
    f = Face()
    f=indicesList[i] 
    print >>fw,f.numSides,
    for j in range (f.numSides):
      print >>fw,f.indices[j],
    print >>fw,''
    
  
if __name__ == "__main__":
    main()
