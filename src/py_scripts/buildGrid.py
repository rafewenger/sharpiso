#!/usr/bin/python 
import sys
import struct

#http://docs.python.org/library/struct.html#struct-format-strings



positions=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
edges= [[0,4],[1,5],[2,6],[3,7],[0,2],[1,3],[4,6],[5,7],[0,1],[2,3],[4,5],[6,7]]
counter=10
level = 0
width = 0
loc = [0,0,0]
fi=open ('vert.off','w')
def main():
  l=[]
  with open(sys.argv[1], "rb") as f:
      byte = f.read(1)
      while byte != "":
          l.append(byte)
          byte = f.read(1) 
  global counter  
  print counter                   
  X = convert("i", l,   4)
  Y = convert("i", l,   4)
  Z = convert("i", l,   4)
  print X, Y,Z, counter
  #node
  node_type = convert ("i", l, 4)
  print 'node_type',node_type
  print 'start'
  readNode(l) 

  
  
   
def readNode(l):
  global counter
  global level
  global loc 
  global width
  origloc = loc
  level = level + 1
  w = 4
  width = w*1.0/2**(level)
  print 'Width',width, 'Level[',level,']'
  
  #node 
  for n in range (8):
    print '----- N ',n,'key positions ',positions[n],
    temploc=[0,0,0]
    width =  w*1.0/2**(level)
    compute_origin(width , n, loc , temploc)
    print 'corner for this cube',temploc,
    node_type=convert("i", l,  4)
    print 'node_type',node_type
 
    if node_type==1:
      isempty=convert("h", l,  2)
      print 'isempty',isempty,
    
    if node_type==2:
      sign_corner=[]
      for i in range(8):
        sign_corner.append(convert("h", l,  2))
      print 'sign_cor',sign_corner  
      
      for i in range(12):  
        ne=convert("i",l,4)
        #print 'ne[',ne,']'
        intersection_data=[] 
        offset=0
        if ne>0:
          print 'Intersection on edge[',i,']',edges[i],
          endPts = edgePoint(i,temploc,width)
          for i in range(2):
            for j in range(3):
              print >>fi,endPts[i][j],
            print >>fi, ''
            
        if ne >=2:
          print 'More than 2 intersections'  
          sys.exit()
        for j in range(ne):
        
          print 'j',j
          a=convert("f", l, 4)
          b=convert("f", l, 4)
          c=convert("f", l, 4)
          d=convert("f", l, 4)
          print '[offset ',a,',',b,',',c,',',d,']'
          intersection_data.append(a)
          intersection_data.append(b)
          intersection_data.append(c)
          intersection_data.append(d)  
 
        #if len(intersection_data)>0:
        #  print intersection_data
        intPoint=[0,0,0]
        if ne==1:
          for d in range(3):
            if endPts[0][d]- endPts[1][d] != 0:
              intPoint[d]= endPts[0][d] + width*intersection_data[0]
            else:
              intPoint[d]= endPts[0][d]
            #intPoint[d]= endPts[0][d] + intersection_data[0]
            #print >>fi,intPoint[d],
          #print >>fi ,' '
        
    if node_type==0:
      print 'calling readNode again'
      loc = temploc
      readNode(l) 
      level = level-1
      loc = origloc

def edgePoint(i, base, w):
  print '{',positions[edges[i][0]],
  print positions[edges[i][1]],'}',
  endPts=[]  
  for j in range(2):
    c = [ w*x for x in positions[edges[i][j]] ]
    #print 'c',c,
    TempendPt=[0,0,0]
    for m in range (3):
      TempendPt[m]=base[m]+c[m]
    #print 'endPt ',TempendPt,
    endPts.append(TempendPt)
  return endPts
  #print base + w*positions[edges[i][0]
         
def convert(typ,l, num):  
  global counter
  s=counter
  counter = counter + num
  return   struct.unpack("<"+typ, ''.join(l[s:counter]))[0]
  

def compute_origin(width , n, origloc, temploc):
  pos = positions[n] # is a 3 tuple
  temploc[0]=origloc[0]+pos[0]*width
  temploc[1]=origloc[1]+pos[1]*width
  temploc[2]=origloc[2]+pos[2]*width

# call the main function
if __name__ == "__main__":
    main()
    
