#!/usr/bin/python
#read the vert file 
import sys
def main():
  points=[]
  name=sys.argv[1]
  oname=sys.argv[1].split('.')[0]+'.vect'
  f= open (name,'r')
  fw=open (oname,'w')
  print>>fw,'VECT'
  for index,line in enumerate(f):
    #print line.split()
    pts=' '.join(line.split())  
    points.append(pts)
  print >>fw,index+1,index+1,'1'
  for i in range(index+1):
    print >>fw, '1',
  print >>fw,''
  print >>fw,'1',
  for i in range(index):
    print >>fw, '0',
  print >>fw,''
  for l in points:
    print >> fw , l
  print >>fw,'1 0 0 0.8' 
    

main()
