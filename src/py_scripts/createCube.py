#!/usr/bin/python
import sys 

def main ():
  '''
  input 1 2 3 origin location 
  input 4 unit length 
  input 5 cube name
  '''
  origin =[0,0,0]
  origin = [float(sys.argv[1+x]) for x in range (3)]
  ul = float(sys.argv[4])
  fname = sys.argv[5]
  f = open (fname, 'w')
  print >>f, 'LINEC \n1 0 0 1'
  print >>f, '8 12'
  for i in range (2):
    for j in range (2):
      for k in range (2):
        print >>f, origin[0]+i*ul, origin[1]+j*ul, origin[2]+k*ul 
  
  print >>f, '0 1\n1 5\n5 4\n4 6\n6 7\n7 5\n7 3\n3 2\n2 0\n6 2\n1 3\n0 4'


if __name__ == "__main__":
    main()
