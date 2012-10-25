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
  print >>f, 'OFF'
  print >>f, '8 12 0'
  for i in range (2):
    for j in range (2):
      for k in range (2):
        print >>f, origin[0] + i*ul, origin[1] + j*ul, origin[2] +  k*ul 
  
  print >>f, '\n3 0 1 2\n3 1 2 3\n3 0 1 4\n3 4 1 5\n3 6 7 3\n3 3 6 2\n3 5 1 7\n3 1 7 3\n3 4 0 6\n3 0 6 2\n3 4 5 6\n3 6 5 7'


if __name__ == "__main__":
    main()
