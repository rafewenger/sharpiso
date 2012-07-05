#!/usr/bin/python 
import sys
import struct

#http://docs.python.org/library/struct.html#struct-format-strings



counter=10
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
  readNode(l) 
  
  
   
def readNode(l):
  global counter
  #node 
  for n in range (8):
    print '----- N ',n,
    node_type=convert("i", l,  4)
    print 'node_type',node_type
    
    if node_type==1:
      isempty=convert("h", l,  2)
      print isempty
    
    if node_type==2:
      sign_corner=[]
      for i in range(8):
        sign_corner.append(convert("h", l,  2))
      #print 'sign_cor',sign_corner  
      
      for i in range(12):  
        n=convert("i",l,4)
        intersection_data=[]      
        for j in range(n):
          intersection_data.append(convert("f", l, 4))
          intersection_data.append(convert("f", l, 4))
          intersection_data.append(convert("f", l, 4))
          intersection_data.append(convert("f", l, 4))  
        if len(intersection_data)>0:
          print intersection_data
        
    if node_type==0:
      readNode(l) 
          
      
      
def convert(typ,l, num):  
  global counter
  s=counter
  counter = counter + num
  return   struct.unpack("<"+typ, ''.join(l[s:counter]))[0]
  



# call the main function
if __name__ == "__main__":
    main()
    
