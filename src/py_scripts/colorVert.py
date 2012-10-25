#!/usr/bin/python
#read in the off file 
import sys
import random
def main():
  outname=sys.argv[1]
  open_file=open(outname,'r')
  fi=open ('c'+outname,'w')
  file_lines=open_file.readlines()
  line1= file_lines[0].strip()  # First Line
  line2 = file_lines[1].strip('')  # Second Line
  print line1
  info=line2.split()
  print info
  n=int(info[0])+2 
  count=0
  for f in file_lines: 
    count=count+1
    if count <= n and count > 2:
      #r=random.randint(0, 255)
      #g=random.randint(0, 255)
      #b=random.randint(0, 255)
      r = 1.0
      g = 0.0
      b = 0.0
      print r,g,b
      print >> fi,f.strip(),r,g,b,'1.0'
    elif count == 1:
      print >>fi, 'COFF'
    else: 
      print >> fi, f.strip()



main()  
  
