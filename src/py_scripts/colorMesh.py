#!/usr/bin/python
#read in the off file 
import sys
def main():
  outname=sys.argv[1]
  open_file=open(outname,'r')
  fi=open ('c'+outname,'w')
  file_lines=open_file.readlines()
  line1= file_lines[0].strip()  # First Line
  line2 = file_lines[1].strip('')  # Second Line
  #print line1
  info=line2.split()
  n=int(info[0])+3
  count=0
  for f in file_lines: 
    count=count+1
    #print >> fi, f.strip(),
    fi.write(f.strip())
    if count > n:
      #print >> fi,' 0  255 0 0.1'
      fi.write(' 0  255 0 0.1\n')
    else:
      #print >> fi,''
      fi.write('\n')





main()  
  
