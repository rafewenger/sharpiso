#run the test_script06 over a bunch of test data sets
import subprocess as sp

args=['python','test_script06.py']
pos=['gradEC','gradCD','gradCDdup','gradNS']


for p in pos:
  output=[]
  print 'running with on ' , p
  ofname= "out_annulus_"+p +".txt"
  fi=open (ofname,'w')
  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-reposition')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  fi.write ('-----reposition-----\n')
  for line in f1:
    print >>fi,line,
    output.append(line)
  f1.close()    

  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-allow_conflict')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  
  fi.write ('-----allow_conflict-----\n')
  for line in f1:
    print >>fi,line,
    output.append(line)
  f1.close()
'''
  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-clamp_far')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  
  fi.write ('-----clamp far-----\n')
  for line in f1:
    print >>fi,line,
  f1.close()
  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-reselectg')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  
  fi.write ('-----reselectg-----\n')
  for line in f1:
    print >>fi,line,

# WITH CGRADIENT 
args=['python','test_script06.py','-cgradient']
for p in pos:
  print 'running with cgradient on ' , p
  ofname= "out_annulus_c"+p +".txt"
  fi=open (ofname,'w')
  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-reposition')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  fi.write ('-----with cgradient reposition-----\n')
  for line in f1:
    print >>fi,line,
  f1.close()    
  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-allow_conflict')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  fi.write ('-----with cgradient allow_conflict-----\n')
  for line in f1:
    print >>fi,line,
  f1.close()
  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-clamp_far')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  fi.write ('-----with cgradient reposition-----\n')
  for line in f1:
    print >>fi,line,
  f1.close()
  temp_args=args[:]
  temp_args.append('-position')
  temp_args.append(p)
  temp_args.append('-reselectg')
  p1 = sp.call (temp_args)
  #read the output.txt
  f1 = open ('output.txt','r')
  fi.write ('-----with cgradient reselectg-----\n')
  for line in f1:
    print >>fi,line,
'''     
print 'output', output
print 'DONE' 
