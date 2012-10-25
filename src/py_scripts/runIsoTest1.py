#!/usr/bin/python 
import sys
import subprocess as sp
pos_set = False 
cmd=['./isodual3D','-trimesh','-o','out.off'] 
res=[]

for n in range (len (sys.argv)):
  if n != 0:
    cmd.append(sys.argv[n])    
#isodual3D
sp.call(cmd)

cmd_findEdge=['findedge','140','out.off']
sp.call (cmd_findEdge)
cmd_findEdgeCount=['findEdgeCount','-fp','out.line']
a=sp.check_output (cmd_findEdgeCount)


res.append(a.split(' ')[1]) 


print res

