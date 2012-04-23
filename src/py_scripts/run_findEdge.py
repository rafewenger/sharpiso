import sys
import os
import list_files
import shlex
import subprocess
import glob
reload (list_files)
for files in glob.glob('*.off'):
  findedge_cmd = ['./findedge']
  findedge_cmd.append('140')
  findedge_cmd.append(files)
  #print findedge_cmd
  p=subprocess.Popen(findedge_cmd) # Success!

'''
for  files in glob.glob('*.line'):
  findEdgeCount_cmd = ['./findEdgeCount']
  findEdgeCount_cmd.append(files)
  p2=subprocess.Popen(findEdgeCount_cmd,stdout=subprocess.PIPE) # Success!
  out= p2.communicate() 
  print out[0]
'''  
