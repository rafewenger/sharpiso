import sys
import os
import list_files
import shlex
import subprocess
import glob
reload (list_files)
#print sys.argv
n=len(sys.argv)
typ = sys.argv[n-1]  # last argument cube3D annulus3D ...
files=list_files.use_glob_to_get_files_of_type (typ)
#isovalue 
isoval=sys.argv[n-2]
#print 'isovalue',isoval
#print sys.argv[1:n-1]
#print os.getcwd()
#  ADD TO THIS LIST FOR THE VARIOUS POSITIONS

position_list = ['gradNS', 'gradEC'] 
#print position_list

for single_file in files:
  cmd = sys.argv[1:n-2]
  # for each position in position list 
  for pos in position_list:
    updated_cmd = cmd[:]
    updated_cmd.append('-position')
    updated_cmd.append(pos)
    updated_cmd.append('-trimesh')
    updated_cmd.append('-o')
    outfile=pos+"."+single_file+".off"
    updated_cmd.append(outfile)
    updated_cmd.append(isoval) # update isovalue
    #add path to data
    file_path = 'data/'+single_file
    updated_cmd.append(file_path)
    updated_cmd.insert(0,'./isodual3D')
    #print 'command ',updated_cmd
    p=subprocess.Popen(updated_cmd) # Success!
    
    

   



