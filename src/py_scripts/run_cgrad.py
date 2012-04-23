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
if n!=2:
  print ' example run : python [typ]\
  {cube3D,annulus3D}' 
#  ADD TO THIS LIST FOR THE VARIOUS POSITIONS

position_list = ['gradNS', 'gradEC'] 
#print position_list

for single_file in files:
  #calculate cgrad file name
  cgradfname_list=single_file.split('.')
  cgradfname_list.insert(2,'cgrad')
  cgradfname=cgradfname_list[0]+\
  '.'+cgradfname_list[1]+'.'+cgradfname_list[2]+\
  '.'+cgradfname_list[3]
  cgrad_command=[]
  cgrad_command.append('./cgradient')
  cgrad_command.append('data/'+single_file)
  cgrad_command.append(cgradfname)
  p2=subprocess.Popen(cgrad_command, stdout=subprocess.PIPE) # Success!
  
    

   



