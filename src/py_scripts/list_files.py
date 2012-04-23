
import glob
import os
import fnmatch

def use_glob():
  cwd =os.getcwd()
  os.chdir("/home/4/wenger/programs/ijk/data/3D/sharp")
  print "new"
  data = raw_input('enter data name example: cube3D, annulus3D ')
  data = data + ".????.nrrd"
  a=[]
  for files in glob.glob(data):
    print files
    a.append(files)
  os.chdir(cwd)
  return a

def use_fn():
  print "here"
  os.chdir("/home/4/wenger/programs/ijk/data/3D/sharp")
  data = raw_input('enter data name ')
  for file in os.listdir('.'):
    if fnmatch.fnmatch(file, data):
        print file    
  return
  
# hello.py
def hello():
    print "Hello World!"
    return
    
"""
  returns all files of type (typ)
  where typ can be, cube3D, annulus3D
"""
def use_glob_to_get_files_of_type (typ):
  cwd =os.getcwd()
  os.chdir("/home/4/wenger/programs/ijk/data/3D/sharp")
  data=typ+".????.nrrd"
  a=[]
  for files in glob.glob(data):
    print files
    a.append(files)
  os.chdir(cwd)
  return a
