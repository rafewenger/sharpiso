#!/usr/bin/python


import glob
files=glob.glob('./temp/*.nrrd')
files2=[]
for l in files:
  if  not(l.endswith('.grad.nrrd') or l.endswith('.cgrad.nrrd')) :
    print l ,'ends with '
    files2.append(l)
print files2
print len(files)
    
str = "this is string example....wow!!!";

suffix = "wow!!!";
print str.endswith(suffix);
print str.endswith(suffix,20);

suffix = "is";
print str.endswith(suffix, 2, 4);
print str.endswith(suffix, 2, 6);
