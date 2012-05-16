
import subprocess as sp

args =['../ijkgenscalar2', '-grad', '-dim', '3', '-asize', '100']
cgrad=['../cgradient']
dirs=['1 0 0 ','1 1 0 ','1 1 1','2 1 0','2 1 1','2 2 1','3 1 0','3 1 1','3 2 0','3 2 1','3 2 2','3 3 1','3 3 2']
cntrs=['50 50 50','50.1 50.2 50.3','50.2 50.3 50.4','50.3 50.4 50.5']
fi=open ('data_sheet.txt','w')

n=1
# annulus data set
for cen in cntrs:
  temp=args[:]
  temp.append('-field')
  temp.append('annulus')
  temp.append('-center')
  temp.append(cen)
  for di in dirs:
    temp2=temp[:]	
    temp2.append('-dir')
    temp2.append(di)
    temp2.append('-ratio')
    temp2.append('1')
    temp2.append('-radius')
    temp2.append('20')
    flname='annulus'+str(n)+'.nrrd'
    temp2.append(flname)
    n=n+1
    print >>fi, flname, 'dir ',str(di),' center ',str(cen)
    p1=sp.call(temp2)
    ctemp=cgrad[:]
    ctemp.append(flname)
    cgrad_fname='annulus'+str(n-1)+'.cgrad'+'.nrrd'
    ctemp.append(cgrad_fname)
    p2=sp.call(ctemp)
    
n=1    
# flange data set    
for cen in cntrs:
  temp=args[:]
  temp.append('-field')
  temp.append('flange')
  temp.append('-center')
  temp.append(cen)
  for di in dirs:
    temp2=temp[:]	
    temp2.append('-dir')
    temp2.append(di)
    temp2.append('-ratio')
    temp2.append('1')
    temp2.append('-radius')
    temp2.append('20')
    flname='flange'+str(n)+'.nrrd'
    temp2.append(flname)    
    n=n+1
    print >>fi, flname, 'dir ',str(di),' center ',str(cen)
    p1=sp.call(temp2)
    ctemp=cgrad[:]
    ctemp.append(flname)
    cgrad_fname='flange'+str(n-1)+'.cgrad'+'.nrrd'
    ctemp.append(cgrad_fname)
    p2=sp.call(ctemp)
    
    
n=1    
cntrs1=['30 35 30 60 55 65', '35.1 40.2 38.3 61.2 60.3 62.4']

# twocubes data set    
for cen in cntrs1:
  temp=args[:]
  temp.append('-field')
  temp.append('cube')
  temp.append('-two')
  temp.append('-center')
  temp.append(cen)
  for di in dirs:
    temp2=temp[:]	
    temp2.append('-dir')
    temp2.append(di)
    temp2.append('-side_dir')
    temp2.append('0 1 1')
    flname='two_cubes'+str(n)+'.nrrd'
    temp2.append(flname)    
    n=n+1
    print >>fi, flname, 'dir ',str(di),' center ',str(cen)
    p1=sp.call(temp2)
    ctemp=cgrad[:]
    ctemp.append(flname)
    cgrad_fname='two_cubes'+str(n-1)+'.cgrad'+'.nrrd'
    ctemp.append(cgrad_fname)
    p2=sp.call(ctemp)
