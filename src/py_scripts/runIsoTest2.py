#!/usr/bin/python 
# Run Iso test 2 
# Test isodual 3D 
# for each gradient try all options 
import subprocess as sp 
import sys as s
import random
import glob
global run_tests  
global set_isov
global print_res
global lst 
Stest=False # default is to run the short test
# set up the names 
opts_name=[]
positions=[]
isovals=[]
loc = './'
'''
OPTIONS
'''
def set_tests():
  global positions
  global types 
  global isovals  
  global lst
  if (Stest==False) :
    l1 = ['gradCD',['clmpC','-clamp_conflict'],['allowC','-allow_conflict']]
    l2 = ['gradEC',['lin&dAllowC','-lindstrom','-allow_conflict']]
    lst  = [l1,l2]
    isovals = ['5.1']
  else:
    l1 = ['gradCD',['-clamp_conflict']]
    l2 = ['gradEC',['-lindstrom','-allow_conflict']]
    lst= [l1,l2]
    isovals =  ['5.1', '5.2']
  return 0    
iso_cmd = "isodual3D"
def_parms = ['-trimesh', '-resolve_ambig', '-s', '-o', 'test.off']




'''
function to set up the calls to isodual3D
'''    
def run_tests():
  row_lists = []
  row = []
  ex = []
  pos_op = [] 
  print 'Isovals', isovals
  files=glob.glob('*.nrrd')
  list_files=[]
  for l in files:
    if  not(l.endswith('.grad.nrrd') or l.endswith('.cgrad.nrrd')) :
      list_files.append(l)
  for filename in list_files:
    print 'Running test on ',filename
    for iso in isovals:
      row=[]
      row.append(filename)
      row.append(iso)
      
      for l in lst:
          #print l
          pos = l[0]
          OPTS = l[1:]
          row.append(pos)
          opts_results=[]
          for op in OPTS:         
            full_name= loc + filename
            ex=[]
            ex=[iso_cmd]+op[1:]+['-position', pos] + def_parms[:] + [iso, full_name]
            sp.call(ex)
            sp.call(['findedge', '140', 'test.off'])
            ot = sp.check_output(['findEdgeCount', '-fp', 'test.line'])
            opts_results.append(ot.split()[1])
            row.append(ot.split()[1])
      row_lists.append(row)          
  return row_lists 
  
 
def print_res2(res):
    fi=open ('iso_compile_res.csv','w')
    print >>fi, 'Filename,', 'isoval,',
    for l in lst:
      print >>fi, l[0],',',
      for ops in l[1:]:
        print >>fi, ops[0],',',
    print >>fi,''
    for row in res:
      str_row = ','.join (row)
      print >>fi,str_row
'''
print_res format 3
'''
def print_res3(res):
    fi=open ('results_large.txt','w')
    for i in range (len(positions)):
        outname='isotest_csv'+positions[i]+'.txt'
        fi=open (outname,'w')
        cnt_list=[0]*len(OPTS)
        num_test=0
        min0=0
        print >>fi,'filename, iso,',
        for n in range(len(opts_name)):
          print >>fi,opts_name[n],',',
        #line break after printing the headers
        print >>fi,' '
        for row in res:
            
            if str(row[2])==positions[i]:
                num_test=num_test+1 
                print >>fi,row[0].split('.')[0],',',
                print >>fi,row[1].ljust(3),',',
                for ele in row[3]:
                    print >>fi,ele,',',          
                vals=row[3]
                vals=map(int,vals)
                minlist = min(vals)
                maxlist = max(vals)
                if minlist == 0:
                    min0=min0+1
                print >>fi,minlist,',',maxlist
                for ind,n in enumerate(vals):
                    if n==minlist:
                       cnt_list[ind]=cnt_list[ind]+1
 

'''
main function 
'''
def main ():
  print 'welcome to isodual tests'
  global Stest
  for arg in range(len(s.argv)):
    if (s.argv[arg]=='-short_test'):
        Stest=True
    if (s.argv[arg]=='-help'):
        print 'options : -short_test {to run a short test}'
        print '          -loc { location to look for files}'
        s.exit(0)
    if (s.argv[arg] == '-loc'):
        loc = argv[arg+1]
        print 'looking for files in ',loc
        
  set_tests()  
  res = run_tests()
  print_res2(res)
  print ('')
  

if __name__ == "__main__":
    main()
    

      
  


    
