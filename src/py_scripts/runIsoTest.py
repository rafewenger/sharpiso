#!/usr/bin/python 
# test isodual 3D 
# for each gradient try all options 
import subprocess as sp 
import sys as s
import random
import glob
global run_tests  
global set_isov
global print_res
Stest=False # default is to run the short test
OPTS = []
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
  if (Stest==False) :
    positions = positions+['gradCD', 'gradEC','gradNS', 'gradCDdup','gradIEDir']
    num=50
    isovals = ['5.1', '5.2','5.4','5.5','5.7']
  else:
    positions = ['gradCD', 'gradEC']
    num = 3
    isovals =  ['5.1', '5.2']
  return 0    
iso_cmd = "isodual3D"
def_parms = ['-trimesh', '-resolve_ambig', '-s', '-o', 'test.off']


'''
set up the isodual3D commend
'''
def setup_isocmd():
  '''
  example 
  
  OPTS.append(['-clamp_conflict'])
  opts_name.append('clampC')
  '''
  OPTS.append(['-clamp_conflict'])
  opts_name.append('clampC')
  OPTS.append(['-centroid_conflict'])
  opts_name.append('centC')
  OPTS.append(['-lindstrom'])
  opts_name.append('lnd')
  return OPTS



'''
function to set up the calls to isodual3D
'''    
def run_tests():
  row_lists = []
  row = []
  ex = []
  pos_op = []
  OPTS = setup_isocmd() 
  print 'Positions: ',positions
  print 'Options: ',OPTS
  print 'Isovals', isovals
  files=glob.glob('*.nrrd')
  list_files=[]
  for l in files:
    if  not(l.endswith('.grad.nrrd') or l.endswith('.cgrad.nrrd')) :
      list_files.append(l)
  for filename in list_files:
    print 'Running test on ',filename
    for iso in isovals:
      for pos in positions:
          row=[]
          row.append(filename)
          row.append(iso)
          row.append(pos)
          opts_list=[]
          for opts in OPTS:
              full_name=loc + filename
              ex=[]
              ex=[iso_cmd]+opts[:]+['-position', pos] + def_parms[:] + [iso, full_name]
              sp.call(ex)
              sp.call(['findedge', '140', 'test.off'])
              ot = sp.check_output(['findEdgeCount', '-fp', 'test.line'])
              opts_list.append(ot.split()[1])                
          row.append(opts_list)
          row_lists.append(row)   
  return row_lists 
  
 
def print_res2(res):
    #fi=open ('results_large.txt','w')
    for i in range (len(positions)):
        outname='isotest_'+positions[i]+'.txt'
        fi=open (outname,'w')
        print >>fi,'filename   iso',
        for ele in range(len(OPTS)):
            print >>fi,opts_name[ele].ljust(8),
        print >>fi,''
        cnt_list=[0]*len(OPTS)
        num_test=0
        min0=0
        for row in res:
            if str(row[2])==positions[i]:
                num_test=num_test+1 
                print >>fi,row[0].split('.')[0].ljust(7),
                print >>fi,row[1].ljust(3),' | ',
                for ele in row[3]:
                    print >>fi,ele.ljust(8),          
                vals=row[3]
                vals=map(int,vals)
                minlist = min(vals)
                maxlist = max(vals)
                if minlist == 0:
                    min0=min0+1
                print >>fi,'|',minlist,maxlist
                for ind,n in enumerate(vals):
                    if n==minlist:
                       cnt_list[ind]=cnt_list[ind]+1
        #print 'cnt_list ',cnt_list
        print >>fi,'\n \n SUMMARY'
        print >>fi,'Test on ',positions[i], 'num of tests ',num_test
        for ele in range(len(cnt_list)):
            print >>fi,'%15s was best in %5s cases'% (OPTS[ele], cnt_list[ele])        
        print >>fi,'the min error was 0 in ',min0,'cases'
        print >>fi,'*************************************\n'           

'''
print_res format 3
'''
def print_res3(res):
    #fi=open ('results_large.txt','w')
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
  print_res3(res)
  print_res2(res)
  print ('created files isotest_csv*.txt and isotest_*.txt')
  

if __name__ == "__main__":
    main()
    

      
  


    
