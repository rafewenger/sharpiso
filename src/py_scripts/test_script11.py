 #!/usr/bin/python 
# test isodual 3D 
# for each gradient try all options 
import subprocess as sp 
import sys as s
import random
global run_tests  
global set_isov
global print_res
Stest=True # default is to run the short test
OPTS = []
# set up the names 
opts_name=[]
types=[]
positions=[]
'''
OPTIONS
'''
def set_tests():
  #set up the types  and num 
  #types = ['annulus', 'two_cubes']   
  global positions
  global types 
  types = ['annulus']
  global num
  if (Stest==False) :
    positions = positions+['gradCD', 'gradEC','gradNS', 'gradCDdup','gradIEDir']
    num=50
  else:
    positions = ['gradCD', 'gradEC','gradNS', 'gradCDdup','gradIEDir']
    num = 2
  return 0    
iso_cmd = "isodual3D"
def_parms = ['-trimesh', '-multi_isov','-sep_pos', '-s', '-o', 'positions.off']
'''
set up the isodual3D commend
'''
def setup_isocmd():
  OPTS.append(['-clamp_conflict'])
  opts_name.append('clmConf')
  OPTS.append(['-reposition'])
  opts_name.append('repo')
  OPTS.append(['-no_reposition'])
  opts_name.append('norepo')
  OPTS.append(['-lindstrom'])
  opts_name.append('lind')
  if (True):
    OPTS.append(['-allow_conflict'])
    opts_name.append('alwCnf')
    OPTS.append(['-centroid_conflict'])
    opts_name.append('cntrdCnf')
    OPTS.append(['-clamp_far'])
    opts_name.append('clpFr')
    OPTS.append(['-clamp_far','-max_dist','1.5'])
    opts_name.append('clpFr1.5')
    
    OPTS.append(['-allow_conflict','-Linf','-max_dist','1.5'])
    opts_name.append('01')

    OPTS.append(['-allow_conflict','-Linf','-max_dist','2.0'])
    opts_name.append('02')

    OPTS.append(['-allow_conflict','-Linf','-max_dist','0.5'])
    opts_name.append('03')
  
    #
    OPTS.append(['-clamp_conflict','-Linf','-max_dist','1.5'])
    opts_name.append('04')

    OPTS.append(['-clamp_conflict','-Linf','-max_dist','2.0'])
    opts_name.append('05')

    OPTS.append(['-clamp_conflict','-Linf','-max_dist','0.5'])
    opts_name.append('06')
  
    #
    # #
    OPTS.append(['-allow_conflict','-no_Linf','-max_dist','1.5'])
    opts_name.append('07')

    OPTS.append(['-allow_conflict','-no_Linf','-max_dist','2.0'])
    opts_name.append('08')

    OPTS.append(['-allow_conflict','-no_Linf','-max_dist','0.5'])
    opts_name.append('09')
  
    #
    OPTS.append(['-clamp_conflict','-no_Linf','-max_dist','1.5'])
    opts_name.append('10')

    OPTS.append(['-clamp_conflict','-no_Linf','-max_dist','2.0'])
    opts_name.append('11')

    OPTS.append(['-clamp_conflict','-no_Linf','-max_dist','0.5'])
    opts_name.append('12')
  
    #
    # #
    OPTS.append(['-reselectg'])
    opts_name.append('rslctG')
    OPTS.append(['-allow_conflict','-reposition'])
    opts_name.append('alwCf&repo')
  
  OPTS.append(['-centroid_conflict','-reposition']) 
  opts_name.append('cntrdCnf&repo')
  
  return OPTS


'''
function to set up the isovalues 
'''
def set_isov(typ):
  if typ == "annulus":
    return ['10.1', '10.2','10.4','10.5','10.7']
    #return ['10.1', '10.7']
  if typ == "two_cubes":
    return ['15.1', '15.2','15.5','15.8']  
  
'''
function to set up the calls to isodual3D
'''    
def run_tests():
  row_lists = []
  row = []
  ex = []
  pos_op = []
  num=50
  OPTS = setup_isocmd()
  print 'in run_tests'  
  print 'positions',positions
  print 'opts',OPTS
  print 'num',num
  for typ in types:
    for n in range(1, num):
      filename = typ + str(n) + '.nrrd'
      print 'Running test on ',filename
      isovals = set_isov(typ)
      for iso in isovals:
        for pos in positions:
            row=[]
            row.append(filename)
            row.append(iso)
            row.append(pos)
            opts_list=[]
            for opts in OPTS:
                full_name='testData2/' + filename
                ex=[]
                ex=[iso_cmd]+opts[:]+['-position', pos] + def_parms[:] + [iso, full_name]
                sp.call(ex)
                sp.call(['findedge', '140', 'positions.off'])
                ot = sp.check_output(['findEdgeCount', '-fp', 'positions.line'])
                opts_list.append(ot.split()[1])                
            row.append(opts_list)
            row_lists.append(row)   
  return row_lists 
  
'''
function to set up the calls to isodual3D
'''    
def run_tests_random():
  row_lists = []
  row = []
  ex = []
  pos_op = []

  OPTS = setup_isocmd()
  print 'in run_tests'  
  print 'positions',positions
  print 'opts',OPTS
  print 'num',num
  for typ in types:
    count=0
    while (count <30):
      count=count+1
      n=random.randint(1,50)
      filename = typ + str(n) + '.nrrd'
      print 'Running test on ',filename
      isovals = set_isov(typ)
      for iso in isovals:
        for pos in positions:
            row=[]
            row.append(filename)
            row.append(iso)
            row.append(pos)
            opts_list=[]
            for opts in OPTS:
                full_name='testData2/' + filename
                ex=[]
                ex=[iso_cmd]+opts[:]+['-position', pos] + def_parms[:] + [iso, full_name]
                sp.call(ex)
                sp.call(['findedge', '140', 'positions.off'])
                ot = sp.check_output(['findEdgeCount', '-fp', 'positions.line'])
                opts_list.append(ot.split()[1])                
            row.append(opts_list)
            row_lists.append(row)   
  return row_lists 
def print_res2(res):
    #fi=open ('results_large.txt','w')
    for i in range (len(positions)):
        outname='isotest_'+positions[i]+'.txt'
        fi=open (outname,'w')
        print >>fi,'filename iso',
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
3
'''
def print_res3(res):
    #fi=open ('results_large.txt','w')
    for i in range (len(positions)):
        outname='isotest_csv'+positions[i]+'.txt'
        fi=open (outname,'w')
        cnt_list=[0]*len(OPTS)
        num_test=0
        min0=0
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
print the results 
'''
def print_res(res):
  fi=open ('positions.txt','w')
  print 'printing the results', len(res)

  for row in res:
    print row 
    '''
    vals=row[3]
    vals = [map(int, x) for x in vals]
    names_of_vals=opt_name[:]
    x=zip(vals,names_of_vals)
    x.sort()
    sortedvals,sortednames=zip(*x)  
    print >>fi,row[0], row[1], row[2], sortednames
    '''

'''
main function 
'''
def main ():
  print 'welcome to isodual tests'
  for arg in range(len(s.argv)):
    print 'input arguments ',s.argv[arg]
    if (s.argv[arg]=='-short_test'):
        Stest=True
    else:
        Stest=False
    if (s.argv[arg]=='-help'):
        print 'options : -short_test {to run a short test}'
        
  set_tests()  
  res = run_tests()
  print_res3(res)
  print_res2(res)
  

if __name__ == "__main__":
    main()
    

      
  


    
