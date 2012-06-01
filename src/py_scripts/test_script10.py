# test isodual 3D
import subprocess as sp 
global run_tests  
global set_isov
global print_res
'''
OPTIONS
'''
#set up the types  and num 
#types = ['annulus', 'two_cubes']
types = ['annulus']
positions = ['gradCD', 'gradEC','gradNS','gradIES','gradIEDir','gradCS','centroid']
#positions = ['gradCD', 'gradEC','gradNS']
iso_cmd = "isodual3D"
def_parms = ['-trimesh', '-multi_isov','-sep_pos', '-s', '-o','out.off']
OPTS = []
'''
set up the isodual3D commend
'''
def setup_isocmd():
  OPTS.append(['-clamp_conflict'])
  OPTS.append(['-reposition'])
  OPTS.append(['-lindstrom'])
  OPTS.append(['-allow_conflict'])
  OPTS.append(['-centroid_conflict'])
  OPTS.append(['-clamp_far'])
  OPTS.append(['-reselectg'])
  OPTS.append(['-allow_conflict','-reposition'])
  return OPTS


'''
function to set up the isovalues 
'''
def set_isov(typ):
  if typ == "annulus":
    return ['10.1', '10.2','10.4','10.7']
    #return ['10.1', '10.8']
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
  num = 30
  OPTS = setup_isocmd()

  for typ in types:
    for n in range(1, num):
      filename = typ + str(n) + '.nrrd'
    
      isovals = set_isov(typ)
      for iso in isovals:
        for opts in OPTS:
        #print 'iso ',iso,' opts ',opts
          #row = [:]
          row=[]
          row.append(filename)
          row.append(iso)
          row.append(opts)
          #print 'row',row 
          #print filename,iso,opts,'\n'
          pos_op=[]
          for pos in positions:
            full_fname = 'testData/' + filename
            ex[:] = []
            ex = [iso_cmd] + opts[:] + ['-position', pos] + def_parms[:] + [iso, full_fname]
            sp.call(ex)  
            sp.call(['findedge', '140', 'out.off'])
            ot = sp.check_output(['findEdgeCount', '-fp', 'out.line'])
            pos_op.append(ot.split()[1])
          row.append(pos_op)
          row_lists.append(row)
          #print 'row_list ',row_lists 
  return row_lists 
  



'''
print the results 
'''
def print_res(res):
  fi=open ('result.txt','w')
  print 'printing the results', len(res)
  print >>fi,'filename'.center(15),'iso','options'.center(30),
  for p in positions:
    print >>fi,p.center(7),
  print >>fi,'\n'
  for row in res:
    #debug 
    print 'row',row
    print >>fi,row[0].ljust(15),row[1].ljust(3),
    temp=''
    for opt in row[2]:
      temp=temp+opt
    print >>fi,temp.ljust(30),
    for val in row[3]:
      print >>fi,val.ljust(7), 
    print >>fi,'\n',
 

    
def print_res2(res):
    ranklist = [0]*len(OPTS)
    fi = open ('results_summary.txt','w')
    
    for i in range(len(positions)):
        ####################
        #setup the counters for the opts 
        clmpconflict_cnt=0
        repo_cnt=0
        allw_conflct_cnt=0
        lind_cnt=0
        cntrd_conflct=0
        clamp_far_cnt=0
        reselectg_cnt=0
        all_conf_repo_cnt=0
        min0=0
        ####################
        for row in res: 
            pos_vals=row[3]
            opts_name=(' ').join(row[2])
            #print 'name',opts_name
            vals = map(int,pos_vals)
            minlist= min(vals)
            maxlist= max(vals)
            if minlist == vals[i] and minlist!=maxlist:
                if (opts_name=='-clamp_conflict'):
                    clmpconflict_cnt=clmpconflict_cnt+1
                elif (opts_name=='-reposition'):
                    repo_cnt=repo_cnt+1
                elif (opts_name=='-allow_conflict'):
                    allw_conflct_cnt = allw_conflct_cnt+1
                elif (opts_name=='-lindstrom'):
                    lind_cnt = lind_cnt+1
                elif (opts_name=='-centroid_conflict'):
                    cntrd_conflct=cntrd_conflct+1
                elif (opts_name=='-clamp_far'):
                    clamp_far_cnt=clamp_far_cnt+1
                elif (opts_name=='-reselectg'):
                    reselectg_cnt=reselectg_cnt+1
                elif (opts_name=='-allow_conflict -reposition'):
                    all_conf_repo_cnt=all_conf_repo_cnt+1
            if minlist == 0:
                min0=min0+1
            
        print >>fi,'for ',positions[i],' out of ',len(res),'tests'
        print >>fi,'-clamp_conflict was at the top for'.ljust(40),clmpconflict_cnt
        print >>fi,'-reposition was best for '.ljust(40),repo_cnt
        print >>fi,'-lindstrom was best for'.ljust(40),lind_cnt
        print >>fi,'-allow_conflict was best for'.ljust(40),allw_conflct_cnt
        print >>fi,'-centroid_conflict was best for'.ljust(40),cntrd_conflct
        print >>fi,'-clamp_far was best for '.ljust(40),clamp_far_cnt
        print >>fi,'-reselectg was best for'.ljust(40),reselectg_cnt
        print >>fi,'-allow_conflict -reposition was best for'.ljust(40),all_conf_repo_cnt
        print >>fi,'the number of errors were 0 in times'.ljust(40),min0
        print >>fi,'************************************************'
            
            
            

'''
main function 
'''
def main ():
  print 'welcome to isodual tests'
  res = run_tests()
  print_res(res)
  print_res2(res)
  

if __name__ == "__main__":
    main()
    


##OBSOLETE     
def print_res3(res):
  fi=open ('result.txt','w')
  print 'printing the results', len(res)
  for row in res:
    print row
    vals=row[3]
    vals = [map(int, x) for x in vals]
    names_of_vals=positions[:]
    x=zip(vals,names_of_vals)
    x.sort()
    sortedvals,sortednames=zip(*x)  
    print >>fi,row[0], row[1].ljust(3), ', '.join(row[2]).ljust(30),
    for y in sortednames:
        print >>fi,"%7s"% y,    
    print >>fi,''    

      
  


    
