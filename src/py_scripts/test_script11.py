# test isodual 3D 
# for each gradient try all options 
import subprocess as sp 
global run_tests  
global set_isov
global print_res
# set up the names 
opt_name=['clmp_cnflct','repo','lind','allow_cnflct','cntrd_cnflct','clmp_far','rslctg'\
  ,'all_cnlf&repo']

'''
OPTIONS
'''
#set up the types  and num 
#types = ['annulus', 'two_cubes']
types = ['annulus']
positions = ['gradCD', 'gradEC','gradNS','gradIES','gradIEDir','gradCS','centroid']

iso_cmd = "isodual3D"
def_parms = ['-trimesh', '-multi_isov','-sep_pos', '-s', '-o', 'positions.off']
'''
set up the isodual3D commend
'''
def setup_isocmd():
  OPTS = []
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
        for pos in positions:
            row=[]
            row.append(filename)
            row.append(iso)
            row.append(pos)
            opts_list=[]
            for opts in OPTS:
                full_name='testData/' + filename
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
print the results 
'''
def print_res(res):
  fi=open ('positions.txt','w')
  print 'printing the results', len(res)

  for row in res:
    vals=row[3]
    vals = [map(int, x) for x in vals]
    names_of_vals=opt_name[:]
    x=zip(vals,names_of_vals)
    x.sort()
    sortedvals,sortednames=zip(*x)  
    print >>fi,row[0], row[1], row[2], sortednames


'''
main function 
'''
def main ():
  print 'welcome to isodual tests'
  res = run_tests()
  print_res(res)
  

if __name__ == "__main__":
    main()
    

      
  


    
