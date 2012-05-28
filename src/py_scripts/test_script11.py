# test isodual 3D 
# for each gradient try all options 
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
#positions = ['gradCD', 'gradEC','gradNS','gradIES','gradIEDir','gradCS','centroid']
positions = ['gradCD', 'gradEC']
iso_cmd = "isodual3D"
def_parms = ['-trimesh', '-s', '-o', 'out.off']
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
    #return ['10.1', '10.2','10.4','10.7']
    return ['10.1']
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
  num = 3
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
                sp.call(['findedge', '140', 'out.off'])
                ot = sp.check_output(['findEdgeCount', '-fp', 'out.line'])
                opts_list.append(ot.split()[1])                
            row.append(opts_list)
            row_lists.append(row)
        '''
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
         '''   
  return row_lists 
  



'''
print the results 
'''
def print_res(res):
  fi=open ('result1.txt','w')
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
    
    
  

'''
main function 
'''
def main ():
  print 'welcome to isodual tests'
  res = run_tests()
  print_res(res)
  

if __name__ == "__main__":
    main()
    

      
  


    
