import subprocess
import def_param as df
#create the commands from args

def create_iso_command(args):
  cmd_list=[]
  #for each file
  for files in args.input_files:
    #for each isovalue
    for iso in args.isovalue:
      #for each positon type
      for pos in args.position:
         cmd=['./isodual3D']
         if args.max_eigen!=None:
          cmd.append('-max_eigen')
          cmd.append(args.max_eigen)
         if args.max_dist!=None:
          cmd.append('-max_dist')
          cmd.append(args.max_dist)
         if args.reposition== True:
          cmd.append('-reposition')
         if args.lindstrom == True:
          cmd.append('-lindstrom')
         if args.allow_conflict ==True:
          cmd.append('-allow_conflict')
         if args.clamp_conflict == True:
          cmd.append('-clamp_conflict')
         if args.clamp_far ==True:
          cmd.append('-clamp_far')
         if args.centroid_far ==True:
          cmd.append('-centroid_far') 
         if args.recompute_eigen2==True:
          cmd.append('-recompute_eigen2')
         if args.cgradient==True:
          cmd.append('-gradient')
          cgrad_file=files.split('.')
          cgrad_file.insert(len(cgrad_file)-1,'cgrad')
          cmd.append('.'.join(cgrad_file))
         cmd.append('-trimesh')
         cmd.append('-s')
         cmd.append('-position')
         cmd.append(pos)
         cmd.append('-o')
         fl=files.split('.')    
         fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.off'
         fl_name=df.temp_loc+fl_name
         cmd.append(fl_name)
         cmd.append(iso)
         print 'data ', df.data_loc+files
         cmd.append(df.data_loc+files)
         print 'cmd',cmd        
         cmd_list.append(cmd)
  return cmd_list

'''
create commands with key and without key 
'''
def create_iso_command_with_key(args,key):
  cmd_list=[]
  #for each file
  for files in args.input_files:
    #for each isovalue
    for iso in args.isovalue:
      #for each positon type
      for pos in args.position:
         cmd=['./isodual3D']
         if args.max_eigen!=None:
          cmd.append('-max_eigen')
          cmd.append(args.max_eigen)
         if args.max_dist!=None:
          cmd.append('-max_dist')
          cmd.append(args.max_dist)
         if args.reposition== True:
          cmd.append('-reposition')
         if args.lindstrom == True:
          cmd.append('-lindstrom')
         if args.allow_conflict ==True:
          cmd.append('-allow_conflict')
         if args.clamp_conflict == True:
          cmd.append('-clamp_conflict')
         if args.clamp_far ==True:
          cmd.append('-clamp_far')
         if args.centroid_far ==True:
          cmd.append('-centroid_far') 
         if args.recompute_eigen2==True:
          cmd.append('-recompute_eigen2')
         if args.cgradient==True:
          cmd.append('-gradient')
          cgrad_file=files.split('.')
          cgrad_file.insert(len(cgrad_file)-1,'cgrad')
          cmd.append('.'.join(cgrad_file))
         cmd.append('-trimesh')
         #add the key 
         cmd.append('-'+key)
         cmd.append('-position')
         cmd.append(pos)
         cmd.append('-o')
         fl=files.split('.')
         fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+key+'.'+iso+'.off'
         fl_name=df.temp_loc+fl_name
         cmd.append(fl_name)
         cmd.append(iso)
         #DEBUGcmd.append(files)
         cmd.append(df.data_loc+files)
         print cmd        
         cmd_list.append(cmd)
  return cmd_list
'''
create the find Edge command [this is OLD code now i generate it based on the cmd_list]
'''
def create_find_edge_command(args):
  cmd_list=[]
  #for each file
  for files in args.input_files:
    #for each isovalue
    for iso in args.isovalue:
      #for each positon type
      for pos in args.position:
        cmd=['./findedge']
        cmd.append('140')
        fl=files.split('.')
        fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.off'
        fl_name=df.temp_loc.fl_name
        cmd.append(fl_name)
        #cmd.append('.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.off')
        print cmd
        cmd_list.append(cmd)
  return cmd_list
'''
create the find Edge command based on command list 
'''
def create_find_edge_command_from_cmd_list(c_list):
  cmd_list=[]
  for f in c_list:
    cmd=['./findedge']
    cmd.append('140')
    cmd.append(f[len(f)-3])
    print cmd
    cmd_list.append(cmd)
  return cmd_list

'''
CREATE THE FIND EDGE COUNT COMMANDS
'''  
def create_find_edge_count_command(args):
  cmd_list=[]
  #for each file
  for files in args.input_files:
    #for each isovalue
    for iso in args.isovalue:
      #for each positon type
      for pos in args.position:
        cmd=['./findEdgeCount']
        cmd.append('-fp')
        fl=files.split('.')
        fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.line'
        fl_name=df.temp_loc.fl_name
        cmd.append(fl_name)
        #cmd.append('.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.line')
        print cmd
        cmd_list.append(cmd)
  return cmd_list
  
'''
FIND EDGE COUNT BASED ON CMD_LIST
'''
def create_find_edge_COUNT_command_from_cmd_list(c_list):
  cmd_list=[]
  for f in c_list:
    cmd=['./findEdgeCount']
    cmd.append('-fp')
    oname=f[len(f)-3].split('.')
    oname=oname[:len(oname)-1]
    oname='.'.join(oname)+'.line'
    cmd.append(oname)
    cmd_list.append(cmd)
  return cmd_list
  
  
'''
Run a list of commands
'''
def run_commands(cmd_list):
  for cmd in cmd_list:
    p2=subprocess.call(cmd)


  
'''
Run a list of commands using the check output function 
'''
def run_commands_out_to_file(cmd_list,fi):
  for cmd in cmd_list:
    p2=subprocess.check_output(cmd)
    print >>fi,p2

