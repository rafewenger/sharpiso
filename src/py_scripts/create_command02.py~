import subprocess
import def_param as df
#create the commands from args
fi=open ('output.txt','w')
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
         if args.create_files==True:
          fl=files.split('.')    
          if args.key!=None:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+args.key+'.'+iso+'.off'
          else:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.off'
          fl_name=df.temp_loc+fl_name
          cmd.append(fl_name)
         else:
          cmd.append('out.off')
         cmd.append(iso)
         print 'data ', df.data_loc+files
         cmd.append(df.data_loc+files)
         print 'cmd',cmd       
         cmd_list.append(cmd)
         p2=subprocess.call(cmd)
         
         # now for findEdge
         cmd=['./findedge']
         cmd.append('140')
         if args.create_files==True:
          fl=files.split('.')
          if args.key!=None:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+args.key+'.'+iso+'.off'
          else:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.off'
          fl_name=df.temp_loc+fl_name
          cmd.append(fl_name)
         else:
          cmd.append('out.off')
         print 'findedge',cmd
         p2=subprocess.call(cmd)
         # findedge end
         
         # now for findEdgeCOunt
         cmd=['./findEdgeCount']
         cmd.append('-fp')
         if args.create_files==True:
          fl=files.split('.')
          if args.key!=None:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+args.key+'.'+iso+'.line'
          else:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso+'.line'
          fl_name=df.temp_loc+fl_name
          cmd.append(fl_name)
         else:
          cmd.append('out.line')
         print 'findedgecount',cmd
         p2=subprocess.check_output(cmd)
         
         fl=files.split('.')
         if args.key!=None:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+args.key+'.'+iso
         else:
            fl_name ='.'.join(fl[:(len(fl)-1)])+'.'+pos+'.'+iso
         fl_name=df.temp_loc+fl_name
         ot=p2.split()
         print >>fi, fl_name,ot[1]
         # findedge count end
         
         
  return cmd_list      
         
