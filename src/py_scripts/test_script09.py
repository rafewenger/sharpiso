# test script 09 for create_commmand 04 
import cmd_line as cl
import create_command04 as cc
import def_param as df
import time 

fi=open ('testData/data_sheet.txt','r')
flist=[]
# set up the default parameters
for files in fi:
  a=files.split()
  flist.append(a[0])
  
args=cl.parse_command_line()

#args.input_files=flist[1:40]  # for annulus 
#args.input_files=flist[105:130]  # for twocubes 
args.input_files=flist[105:106]  # for twocubes

#dt=str(dt.datetime.now())
localtime=time.localtime(time.time())
dt=str(localtime[2])+"-"+str(localtime[1])\
+"-"+str(localtime[3])+"-"+str(localtime[4])
tab_width=6


op=cc.create_iso_command(args)
for iv,out in enumerate(op):
  fname='out'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv],
  print >>fi,'**************ISODUAL3D with perfect gradients **************'
  print >>fi,df.TEST_HEADER_STRING1
  print >>fi,df.TEST_HEADER_STRING2
  for out_sub in out:
    if out_sub != 'break':
      print >>fi,out_sub.ljust(tab_width,' '),
    else:
      print >>fi,'\n',
  print 'output written to',fname,']' 
  


   
################################## cgradient
op=cc.create_iso_command_cgrad(args)
for iv,out in enumerate(op):
  fname='outc'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv],
  print >>fi,'**************ISODUAL3D with cgrad **************'
  print >>fi,df.TEST_HEADER_STRING1
  print >>fi,df.TEST_HEADER_STRING2
  for out_sub in out:
    if out_sub != 'break':
      print >>fi,out_sub.ljust(tab_width,' '),
    else:
      print >>fi,'\n',
  print 'output written to',fname,']' 
  

################################## negative gradient
import create_command05 as cc5
args.position =  ['gradCD','gradIES','gradIEDir']
args.isovalue =  ['-10.1','-10.3','-10.6']
op=cc5.create_iso_command_negate_sep_neg(args) 
for iv,out in enumerate(op):
  fname='outsetneg'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv],
  print >>fi,'**************ISODUAL3D with negative gradients separate negatives **************'
  print >>fi,df.TEST_HEADER_STRING1
  print >>fi,df.TEST_HEADER_STRING2
  for out_sub in out:
    if out_sub != 'break':
      print >>fi,out_sub.ljust(tab_width,' '),
    else:
      print >>fi,'\n',
  print 'output written to',fname,']' 
  
###############
op=cc5.create_iso_command_negate_sep_pos(args) 
for iv,out in enumerate(op):
  fname='outsetpos'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv],
  print >>fi,'**************ISODUAL3D with negative gradients separate positives **************'
  print >>fi,df.TEST_HEADER_STRING1
  print >>fi,df.TEST_HEADER_STRING2
  for out_sub in out:
    if out_sub != 'break':
      print >>fi,out_sub.ljust(tab_width,' '),
    else:
      print >>fi,'\n',
  print 'output written to',fname,']' 
  

  