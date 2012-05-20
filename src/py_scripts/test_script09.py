# test script 09 for create_commmand 04 
import cmd_line as cl
import create_command04 as cc
import def_param as df
import datetime as dt
import datetime

args=cl.parse_command_line()
args.position=['gradEC','gradNS']
op=cc.create_iso_command(args)
dt=str(dt.datetime.now())


print '=====output==='

for iv,out in enumerate(op):
  fname='out'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv]
  print >>fi,df.tests
  for out_sub in out:
    if out_sub != 'break':
      print >>fi,out_sub.ljust(5,' '),
    else:
      print >>fi,'\n',
  print 'output written to',fname,']' 
  


      
################################## cgradient
op=cc.create_iso_command_cgrad(args)
print '=====output cgradient==='

for iv,out in enumerate(op):
  fname='outc'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv]
  print >>fi,df.tests
  for out_sub in out:
    if out_sub != 'break':
      print >>fi,out_sub.ljust(5,' '),
    else:
      print >>fi,'\n',
  print 'output written to',fname,']' 