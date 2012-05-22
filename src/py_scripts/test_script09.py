# test script 09 for create_commmand 04 
import cmd_line as cl
import create_command04 as cc
import def_param as df
import datetime as dt
import datetime

args=cl.parse_command_line()

dt=str(dt.datetime.now())

tab_width=6

#debug 
print 'original options', args.position
args.position=['gradEC', 'gradNS']
print '=====output==='
op=cc.create_iso_command(args)
for iv,out in enumerate(op):
  fname='out'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv]
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
print '=====output cgradient==='

for iv,out in enumerate(op):
  fname='outc'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv]
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
op=cc5.create_iso_command_negate(args)
print '=====output negative gradient==='
args.position =  ['gradCD','gradIES','gradIEDir']
args.isovalue =  ['-10.1','-10.6'] 

for iv,out in enumerate(op):
  fname='outneg'+args.position[iv]+dt+'.txt'
  fi=open (fname,'w')
  print >>fi,args.position[iv]
  print >>fi,df.TEST_HEADER_STRING1
  print >>fi,df.TEST_HEADER_STRING2
  for out_sub in out:
    if out_sub != 'break':
      print >>fi,out_sub.ljust(tab_width,' '),
    else:
      print >>fi,'\n',
  print 'output written to',fname,']' 
  


  
'''
################################# neg 
import create_command04 as cc
print '=====output negative==='
fname='outneg'+dt+'.txt'
op=cc.create_iso_command_neg(args) 
fi=open (fname,'w')
for out_sub in op:
  if out_sub !='break':
    print >>fi,out_sub.ljust(tab_width,' '),
  else:
    print >>fi,'\n',
print 'output written to',fname,']'
'''  