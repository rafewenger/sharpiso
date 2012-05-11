#!/home/3/bhattaca/bin/python
# not the latest version 

import cmd_line
import create_command02
import subprocess
import def_param as df

reload (cmd_line)
reload (create_command02)
cl=cmd_line

cc=create_command02
#parse the commmandline arguments
args=cl.parse_command_line()
if args.long==True:
  print'Run more extensive tests ...'
  args.position = df.POS_DEF_LONG
  args.isovalue = df.ISO_DEF_LONG
  args.input_files = df.DATA_DEF_LONG
# run the commands 
# If the key is on then we run it twice
if args.key !=None:
  cmd_list=cc.create_iso_command(args)
  args.key=None
  
cmd_list=cc.create_iso_command(args)
  

