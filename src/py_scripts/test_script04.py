#test script4 this is the latest version 
 #!/home/3/bhattaca/bin python

import cmd_line
import create_command
import subprocess
import chart_data

reload (cmd_line)
reload (chart_data)
reload (create_command)

cl=cmd_line
cc=create_command
#parse the commmandline arguments
args=cl.parse_command_line()

#if cgradient then compute the central gradient 
if args.cgradient == True:
  #cgrad_cmd=cc.create_cgrad_command(args)
  #run the cgrad
  print 'create the grad file'

'''
RUN ISODUAL3D ON THE FILES
'''

cmd_list=cc.create_iso_command(args)


if args.key!=None:
  '''
  RUN ISODUAL3D ON THE FILES WITH THE KEY
  '''  
  cmd_list_with_key=cc.create_iso_command_with_key(args,args.key)
  cmd_list=cmd_list + cmd_list_with_key
  
#print 'files',args.input_files
cc.run_commands(cmd_list)

#RUN FIND EDGE ON THE FILES
fe_cmd_list=cc.create_find_edge_command_from_cmd_list(cmd_list) 
cc.run_commands(fe_cmd_list)


#RUN FIND EDGE COUNT ON THE FILES

#print to a file
f=open ('output.txt','w')
fec_cmd_list=cc.create_find_edge_COUNT_command_from_cmd_list(cmd_list)
cc.run_commands_out_to_file(fec_cmd_list,f)



#CHART DATA

f=open('output.txt','r')
print 'reading from file ...'
data=f.readlines()
chart_data.plot_data(data)

