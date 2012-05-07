# program to parse the commmand line arguments
import argparse
import def_param as df

def parse_command_line():
  parser = argparse.ArgumentParser()
  
  parser.add_argument('-position','-pos',nargs='*', default=df.POS_DEF)
  parser.add_argument('-cgradient',  action='store_true')
  parser.add_argument('-max_eigen')
  parser.add_argument('-max_dist')
  parser.add_argument('-create_files',action='store_true')
  parser.add_argument('-key')
  parser.add_argument('-reposition',action='store_true')
  parser.add_argument('-lindstrom',action='store_true')
  parser.add_argument('-allow_conflict',action='store_true')
  parser.add_argument('-clamp_conflict',action='store_true')
  parser.add_argument('-clamp_far',action='store_true')
  parser.add_argument('-recompute_eigen2',action='store_true')
  parser.add_argument('-isovalue','-isoval',nargs='*', default=df.ISO_DEF)
  parser.add_argument('-input_files', nargs='+', default=df.DATA_DEF)
  
  
  args=parser.parse_args()
  return args
