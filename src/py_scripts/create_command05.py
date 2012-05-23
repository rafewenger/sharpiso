# create command05 this is to run the negate commands 
import subprocess as sp

def create_iso_command_negate_sep_neg(args):
  outlist=[]
  #read the data sheet info 
  fda=open ('testData/data_sheet.txt','r')
  print '**************ISODUAL3D with negative gradients separate negatives **************'
  for pos in args.position:
    output=[]
    print 'position ',pos
    for infile in args.input_files:
      #output.append(infile)
      tmp_fda=fda.readline()
      tmp_fda=tmp_fda.split(" ")
      dir=tmp_fda[3]+tmp_fda[4]+tmp_fda[5]
      for iso in args.isovalue:
        infile_temp=infile.split('.')
        output.append(infile_temp[0])
        output.append(dir)
        output.append(iso)
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #repostion
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-reposition', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #lindstrom
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-lindstrom', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow_conflict
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-allow_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #centroid_conflcit
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-centroid_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp_far
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-clamp_far', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #reselectg
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-reselectg', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow conflict + repos
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-reposition','-allow_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp conflict + repos
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_neg','-s','-position',pos,'-o', 'out.off','-reposition','-clamp_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        output.append('break')
      print '--',infile,'done'
    outlist.append(output)
  return outlist

############### seperate the positives
def create_iso_command_negate_sep_pos(args):
  outlist=[]
  #read the data sheet info 
  fda=open ('testData/data_sheet.txt','r')
  print '**************ISODUAL3D with negative gradients separate positives **************'
  for pos in args.position:
    output=[]
    print 'position ',pos
    for infile in args.input_files:
      #output.append(infile)
      tmp_fda=fda.readline()
      tmp_fda=tmp_fda.split(" ")
      dir=tmp_fda[3]+tmp_fda[4]+tmp_fda[5]
      for iso in args.isovalue:
        infile_temp=infile.split('.')
        output.append(infile_temp[0])
        output.append(dir)
        output.append(iso)
        sp.call(['./isodual3D','-multi_isov','-sep_pos','-trimesh','-s','-position',pos,'-o', 'out.off', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #repostion
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-reposition', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #lindstrom
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-lindstrom', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow_conflict
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-allow_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #centroid_conflcit
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-centroid_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp_far
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-clamp_far', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #reselectg
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-reselectg', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow conflict + repos
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-reposition','-allow_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp conflict + repos
        sp.call(['./isodual3D','-trimesh','-multi_isov','-sep_pos','-s','-position',pos,'-o', 'out.off','-reposition','-clamp_conflict', iso,\
                 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        output.append('break')
      print '--',infile,'done'
    outlist.append(output)
  return outlist
############### run with cgradient
def create_iso_command_negate_cgrad(args):
  outlist=[]
  #read the data sheet info 
  fda=open ('testData/data_sheet.txt','r')
  print 'create_iso_command_cgrad'
  for pos in args.position:
    output=[]
    print 'position cgradient',pos
    for infile in args.input_files:
      #output.append(infile)
      tmp_fda=fda.readline()
      tmp_fda=tmp_fda.split(" ")
      dir=tmp_fda[3]+tmp_fda[4]+tmp_fda[5]
      for iso in args.isovalue:
        infile_temp=infile.split('.')
        output.append(infile_temp[0])
        output.append(dir)
        output.append(iso)
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-gradient',\
                  'testData/'+infile_temp[0]+'neg.cgrad.nrrd', iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #repostion
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #lindstrom
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-lindstrom','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow_conflict
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-allow_conflict','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #centroid_conflcit
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-centroid_conflict','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp_far
        sp.call(['./isodual3D','-trimesh','-multi_isov','-s','-position',pos,'-o', 'out.off','-clamp_far','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #reselectg
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-reselectg','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow conflict + repos
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-allow_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp conflict + repos
        sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-clamp_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        
        output.append('break')
      print infile,'done'
    outlist.append(output)
  return outlist

'''
def create_iso_command_neg(args):
  output=[]
  isovalue=['-10.1','-10.2']
  for infile in args.input_files:
    #output.append(infile)
    for iso in isovalue:
      infile_temp=infile.split('.')
      output.append(infile_temp[0]+'neg')
      output.append(iso)
      
      sp.call(['./isodual3D', '-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off',\
                iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      
      #repostion
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-reposition',\
                iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      #lindstrom
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-lindstrom',\
               iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      #allow_conflict
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-allow_conflict',\
                iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      #centroid_conflcit
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-centroid_conflict',\
                iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      #clamp_far
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-clamp_far',\
               iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      #reselectg
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-reselectg',\
                iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      #allow conflict + repos
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-reposition',\
                '-allow_conflict',iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      #clamp conflict + repos
      sp.call(['./isodual3D','-multi_isov','-trimesh','-s','-position','gradCD','-o', 'out.off','-reposition',\
                '-clamp_conflict',iso, 'testData/'+infile_temp[0]+'neg.nrrd'])
      sp.call(['./findedge','140','out.off'])
      op=sp.check_output(['./findEdgeCount','-fp','out.line'])
      opsplit=op.split()
      output.append(opsplit[1])
      
      output.append('break')
  return output
 '''       