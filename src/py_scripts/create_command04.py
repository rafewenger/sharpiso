import subprocess as sp

def create_iso_command(args):
  outlist=[]
  #read the data sheet info 
  fda=open ('testData/data_sheet.txt','r')
  print '**************ISODUAL3D with perfect gradients**************'
  for pos in args.position:
    output=[]
    print 'position ',pos
    for infile in args.input_files:
      #output.append(infile)
      tmp_fda=fda.readline()
      tmp_fda=tmp_fda.split(" ")
      if len(tmp_fda)>=5:
        dir=tmp_fda[3]+tmp_fda[4]+tmp_fda[5]
      else:
        dir='-1'+'-1'+'-1'
      
      for iso in args.isovalue:
        infile_temp=infile.split('.')
        output.append(infile_temp[0])
        output.append(dir)
        output.append(iso)
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #repostion
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #lindstrom
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-lindstrom', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow_conflict
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-allow_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #centroid_conflcit
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-centroid_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp_far
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-clamp_far', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #reselectg
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reselectg', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow conflict + repos
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-allow_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp conflict + repos
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-clamp_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        output.append('break')
      print '--',infile,'done'
    outlist.append(output)
  return outlist


############### run with cgradient
def create_iso_command_cgrad(args):
  outlist=[]
  #read the data sheet info 
  fda=open ('testData/data_sheet.txt','r')
  print '**************ISODUAL3D with central gradients**************'
  for pos in args.position:
    output=[]
    print 'position cgradient',pos
    for infile in args.input_files:
      #output.append(infile)
      tmp_fda=fda.readline()
      tmp_fda=tmp_fda.split(" ")
      if len(tmp_fda)>=5:
        dir=tmp_fda[3]+tmp_fda[4]+tmp_fda[5]
      else:
        dir='-1'+'-1'+'-1'
      
      for iso in args.isovalue:
        infile_temp=infile.split('.')
        output.append(infile_temp[0])
        output.append(dir)
        output.append(iso)
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #repostion
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #lindstrom
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-lindstrom','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow_conflict
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-allow_conflict','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #centroid_conflcit
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-centroid_conflict','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp_far
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-clamp_far','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #reselectg
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reselectg','-gradient',\
                 'testData/'+infile_temp[0]+'.cgrad.nrrd', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #allow conflict + repos
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-allow_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        #clamp conflict + repos
        sp.call(['./isodual3D','-trimesh','-s','-position',pos,'-o', 'out.off','-reposition','-clamp_conflict', iso, 'testData/'+infile])
        sp.call(['./findedge','140','out.off'])
        op=sp.check_output(['./findEdgeCount','-fp','out.line'])
        opsplit=op.split()
        output.append(opsplit[1])
        
        output.append('break')
      print '--',infile,'done'
    outlist.append(output)
  return outlist


     