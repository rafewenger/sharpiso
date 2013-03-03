#!/usr/bin/python 
import sys
import subprocess as sp


'''
the first argument is the name of the file {has to be a .stl}
the second argument is the the depth of the octree 
the third argument is the scale

'''

octtree_depth = sys.argv[2]
scale = sys.argv[3]

oname=''.join([sys.argv[1].split('.')[0],'-dc-',sys.argv[2],'-',sys.argv[3],'.ply'])
onameDCF=''.join([sys.argv[1].split('.')[0],'-dc-',sys.argv[2],'-',sys.argv[3],'.dcf'])

onameOFF= ''.join([sys.argv[1].split('.')[0],'-dc-',sys.argv[2],'-',sys.argv[3],'.off'])
onameOFF1= ''.join([sys.argv[1].split('.')[0],'-dc-',sys.argv[2],'-',sys.argv[3],'-clean.off'])
onameOFFTRI= ''.join([sys.argv[1].split('.')[0],'-dc-tri-',sys.argv[2],'-',sys.argv[3],'.off'])

onameLine= ''.join([sys.argv[1].split('.')[0],'-dc-tri.line'])




wine=['wine','PolyMender-dc',sys.argv[1],octtree_depth, scale ,oname]
sp.call (wine)
print 'PolyMender-dc ran successfully'
winedcf=['wine','PolyMender-dc',sys.argv[1],octtree_depth, scale, onameDCF]
sp.call (winedcf)
print 'PolyMender-dc to create .dcf ran successfully'
'''
meshconv=['meshconv',oname,'-c','off']
'''
meshconv =['python','readPLY.py',oname, onameOFF]
sp.call(meshconv)

print 'readPLY ran successfully'
meshclean =['/usr/local/bin/ijkmesh','deldup_polyvert',onameOFF, onameOFF1]
sp.call(meshclean)

meshclean =['/usr/local/bin/ijkmesh','deldup_polyvert',onameOFF1, onameOFF1]
sp.call(meshclean)
print 'clean mesh'
tri=['ijkmesh','triangulate','-split_max_angle',onameOFF1, onameOFFTRI]
sp.call(tri)

print 'triangulation done'

findedge=['findsharp','140',onameOFFTRI]
sp.call(findedge)
#show=['meshview',onameOFF,onameLine]
#sp.call(show)
print 'findedge done'
print 'fname',sys.argv[1]

print 'oname',oname

print 'onameOFF',onameOFF

print 'onameOFFTRI',onameOFFTRI

print 'onameLine',onameLine
