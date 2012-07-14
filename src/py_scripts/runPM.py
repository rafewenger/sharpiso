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

oname=''.join([sys.argv[1].split('.')[0],'-dc.ply'])
onameDCF=''.join([sys.argv[1].split('.')[0],'-dc.dcf'])

onameOFF= ''.join([sys.argv[1].split('.')[0],'-dc.off'])

onameOFFTRI= ''.join([sys.argv[1].split('.')[0],'-dc-tri.off'])

onameLine= ''.join([sys.argv[1].split('.')[0],'-dc-tri.line'])




wine=['wine','PolyMender-dc',sys.argv[1],octtree_depth, scale ,oname]
sp.call (wine)

winedcf=['wine','PolyMender-dc',sys.argv[1],octtree_depth, scale,onameDCF]
sp.call (winedcf)


meshconv=['meshconv',oname,'-c','off']
sp.call(meshconv)
tri=['ijkmesh','triangulate','-uniform',onameOFF, onameOFFTRI]
sp.call(tri)
findedge=['findedge','140',onameOFFTRI]
sp.call(findedge)
show=['meshview',onameOFF,onameLine]
sp.call(show)

print 'fname',sys.argv[1]

print 'oname',oname

print 'onameOFF',onameOFF

print 'onameOFFTRI',onameOFFTRI

print 'onameLine',onameLine
