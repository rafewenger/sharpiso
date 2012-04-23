'''
FindEdgeCount implementation in Python
Reads all the .line files in the folder and runs findEdgeCount on them

'''
import read_data_module
import plot_results
import sys

import sys
import os
import shlex
import subprocess
import glob
print sys.argv
graph=0
show_pts=0
for  opts in sys.argv:
	if opts == '-help':
		print '  in help '
	if opts =='-graph':
		graph=1
		print ' plot draw_graph'
	if 	opts =='-show_pts':
		show_pts=1


graph_data=[]
labels=[]
for files in glob.glob('*.line'):

	fname=files
	#read the line file
	data=read_data_module.readDatafrom (fname)
	#print 'data 0 ',data[0]
	if data[0]!=['LINEC']:
		print 'This is not a LINEC file.'
		sys.exit()

	#num of points and lines
	num_points=int(data[2][0])
	num_line=int(data[2][1])

	# array of vertex degrees
	deg=[0]*num_points

	#print deg
	#find Degrees of each point
	for l in data[3+num_points+1:]:
		deg[int(l[0])]+=1
		deg[int(l[1])]+=1

	#print 'deg',deg
	vert_deg0=[]
	vert_deg1=[]
	vert_deg2=[]
	vert_deg3=[]
	vert_deg_grt_thn_3=[]
	for i in range(num_points):
		if deg[i]==0:
			vert_deg0.append(i)
		if deg[i]==1:
			vert_deg1.append(i)
		if deg[i]==2:
			vert_deg2.append(i)
		if deg[i]==3:
			vert_deg3.append(i)
		if deg[i]>3:
			vert_deg_grt_thn_3.append(i)
	if show_pts==1:
		print fname
		print 'vert_deg1',vert_deg1
		print 'vert_deg3',vert_deg3
		print 'vert_deg>3',vert_deg_grt_thn_3

	#print to console
	if int(graph)==0 and int(show_pts)==0:
		print 'filename',fname
		print 'vertices with deg 0',len(vert_deg0)
		print 'vertices with 1',len(vert_deg1)
		print 'vertices with 2',len(vert_deg2)
		print 'vertices with 3',len(vert_deg3)
		print 'vertices with >3',len(vert_deg_grt_thn_3)
		print 'vert with deg one three or more ',len(vert_deg1)\
		+len(vert_deg3)+len(vert_deg_grt_thn_3)

	if (int(graph)==1):
		out_all=fname.split('.')
		out_fname=out_all[0]+'.'+out_all[2]
		labels.append(out_fname)
		graph_data.append(len(vert_deg_grt_thn_3)\
		+len(vert_deg3)+len(vert_deg1))

#print 'graph_data',graph_data
if graph==1:
	plot_results. plot_direct(graph_data,labels)


