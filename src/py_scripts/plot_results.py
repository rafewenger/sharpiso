'''
	script to read data and plot it
'''
#!/usr/bin/env python
import read_data_module
import matplotlib.pyplot as plt
import numpy

def plot_from_file(fname):
	print 'fname',fname
	# read data from the file
	data=read_data_module.readDatafrom (fname)
	print 'data',data
	# select a column
	col=6
	col_data=[]
	for f in data:
		col_data.append(f[col])
	print 'col_data',col_data

	# FIND THE LABELS
	LABELS=[]
	for f in data:
		LABELS.append(f[0])
	print 'LABELS',LABELS

	# plot a BAR GRAPH
	x_len=len(col_data)
	ind=numpy.arange(x_len)
	#convert from string to int
	col_data = map(float, col_data)
	WIDTH=0.35
	rects = plt.bar(ind, col_data, WIDTH)
	plt.xticks(ind+WIDTH/2., LABELS)
	plt.ylabel('y label')
	plt.xlabel('x label')
	plt.show()

def plot_direct(data,LABELS):
	# plot a BAR GRAPH
	x_len=len(data)
	
	#convert from string to int
	data = map(float, data)
	WIDTH=0.35
	print 'ind',ind
	rects = plt.bar(ind, data, WIDTH)
	plt.xticks(ind,LABELS)
	plt.ylabel('y label')
	plt.xlabel('x label')
	plt.show()
