'''
read data form file
'''
import matplotlib.pyplot as P
from matplotlib import collections, axes, transforms
from matplotlib.colors import colorConverter
import numpy as N


# Function to read and return  data from file
def readDatafrom (fname):
	# read the data into a file
	f = open(fname, 'r')
	list_data=[]
	# read data and append to list_data
	for line in f:
		list_data.append(line.split())
	return  list_data


