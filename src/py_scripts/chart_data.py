'''
functions to print out charts 
'''

import numpy as np
import matplotlib.pyplot as plt

from pylab import *


def plot_data(data):
  names=[]
  y=[]
  for f in data:
    fval = f.split()
    if len(fval) > 1:
      names.append(fval[0])
      y.append(fval[len(fval)-3])
  pos = np.arange(len(names))
  #debug 
  print 'names',names
  print 'y', y 
  y = map(float, y)
  
  val = 3+10*rand(5)    # the bar lengths
  pos = np.arange(5)+.5    # the bar centers on the y axis
  ind=np.arange(len(names))
  plt.barh(pos, val, align='center') 
  plt.yticks(ind,names)
  plt.show()
  

