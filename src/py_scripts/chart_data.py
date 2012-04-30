'''
functions to print out charts 
'''

import numpy as np
import matplotlib.pyplot as plt




def plot_data(data):
  names=[]
  y=[]
  for f in data:
    fval = f.split()
    if len(fval) > 1:
      names.append(fval[0])
      y.append(fval[len(fval)-3])
  print names
  print y
  pos = np.arange(len(names))
  y = map(float, y)
  ind=np.arange(len(names))
  plt.barh(pos, y, align='center') 
  plt.yticks(ind,names)
  plt.show()
  

