# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 13:20:42 2015

@author: mhumbert
"""

import numpy as np
#import pylab
import matplotlib.pyplot as plt
temp = np.loadtxt('temp.csv')
time = np.loadtxt('time.csv')
vacf = np.loadtxt('vacf.csv')
radius = np.loadtxt('radius.csv')
g = np.loadtxt('pairdistrobution.csv')

plt.figure()
plt.plot(radius,g)
plt.show()