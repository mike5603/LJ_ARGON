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
g = np.loadtxt('pairdistribution.csv')

plt.figure()
plt.plot(radius,g)
plt.xlabel('radius')
plt.ylabel('g(r)')
plt.title('Pair Distribution Function')
plt.show()

plt.figure()
plt.plot(time,vacf)
plt.xlabel('time')
plt.ylabel('VACF')
plt.title('Velocity Autocorrelation Function')
plt.show()