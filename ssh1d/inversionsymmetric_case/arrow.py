#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

y = np.linspace(0,100,100)
x = np.cos(y/5.)

# use masked arrays
x1 = np.ma.masked_array(x[:-1], np.diff(x)>=0)
x2 = np.ma.masked_array(x[:-1], np.diff(x)<=0)

# print the line and the markers in seperate steps
plt.plot(x, y, 'k-')
plt.plot(x1, y[:-1], 'k<')
plt.plot(x2, y[:-1], 'k>')
plt.show()
