#!/usr/bin/env python
from scipy import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D






if __name__=='__main__':
	kl=linspace(0,pi,50)
	t=1.0;dlt=-0.2;m=0.1
	hk_x=[];hk_y=[];hk_z=[]
	for k in kl:
	    hk_x.append(2*t*cos(k)**2+2*dlt*sin(k)**2)
	    hk_y.append((t-dlt)*sin(2*k))
	    hk_z.append(m)
	# use masked arrays
	hk_x=array(hk_x)
	hk_y=array(hk_y)
	hk_z=array(hk_z)
	x1 = np.ma.masked_array(hk_x[:-1], np.diff(hk_x)>=0)
	x2 = np.ma.masked_array(hk_x[:-1], np.diff(hk_x)<=0)
	fig = plt.figure(1)
	ax = fig.gca(projection='3d')
	arrowr = u'$\u2192$'
	arrowl = u'$\u2190$'
	ax.plot(hk_x,hk_y,hk_z,'r-')
#	ax.plot(x1,hk_y[:-1],hk_z[:-1],'r',marker=arrowr,markersize=20)
	#ax.plot(x2,hk_y[:-1],hk_z[:-1],'r',marker=arrowl,markersize=20)
	grid()
	#axvline(x=0)
	#axhline(y=0)
	ax.set_xlabel(r"$\sigma_x$")
	ax.set_ylabel(r"$\sigma_y$")
	ax.set_zlabel(r"$\sigma_z$")
	savefig("windingnum.png")
	show()
