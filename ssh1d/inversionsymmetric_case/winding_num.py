#!/usr/bin/env python
from scipy import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt





if __name__=='__main__':
	kl=linspace(0,pi,100)
	t=1.0
	for dlt in [-0.2,0.2]:
	  hk_x=[];hk_y=[]
	  hk_xnml=[];hk_ynml=[]
	  for k in kl:
	      hkx=2*t*cos(k)**2+2*dlt*sin(k)**2
	      hky=(t-dlt)*sin(2*k)
	      norm=sqrt(hkx**2+hky**2)
	      hk_x.append(hkx)
	      hk_y.append(hky)
	      hk_xnml.append(hkx/norm)
	      hk_ynml.append(hky/norm)
	  # use masked arrays
	  hk_x=array(hk_x)
	  hk_y=array(hk_y)
	  hk_xnml=array(hk_xnml)
	  hk_ynml=array(hk_ynml)
	  x1 = np.ma.masked_array(hk_x[:-1], np.diff(hk_x)>=0)
	  x2 = np.ma.masked_array(hk_x[:-1], np.diff(hk_x)<=0)
	  figure(1)
	  plt.plot(hk_x,hk_y,'r-',label=r'$\delta t=$'+str(dlt))
	  plt.plot(x1,hk_y[:-1],'r<')
	  plt.plot(x2,hk_y[:-1],'r>')
	  grid()
	  axvline(x=0)
	  axhline(y=0)
	  xlabel(r"$\sigma_x$")
	  ylabel(r"$\sigma_y$")
	  legend(loc=0)
	  savefig("wdnum.png")
	  figure(2)
	  plot(hk_xnml,hk_ynml,lw=2,label=r'$\delta t=$'+str(dlt))
	  grid()
	  axvline(x=0)
	  axhline(y=0)
	  xlabel(r"$\sigma_x$")
	  ylabel(r"$\sigma_y$")
	  legend(loc=0)
	  savefig("wdnum_normalized.png")
	show()
