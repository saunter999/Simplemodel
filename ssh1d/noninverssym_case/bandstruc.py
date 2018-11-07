#!/usr/bin/env python
"""Plotting the band structure of ssh in 1d.
   Parameter:lattice constant(A-B distance)->a=1
	     hopping: t->1
	     perturbation strength: dlt
""" 
from scipy import *
from pylab import *

if __name__=="__main__":
	a=1
	t=1
	Nk=500
	mls=[0.0,0.2*t]
	dltls=[0.0,0.2*t]
	klist=linspace(-pi,pi,Nk)
	for dlt in dltls:
	  for m in mls:
	      Ek=[]
	      for k in klist:
		  Ek.append(2*sqrt(m**2/4.0+(t*cos(k))**2+(dlt*sin(k))**2 ))	
	      Ek=array(Ek)
	      plot(klist,Ek,label=r'$\frac{\delta t}{t}=$'+str(dlt)+','+r'$\frac{m}{t}=$'+str(m)) 
	      plot(klist,-Ek,label=r'$\frac{\delta t}{t}=$'+str(dlt)+','+r'$\frac{m}{t}=$'+str(m))
	legend(bbox_to_anchor=(0.85, 1), loc=2, borderaxespad=0.0,labelspacing=0.5,handlelength=1,prop={'size':14})
	savefig("bandstru.png")
	show()
