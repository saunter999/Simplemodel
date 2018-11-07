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
	dltls=[0.0,0.1*t,0.2*t,0.3*t]
	klist=linspace(-pi,pi,Nk)
	for dlt in dltls:
	    Ek=[]
	    for k in klist:
	        Ek.append(2*sqrt( (t*cos(k))**2+(dlt*sin(k))**2 ))	
            Ek=array(Ek)
	    plot(klist,Ek,label=r'$\delta t=$'+str(dlt)) 
	    plot(klist,-Ek,label=r'$\delta t=$'+str(dlt))
	legend(bbox_to_anchor=(0.85, 1), loc=2, borderaxespad=0.0)
	savefig("bandstru.png")
	show()
