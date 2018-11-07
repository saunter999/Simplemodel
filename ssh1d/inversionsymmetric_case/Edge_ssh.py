#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA




if __name__=='__main__':
	"""
	Edge states in 1d ssh model with open boundary condition:
	1)dlt<0,it is topological nontrival,there is a zero mode,which is localized at the two boundaries
	(by inspecting the wavefucntion distribution);
	2)dlt>0,it is topological trival.   
	"""
	N=50            #number of unit cell.
	t=1.0 
	dlt=-0.1        
	ham=zeros((2*N,2*N))
	x=zeros(2*N)
	for i in range(2*N):
		if(i==0):
		   ham[i,i+1]=t+dlt
		if(i==2*N-1):
		   ham[i,i-1]=t+dlt
		if(i>0) and (i<2*N-1):	
		  if i%2==1:
		     ham[i,i-1]=t+dlt
		     ham[i,i+1]=t-dlt
		  if i%2==0:
		     ham[i,i-1]=t-dlt
		     ham[i,i+1]=t+dlt
	Eig,v=LA.eig(ham)
	idx=Eig.argsort()
	Eig=Eig[idx]
	v=v[:,idx]
	figure(1)
	scatter(x,Eig,color='r',s=10)
	ylabel("E",rotation=0)
	savefig('spectrum_ssh_obc.png')
	#grid()
	figure(2)
	for idx,E in enumerate(Eig):
	    if abs(E)<abs(dlt):
	      plot(range(2*N),v[:,idx]**2,'r-',marker='s')
	      savefig('localized_edge.png')
	xlabel("n:lat_coordinate")
	ylabel(r"$|\psi_n|^2$",rotation=0)
	show()
