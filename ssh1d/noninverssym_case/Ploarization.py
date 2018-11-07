#!/usr/bin/env python
from scipy import *
from pylab import *
import numpy as np
from numpy import linalg as LA
from scipy.integrate import quad
from mpl_toolkits.mplot3d import Axes3D

if __name__=='__main__':
	"""
	  Evaulation of the Berry phase of ssh model using the Bloch wavefunction of the Bloch Hamiltonian(or the Kernal of Hamiltonian).
	  For dlt<0,the berry phase is pi;dlt>0,the berry phase is zero.
	  so the former is topological nontrival while the latter is trival.
	  This is in accord with the edge state calculation.
	"""
	t=1
	Nk=100
	Npara=20
	kl=linspace(0,pi,Nk)	
	dk=pi/Nk
	Hamk=zeros((2,2),dtype=complex)
	uk_ocp=zeros((2,Nk),dtype=complex)  #used to store the eigenvector of the Bloch Hamiltonian at each k point. 
	uk_unocp=zeros((2,Nk),dtype=complex)
	Pol_ocp=zeros((Npara,Npara)) 
	mls=linspace(-t,t,Npara)
	dltls=linspace(-t,t,Npara)
	f=open('Polarz.dat','w')
	print >>f,"#m","t","Polarzation"
	for indm,m in enumerate(mls):
	  for indt,dlt in enumerate(dltls):
	    for idk,k in enumerate(kl):
	       #set the Bloch Hamiltonian.
	       Hamk[0,0]=m
	       Hamk[0,1]=2.*t*cos(k)**2+2.*dlt*sin(k)**2+1j*(dlt-t)*sin(2.0*k)
	       Hamk[1,0]=conjugate(Hamk[0,1])
	       Hamk[1,1]=-Hamk[0,0]
	       Eig,v=LA.eig(Hamk)
	       idx=Eig.argsort() 
	       Eig=Eig[idx]
	       v=v[:,idx]
	       uk_ocp[:,idk]=v[:,0]
	       uk_unocp[:,idk]=v[:,1]
#Using the uk evaulated numerically and the equivalent form introduced by D.Vanderbilt to evaluate it,which is numerically stable.
	    umul1=1.0
	    umul2=1.0
	    for i in range(1,Nk):
		uovlap1=np.vdot(uk_ocp[:,i-1],uk_ocp[:,i])
		uovlap2=np.vdot(uk_unocp[:,i-1],uk_unocp[:,i])
		umul1=umul1*uovlap1	
		umul2=umul2*uovlap2	
	    print >>f,m,dlt,(np.log(umul1)).imag/(2*pi) 
	f.close()
	path='./Polarz.dat'
	dat=loadtxt(path).transpose()
	Polar=dat[2]
	dim=int(sqrt(len(Polar)))
	mapPolar=zeros((dim,dim))
	for j in range(dim):
	    for i in range(dim):
		mapPolar[i,j]=Polar[i+j*dim]
	imshow(mapPolar,extent=[-t,t,-t,t],cmap=cm.jet,origin='lower')
	colorbar(shrink=1.0)
	xlabel(r'$\frac{m}{t}$',size=15)
	ylabel(r'$\frac{\delta t}{t}$',size=15,rotation=0)
	title('map of ' r'$ Polarization$')
	savefig("Polarz.png")
	show()
