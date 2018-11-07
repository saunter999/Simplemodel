#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
import numpy as np


def set_ham(kx,Nj):
	ham=zeros((Nj,Nj),dtype=complex)
	for j in range(Nj):
		ham[j,j]=-2.*t*cos(kx-alpha*j)
		if j==0:
		    ham[0,1]=-t
		if j==(Nj-1):
		    ham[Nj-1,Nj-2]=-t
		if j!=0 and j!=(Nj-1):
		    ham[j,j-1]=-t
		    ham[j,j+1]=-t 
	return ham
      
def Ham_solver(H):
	(Eig,v)=LA.eig(H)
	idx=Eig.argsort() #This one is needed as the eigenvaules returned by LA.eig aren't in the order of -/+ levels.
	Eig=Eig[idx]
	v=v[:,idx]
	return (Eig,v)

if __name__=='__main__':
	### Model parameter ###
	t=1.0;nflux=2./5.;   ##nflux->flux quantum 
	alpha=2.*pi*nflux
        ##########################
	itask1="1"
	if itask1=="1":
            print "--------------------------------------------------"
	    print "task1:calculating spectrum with open boundary in y direction:running"
	    Nk=100;kp=linspace(-pi,pi,Nk)
	    Nj=151;  #number of chains in y direction
	    print "model parameter:"
	    print "t=",t,",flux quantum",nflux,",alpha=",alpha
	    print "Nk=",Nk
	    print "Nj=",Nj
	    Eigkx=zeros((Nj,Nk))
	    for indk,kx in enumerate(kp):
		ham=zeros((Nj,Nj),dtype=complex)
		ham=set_ham(kx,Nj)
		print "indk=",indk,",","kx=",kx
		#print ham
		Eig,v=Ham_solver(ham)
		###we need to transpose v to make that v[i,:] is the eigenvector of eigvalue w[i]. 
		v=v.transpose()
		for j in range(Nj):
		    Eigkx[j,indk]=Eig[j].real
	    figure(1)
	    plot(kp,Eigkx[0,:],"*")
	    plot(kp,Eigkx[1,:],"8")
	    plot(kp,Eigkx[2,:],"^")
	    plot(kp,Eigkx[3,:],"o")
	    #for j in range(Nj):
	     #   plot(kp,Eigkx[j,:])
	    grid()
	else:
            print "--------------------------------------------------"
	    print "task1:calculating spectrum with open boundary in y direction:NOT running"
	show()
