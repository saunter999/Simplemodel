#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA


def Hinc_matrix(k):
	Hinc[0,0]=-delta-2.0*tprime*cos(phi-sqrt(3)*k)
	Hinc[1,1]=delta-2.0*tprime*cos(phi+sqrt(3)*k)
	Hinc[0,1]=-2.0*t*cos(sqrt(3)/2.0*k)
	Hinc[1,0]=conjugate(Hinc[0,1])

def HaccR_matrix(k):
	HaccR[0,0]=-2.0*tprime*cos(phi+sqrt(3)/2.0*k)
	HaccR[1,1]=-2.0*tprime*cos(phi-sqrt(3)/2.0*k)
	HaccR[0,1]=1.0
	HaccR[1,0]=0.0

def HaccL_matrix(k):
	HaccL[0,0]=-2.0*tprime*cos(phi+sqrt(3)/2.0*k)
	HaccL[1,1]=-2.0*tprime*cos(phi-sqrt(3)/2.0*k)
	HaccL[0,1]=0.0
	HaccL[1,0]=1.0

if __name__=='__main__':
        N=100
	kp=linspace(-pi,pi,N)
	Nj=100                    #number of chains in y direction
	ham=zeros((2*Nj,2*Nj),dtype=complex)   #kernal of the Hamiltonian
        kxedge=zeros(2*Nj)
	Hinc=zeros((2,2),dtype=complex)        #Hinc:H_in_chain hamiltonian
	HaccR=zeros((2,2),dtype=complex)        #Hacc:H_ac_chain hamiltonian
	HaccL=zeros((2,2),dtype=complex)        #Hacc:H_ac_chain hamiltonian	
	delta=0.
	t=1.0
	tprime=0.1
	phi=pi/2.0
	Eigk=zeros((N,2*Nj))
	for indk,kx in enumerate(kp):
	    Hinc_matrix(kx)
	    HaccR_matrix(kx)
	    HaccL_matrix(kx)
	    for i in range(Nj):
		ham[2*i,2*i]=Hinc[0,0]
		ham[2*i,2*i+1]=Hinc[0,1]
		ham[2*i+1,2*i]=Hinc[1,0]
		ham[2*i+1,2*i+1]=Hinc[1,1]
                if (i !=0) and (i !=(Nj-1)):
			ham[2*i,2*i-2]=HaccL[0,0]
			ham[2*i,2*i-1]=HaccL[0,1]
			ham[2*i+1,2*i-2]=HaccL[1,0]
			ham[2*i+1,2*i-1]=HaccL[1,1]
			ham[2*i,2*i+2]=HaccR[0,0]
			ham[2*i,2*i+3]=HaccR[0,1]
			ham[2*i+1,2*i+2]=HaccR[1,0]
			ham[2*i+1,2*i+3]=HaccR[1,1]
		if i==0:
			ham[0,2]=HaccR[0,0]
			ham[0,3]=HaccR[0,1]
			ham[1,2]=HaccR[1,0]
			ham[1,3]=HaccR[1,1]
		if i==Nj-1:
			ham[2*i,2*i-2]=HaccL[0,0]	
			ham[2*i,2*i-1]=HaccL[0,1]	
			ham[2*i+1,2*i-2]=HaccL[1,0]	
			ham[2*i+1,2*i-1]=HaccL[1,1]	
	    Eig,v=LA.eig(ham)		
	    idx=Eig.argsort()
	    Eig=Eig[idx]
	    v=v[:,idx]
	    for i in range(2*Nj):
		Eigk[indk,i]=Eig[i].real
	#for i in range(2*Nj):
	for i in linspace(Nj,Nj+10):
 		plot(kp,Eigk[:,i]) 
        #The plotting above is nontrival:we are plotting 2*Nj lines-- 
        #--each of which is a function of k points.The i is to index--
        #--the eigvaule,which the : in Eigk[:,i] is the k points.
	grid()
	savefig("edge_nnn.png")
	show()
