#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA

def Banddisp(qx,qy):
      #Hk will be set to be the kernal of the Hamiltonian#
      Hk=zeros((2,2),dtype=complex)
      kvec=array([qx,qy])
      Hk[0,0]=0.0
      Hk[1,1]=0.0
      for ev in evct:
	Hk[0,1]+=exp( 1.0j*dot(ev,kvec) )
      Hk[0,1]=-t*Hk[0,1]
      Hk[1,0]=conjugate(Hk[0,1])
      Eig,v=LA.eig(Hk)
      return (Eig,v)

if __name__=="__main__":
     t=1.0;lamd=0.13 
     ##Define three NN vectors.
     e1=array([0,1])
     e2=array([-sqrt(3)/2.0,-0.5])
     e3=array([sqrt(3)/2.0,-0.5])
     evct=array([e1,e2,e3])
     ##----Plotting Fermi surface for some ne(or mu)-----##
     alpha=80.
     mu=-0.2
     Nk=100
     FSwtdn=zeros((Nk,Nk))
     #FSwtup=zeros((Nk,Nk))
     kp=linspace(-pi,pi,Nk)
     for idx,kx in enumerate(kp):
	  for idy,ky in enumerate(kp):
	      (Ek,v)=Banddisp(kx,ky)
	      FSwtdn[idx,idy]=exp(-alpha*(-Ek[0].real-mu)**2)
	      #FSwtup[idx,idy]=exp(-alpha*(banddisp(kx,ky).real-mu)**2)
     title(r"$\mu=-0.2$")
     imshow(FSwtdn,origin='lower', extent=[-pi, pi, -pi, pi],interpolation='bicubic')
     #imshow(FSwtup,interpolation='bicubic')
     savefig("FSgraphene.png")
     exit()
     ##----Drawing Brillouion Zone------##
     K0=array([4.*pi/(3*sqrt(3)),0.])
     K1=array([2.*pi/(3*sqrt(3)),2.*pi/3.])
     K2=array([-2.*pi/(3*sqrt(3)),2.*pi/3.])
     K3=array([-4.*pi/(3*sqrt(3)),0.])
     K4=array([-2.*pi/(3*sqrt(3)),-2.*pi/3.])
     K5=array([2.*pi/(3*sqrt(3)),-2.*pi/3.])
     KBZvtx=[K0,K1,K2,K3,K4,K5]
     Edge=[K1-K0,K2-K1,K3-K2,K4-K3,K5-K4,K0-K5]
     KBZ=[]
     for idk,kp in enumerate(KBZvtx):	
	 for i in arange(Nk):
  	    k=kp+Edge[idk]*float(i)/Nk
	    KBZ.append(k)
     KBZ=array(KBZ)	
     KBZx=[];KBZy=[]
     for k in KBZ:
	KBZx.append(k[0])
	KBZy.append(k[1])
     KBZx=array(KBZx)
     KBZy=array(KBZy)
     #plot(KBZx,KBZy,"r--")
     show()
