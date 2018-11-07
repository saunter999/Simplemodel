#!/usr/bin/env python
"""
Chern_nunmber calculation of Haldane model
"""
from scipy import *
from numpy import linalg as LA
from pylab import *
import numpy as np

def Momentum2D(Nk):
      G1=array([2*pi*sqrt(3.)/3.0,-2*pi/3.0])
      G2=array([0.,4*pi/3.0])
      klist=zeros((Nk,Nk,2))
      for i in range(Nk):
	  for j in range(Nk):
	      klist[i,j,:]=(i*G1+j*G2)/(Nk+0.0)
      return klist



def wf2D(klist):
      Nk=len(klist)
      wf=zeros((Nk,Nk,2),dtype=complex) ## 2-component spinor since it is a two-band model
      for i in range(Nk):
	for j in range(Nk):
            Hk=zeros((2,2),dtype=complex)
            kvec=klist[i,j,:]
	    #----assignment of the diagonal terms of Hamiltonian kernal h(k):the complex NNN hopping term---#
	    for vv in vvct:
	      Hk[0,0]+=cos(phi-dot(vv,kvec))
	      Hk[1,1]+=cos(phi+dot(vv,kvec))
	    Hk[0,0]=-2.0*tprime*Hk[0,0]+mos
	    Hk[1,1]=-2.0*tprime*Hk[1,1]-mos
	    #----assignment of the off-diagonal terms of Hamiltonian kernal h(k):the NN hopping term---#
	    for ev in evct:
	      Hk[0,1]+=exp(1j*dot(ev,kvec))
	    Hk[0,1]=-t*Hk[0,1]
	    Hk[1,0]=conjugate(Hk[0,1])
	    Eig,v=LA.eig(Hk)
	    #----sorting the eigvaule and its corresponding eigvector in the order of increasing magnitude of eigvaule.---#
	    idx=Eig.argsort()
	    Eig=Eig[idx]
	    v=v[:,idx] ## column vector stores the eigenvector
	    wf[i,j,:]=v[:,0]
      return wf



def chern_number(wf):
    Nk=len(wf)
    fs=zeros((Nk,Nk),dtype=complex)
    for i in range(Nk):
	ir=(i+1) % Nk
	for j in range(Nk):
	    ju=(j+1) % Nk
	    inprd=np.vdot(wf[i,j,:],wf[ir,j,:])
	    U1x=inprd/abs(inprd) 
	    inprd=np.vdot(wf[ir,j,:],wf[ir,ju,:])
	    U2y=inprd/abs(inprd) 
	    inprd=np.vdot(wf[ir,ju,:],wf[i,ju,:])
	    U3x=inprd/abs(inprd) 
	    inprd=np.vdot(wf[i,ju,:],wf[i,j,:])
	    U4y=inprd/abs(inprd) 
	    fs[i,j]=log(U1x*U2y*U3x*U4y)
    return fs
	    
if __name__=='__main__':
    """
    Parameters specification:
    Nk->Number of k points in one direction.
    t -> nearest hopping.
    tprime ->next-nearest hopping 
    phi -> phase factor affliated with tprime due to the magnetic flux
    evct ->three vectors connecting nearest neighbours of the Honeycomb lattice
    vvct ->three vectors connecting next nearest neighbours of Honeycomb lattice
    """
    "--------Haldane Model parameters----------"
    t=1.0;tprime=1.0*t;
#    print "t=",t,"phi=",phi,"tprime=",tprime,"mos=",mos
    e1=array([0,1])
    e2=array([-sqrt(3)/2.0,-0.5])
    e3=array([sqrt(3)/2.0,-0.5])
    evct=array([e1,e2,e3])   ##NN bonds
    v1=array([sqrt(3.0),0.0])
    v2=array([-sqrt(3)/2.0,1.5])
    v3=array([-sqrt(3)/2.0,-1.5])   
    vvct=array([v1,v2,v3])  ##NNN bonds
    Nk=20 
    klist=Momentum2D(Nk)

    Nph=50;Nm=100
    phils=linspace(-pi,pi,Nph);
    mb=6.
    mosls=linspace(-mb,mb,Nm)*tprime
    cn=zeros((Nph,Nm))
    fout=open("chernum.txt",'w')
    for iph,phi in enumerate(phils):
	for im,mos in enumerate(mosls):
	  print "m=",mos,'phi=',phi
	  wf=wf2D(klist)
	  fs=chern_number(wf)
	  for i in range(Nk):
	    for j in range(Nk):
	        cn[iph,im]+=fs[i,j].imag 
	  cn[iph,im]=cn[iph,im]/(2*pi)
	  print>>fout,phi,mos,cn[iph,im]

    imshow(cn.transpose(), interpolation='nearest',origin='lower',extent=[-pi,pi,-mb,mb],aspect='equal')
    colorbar()
    plot(phils,3*sqrt(3)*sin(phils),'k--',lw=2)
    plot(phils,-3*sqrt(3)*sin(phils),'k--',lw=2)
    tight_layout()
    savefig("chernum.png")
    fout.close()
   
    show()
