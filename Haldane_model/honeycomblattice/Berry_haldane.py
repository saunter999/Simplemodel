#!/usr/bin/env python
"""
Generate the band structure of Haldane model with magnetic flux on square lattice.
"""
from scipy import *
from numpy import linalg as LA
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def Eig_solver(qx,qy):
      Hk=zeros((2,2),dtype=complex)
      kvec=array([qx,qy])
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
      v=v[:,idx]
      return (Eig,v)



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
    f=open("taskcont.dat","r")
    lines=f.readlines()
    for l in lines:
        taskcont=l.split()
        itask1=taskcont[0]
        itask2=taskcont[1]
        itask3=taskcont[2]
        itask4=taskcont[3]
        GrapheneSwitch=taskcont[4]
    if GrapheneSwitch=="G":
      print "--------Graphene Switch on------"
      t=-1.0;phi=0.0;tprime=0.0*t;mos=0.*t
      print "t=",t,"phi=",phi,"tprime=",tprime,"mos=",mos
    else:
      print "--------Haldane Model parameters----------"
      t=-1.0;phi=pi/4.0;tprime=0.0*t;mos=5.*t
      print "t=",t,"phi=",phi,"tprime=",tprime,"mos=",mos
    e1=array([0,1])
    e2=array([-sqrt(3)/2.0,-0.5])
    e3=array([sqrt(3)/2.0,-0.5])
    evct=array([e1,e2,e3])
    v1=array([sqrt(3.0),0.0])
    v2=array([-sqrt(3)/2.0,1.5])
    v3=array([-sqrt(3)/2.0,-1.5])
    vvct=array([v1,v2,v3])
    
   #1)Plotting the Bandstructure across the 2d square BZ
    Nk=50;kp=linspace(-pi,pi,Nk)
    Eigband=zeros((2,Nk,Nk))
    qx=kp;qy=kp;qx,qy=meshgrid(qx,qy)
    if(itask1=="1"):
      print "---------------------------------------------------------"
      print "task1--Plotting Bandstruc across 2d square BZ:running"
      for ind1,q1 in enumerate(kp):
	  for ind2,q2 in enumerate(kp):
		  (Eig,v)=Eig_solver(q1,q2)
		  for i in range(len(Eig)):
		      Eigband[i,ind1,ind2]=Eig[i].real
      ax = plt.axes(projection='3d')
      figure(1)
      for i in range(len(Eig)):
	  ax.plot_surface(qx, qy, Eigband[i], rstride=1,cstride=1,linewidth=0,cmap=cm.coolwarm, antialiased=False)
      savefig("bands_BZ.png")
    else:
      print "---------------------------------------------------------"
      print "task1--Plotting Bandstruc across 2d square BZ:NOT running"
	
    #2):plot bandstructure along high symmetry k lines
    #In the case of honeycomb lattice,the kpath generated below is the high symmetry lines 
    #along the conventional BZ,not the Wignersize BZ of honeycomb lattice,so we don't see the gappless
    #points at the K and Kprime point
    if(itask2=="1"):
      print "---------------------------------------------------------"
      print "task2--Plotting Bandstruc along high symmetry lines in 1st BZ:running"
      Gam=array([0.,0.])
      K=array([4.*pi/(3*sqrt(3)),0.]) 
      Kprm=array([2.*pi/(3*sqrt(3)),2.*pi/3.])
      M=array([0.,2.*pi/3.])
      hsymK=[Gam,K,Kprm,M]
      Dir=[K-Gam,Kprm-K,M-Kprm,Gam-M]
      Eupk=[];Ednk=[];
      for idK,kp in enumerate(hsymK):
         for i in arange(Nk):
 	     k=kp+Dir[idK]*float(i)/Nk	
	     (Eig,v)=Eig_solver(k[0],k[1])
	     idx=Eig.argsort()
	     Eig=Eig[idx]
	     v=v[:,idx]
	     Ednk.append(Eig[0].real)
	     Eupk.append(Eig[1].real)
      Eupk=array(Eupk);Ednk=array(Ednk)
      knum=arange(Eupk.shape[0])
      figure(2)
      plot(knum,Ednk,label="- band")
      plot(knum,Eupk,label="+ band")
      ytext=-3.3
      axvline(x=0*Eupk.shape[0]/4.0,c='k')
      text(0*Eupk.shape[0]/4.0-2, ytext, r"$\Gamma$", fontsize=20)
      axvline(x=1*Eupk.shape[0]/4.0,c='k')
      text(1*Eupk.shape[0]/4.0-2, ytext, r"$K$", fontsize=20)
      axvline(x=2*Eupk.shape[0]/4.0,c='k')
      text(2*Eupk.shape[0]/4.0-2, ytext, r"$K^\prime$", fontsize=20)
      axvline(x=3.*Eupk.shape[0]/4.0,c='k')
      text(3*Eupk.shape[0]/4.0-2, ytext, r"$M$", fontsize=20)
      axvline(x=4*Eupk.shape[0]/4.0,c='k')
      text(4*Eupk.shape[0]/4.0-2, ytext, r"$\Gamma$", fontsize=20)
      legend(loc=0)
      grid()
      savefig("bands_highsymline.png")
    else:
      print "---------------------------------------------------------"
      print "task2--Plotting Bandstruc along high symmetry lines in 1st BZ:NOT running"
	
    #3) Evaluation of Berry phase as a function of k1=(0->1) as a reduced coord of one of two k unit vectors.
       #i.e. we integrate over along a closed path in k2 direction
      
       ##TO DO...still some problem(discontinuity of berry phase at middle point of k1--due to the Dirac point?)
    if(itask3=="1"):
      Nk1=80;Nk2=280;
      print "---------------------------------------------------------"
      print "task3--Evaluation of Berry phase as a function of k1:running"
      G1=array([2*pi*sqrt(3.)/3.0,-2*pi/3.0])
      G2=array([0.,4*pi/3.0])
      G1mesh=[] 
      G2mesh=[]
      for i in arange(Nk1):
     	G1mesh.append( (float(i)/Nk1)*G1 )
      for i in arange(Nk2):
	G2mesh.append( (float(i)/Nk2)*G2 )
      G1mesh=array(G1mesh)
      G2mesh=array(G2mesh)
      ukdn=zeros((2,Nk2),dtype=complex)
      ukup=zeros((2,Nk2),dtype=complex)
      Bphup=zeros(Nk1)
      Bphdn=zeros(Nk1)
      for idk1,k1 in enumerate(G1mesh):
	  for idk2,k2 in enumerate(G2mesh):
	      k=k1+k2 
	      (Eig,v)=Eig_solver(k[0],k[1])
	      idx=Eig.argsort() #This one is needed as the eigenvaules returned by LA.eig aren't in the order of -/+ levels.
	      Eig=Eig[idx]
	      v=v[:,idx]
	      ukdn[:,idk2]=v[:,0]
	      ukup[:,idk2]=v[:,1]
          umuldn=1.;umulup=1.
          for i in range(Nk2):
	     if i==0:
	       uprddn=np.vdot(ukdn[:,Nk2-1],ukdn[:,0])
	       uprdup=np.vdot(ukup[:,Nk2-1],ukup[:,0])
	     else:
	       uprddn=np.vdot(ukdn[:,i-1],ukdn[:,i])
	       uprdup=np.vdot(ukup[:,i-1],ukup[:,i])
	     umuldn=umuldn*uprddn	
	     umulup=umulup*uprdup	
	  Bphdn[idk1]=(np.log(umuldn)).imag
	  Bphup[idk1]=(np.log(umulup)).imag
      k1list=linspace(0.,1.,Nk1)
    #  for idx,bp in enumerate(Bphdn):
#	if bp<0.1:
#	   Bphdn[idx]=bp+2.*pi
 #     for idx,bp in enumerate(Bphup):
#	if bp>0.1:
#	   Bphup[idx]=bp-2.*pi
      plot(k1list,(Bphdn/pi)%2,'*',label="dn")
      plot(k1list,(Bphup/pi)%2,'o',label="up")
      #plot(k1list,((Bphdn+Bphup)/pi)%2,'s',label="TOT")
      xlabel("k1 dirction(0->1)")
      ylabel(r"$\phi_{BR}$")
      legend(loc=0)
      grid()
      savefig("Berry(Haldane_model).png")
      print "winding number=",Bphdn[Nk1-1]-Bphdn[0],Bphup[Nk1-1]-Bphup[0]
      ###modified berry phase(ensuring continuity of berry phase...)
      #for bp in Bphdn:
	
    else:
      print "---------------------------------------------------------"
      print "task3--Evaluation of Berry phase as a function of k1:NOT running"

    #4)Evaluating the Berry phase around "K=array([4.*pi/(3*sqrt(3)),0.])" 
    #checked for the pi Berry phase around K for graphene model Hamiltonian.
    if(itask4=="1"):
      print "---------------------------------------------------------"
      print "task4--Evaulating Berry phase around K=array([4.*pi/(3*sqrt(3)),0.]):running"
      K=array([4.*pi/(3*sqrt(3)),0.]) 
      Nagl=10;rad=0.1 #radius of the circle
      Theta=linspace(0,2*pi,Nagl)
      Kcir=[]
      ukdn=zeros((2,Nagl),dtype=complex)
      ukup=zeros((2,Nagl),dtype=complex)
      for agl in Theta:
	k=K+array([rad*cos(agl),rad*sin(agl)])
	Kcir.append(k)
      Kcir=array(Kcir)
      for idK,k in enumerate(Kcir):
	      (Eig,v)=Eig_solver(k[0],k[1])
	      idx=Eig.argsort() #This one is needed as the eigenvaules returned by LA.eig aren't in the order of -/+ levels.
	      Eig=Eig[idx]
	      v=v[:,idx]
	      ukdn[:,idK]=v[:,0]
	      ukup[:,idK]=v[:,1]
      umuldn=1.;umulup=1.
      for i in range(Nagl):
	     if i==0:
	       uprddn=np.vdot(ukdn[:,Nagl-1],ukdn[:,0])
	       uprdup=np.vdot(ukup[:,Nagl-1],ukup[:,0])
	     else:
	       uprddn=np.vdot(ukdn[:,i-1],ukdn[:,i])
	       uprdup=np.vdot(ukup[:,i-1],ukup[:,i])
	     umuldn=umuldn*uprddn	
	     umulup=umulup*uprdup	
      print "Berry phase of dn band around K is:",(np.log(umuldn)).imag
      print "Berry phase of up band around K is:",(np.log(umulup)).imag
	
    else:
      print "---------------------------------------------------------"
      print "task4--Evaulating Berry phase around K=array([4.*pi/(3*sqrt(3)),0.]):NOT running"
    show()
