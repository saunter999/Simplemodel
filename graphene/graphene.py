#!/usr/bin/env python
from scipy import *
from numpy import linalg as LA
from pylab import *
from scipy import integrate
from scipy import optimize
#from mpl_toolkits.mplot3d import Axes3D

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

def muequ(mu,ne):
	wmu=[]
	Awupmu=[]
	Awdnmu=[]
	for idx,w in enumerate(wmesh):
	   if w<mu:
		wmu.append(w)
		Awupmu.append(Awup[idx])
		Awdnmu.append(Awdn[idx])
	wmu=array(wmu);Awupmu=array(Awupmu);Awdnmu=array(Awdnmu)
	Occ=2.0*(integrate.simps(Awupmu,wmu)+integrate.simps(Awdnmu,wmu))
	return Occ-ne
		
		

if __name__=='__main__':
    """
	This code is implmenation of the band structure and dos calculation of 
	NN Tight binding model of graphene.Using the dos,we calcuate the chemical
	potential for a given "ne" of the unit cell.
	Q.H 2015.11.07
    """
    Nk=30
    t=1.0   
    ##Define three NN vectors.
    e1=array([0,1])
    e2=array([-sqrt(3)/2.0,-0.5])
    e3=array([sqrt(3)/2.0,-0.5])
    evct=array([e1,e2,e3])
    ##------------ Plotting the bands ----------------------------##
    ##Specifying the position of these special K points in the 1st BZ.
    Gam=array([0.,0.])
    K=array([4.*pi/(3*sqrt(3)),0.]) 
    Kprm=array([2.*pi/(3*sqrt(3)),2.*pi/3.])
    M=array([0.,2.*pi/3.])
    hsymK=[Gam,K,Kprm,M]
    Dir=[K-Gam,Kprm-K,M-Kprm,Gam-M]
    Eupk=[];Ednk=[];
    for idx,kp in enumerate(hsymK):
      for i in arange(Nk):
 	k=kp+Dir[idx]*float(i)/Nk	
	(Ek,v)=Banddisp(k[0],k[1])
	Eupk.append(Ek.real)
	Ednk.append(-Ek.real)
    Eupk=array(Eupk);Ednk=array(Ednk)
    knum=arange(Eupk.shape[0])
    plotsg=1
    if (plotsg==1):
       ytext=-4.5
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
       plot(knum,Ednk,label="- band")
       plot(knum,Eupk,label="+ band")
       legend(loc=0)
       grid()
       savefig("graphenebandstruct.png")
    
    ##-------------Evaluating the dos --------------------------------##
#    figure(2)
    Nk=1200
    Ne=300
    dlt=0.01
    G1=array([2*pi*sqrt(3.)/3.0,-2*pi/3.0])
    G2=array([0.,4*pi/3.0])
    G1mesh=[] 
    G2mesh=[]
    for i in arange(Nk):
	G1mesh.append( (float(i)/Nk)*G1 )
	G2mesh.append( (float(i)/Nk)*G2 )
    G1mesh=array(G1mesh)
    G2mesh=array(G2mesh)
    Eupk=[]
    Ednk=[]
    for k1 in G1mesh:
	for k2 in G2mesh:
		k=k1+k2 ##k lies in the 1st BZ defined by G1 and G2.
		(Eig,v)=Banddisp(k[0],k[1])
		Eupk.append(Eig[0].real)
		Ednk.append(Eig[1].real)
    Eupk=array(Eupk)
    Ednk=array(Ednk)
    wmesh=linspace(-4.*t,4.*t,Ne)
    wmesh=(list(wmesh))
    wmesh.append(0)
    wmesh.sort()
    wmesh=array(wmesh)
    Awup=[];Awdn=[]
    fdos=open("dosgraphene.dat",'w')
    for w in wmesh:
        Guploc=0.0j
        Gdnloc=0.0j
	for E in Eupk:
	    Guploc+=1./(w-E+1j*dlt)
	for E in Ednk:
	    Gdnloc+=+1./(w-E+1j*dlt)
	Guploc=Guploc/Nk**2
	Gdnloc=Gdnloc/Nk**2
	Aup=-1.0/pi*Guploc.imag
	Adn=-1.0/pi*Gdnloc.imag
	Awup.append(Aup)
	Awdn.append(Adn)
    Awup=array(Awup)
    Awdn=array(Awdn)
    I1=integrate.simps(Awup,wmesh)  #Integration over the wmesh,which is used to normalize the dos.
    I2=integrate.simps(Awdn,wmesh)
    print "I1=",I1,"I2=",I2
    Awup=Awup/sqrt(I1)
    Awdn=Awdn/sqrt(I2)
    for i in arange(Ne+1):
	    print>>fdos,wmesh[i],Awup[i],Awdn[i],Awup[i]+Awdn[i]
    fdos.close()
   ###---To submit the job,we need to comment out the plot---### 
   # plot(wmesh,Awup)    
   # plot(wmesh,Awdn)
    #plot(wmesh,Awup+Awdn)
   # xlabel(r"$\omega$")
   # ylabel("DOS")
   # xlim([-4.0*t,4.0*t])
   # grid()
   # savefig("dosgraphene.png")
   # show()
   ###-----------------------###
    ##------------Determining mu for a given ne----------##
    nemesh=linspace(1.0,3.0,20)
    nemesh=list(nemesh)
    nemesh.append(2.0)
    nemesh.sort()
    nemesh=array(nemesh)
    mune=[]
    fmu=open("occ_mu.dat","w")
    for ne in nemesh:
    	mu=optimize.brentq(muequ,-3.89*t,4.0*t,args=(ne,))
	mune.append(mu)
	print>>fmu,ne,mu
    mune=array(mune) 
    fmu.close()
