#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
from scipy import integrate
from scipy import optimize

def banddisp(qx,qy):
	kvec=array([qx,qy])
	fk=0.0j;mk=0.
	for ev in NNvct:
	    fk+=exp(-1.0j*dot(ev,kvec))
	mk=mk+sin(dot(v1,kvec))-sin(dot(v2,kvec))+sin(dot(v3,kvec))
	return sqrt(t**2*fk*conjugate(fk)+4.*lamd**2*mk**2)

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
		
if __name__=="__main__":
   Nk=30
   t=1.0;lamd=0.13 
   ##Define three NN vectors.
   e1=array([0,1])
   e2=array([-sqrt(3)/2.0,-0.5])
   e3=array([sqrt(3)/2.0,-0.5])
   NNvct=array([e1,e2,e3])
   ##Define three NNN vectors
   v1=e3-e2;v2=e1-e2;v3=e1-e3;
   NNNvct=array([v1,v2,v3])
   ##-------Plotting the band structure along gamma-K-K'-M-gamma-----#
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
	Ek=banddisp(k[0],k[1]).real
	Eupk.append(Ek)
	Ednk.append(-Ek)
   Eupk=array(Eupk);Ednk=array(Ednk)
   knum=arange(Eupk.shape[0])
   plotsg=0
   if (plotsg==1):
     title(r"$\lambda$="+str(lamd))
     ytext=-3.5
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
     axhline(y=3.*sqrt(3)*lamd,c='k',ls='--')
     axhline(y=-3.*sqrt(3)*lamd,c='k',ls='--')
     text(4*Eupk.shape[0]/4.0-1.,3.*sqrt(3)*lamd, r"$y=\frac{E_g(K)}{2}$", fontsize=16)
     plot(knum,Ednk,label="- band")
     plot(knum,Eupk,label="+ band")
     legend(loc=0)
     grid()
     savefig("bandstruct.png")

  
   ###----density of states----###
   Nk=1200
   Ne=200
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
	      Eigk=banddisp(k[0],k[1]).real
	      Eupk.append(Eigk)
	      Ednk.append(-Eigk)
   Eupk=array(Eupk)
   Ednk=array(Ednk)
   wmesh=linspace(-4.*t,4.*t,Ne)
   Awup=[];Awdn=[]
   fdos=open("doskanemele.dat",'w')
   print>>fdos,"#omega,Awup,Awdn,Awtot"
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
   print "I1_upband=",I1,"I2_dnband=",I2
   Awup=Awup/sqrt(I1)
   Awdn=Awdn/sqrt(I2)
   for i in arange(Ne):
	  print>>fdos,wmesh[i],Awup[i],Awdn[i],Awup[i]+Awdn[i]
   fdos.close()

   ##---Determining mu for a given ne--##
   nemesh=linspace(0.1,3.99,200)
   nemesh=list(nemesh)
   nemesh.append(2.0)
   nemesh.sort()
   nemesh=array(nemesh)
   mune=[]
   fmu=open("occ_mu.dat","w")
   print>>fmu,"#mu,ne_occ"
   for ne in nemesh:
    	mu=optimize.brentq(muequ,-3.99*t,4.0*t,args=(ne,))
	mune.append(mu)
	print>>fmu,mu,ne
   mune=array(mune) 
   fmu.close()
