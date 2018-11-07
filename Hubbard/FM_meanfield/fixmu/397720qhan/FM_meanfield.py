#!/usr/bin/env python
from scipy import *
from scipy import optimize
from pylab import *

def DOS():
    fdos=open("dos.txt","w")
    Emesh=linspace(-8.0*t,8.0*t,Ne)
    Aw=[]
    for omg in Emesh:
	Gloc=0.0j
	for ek in ekls:
	   Gloc+=1./(omg-ek+1.0j*dlt)
	Gloc=Gloc/Nk**2
	Aw.append(-1./pi*Gloc.imag)
	print>>fdos,omg,-1./pi*Gloc.imag
    #figure(1)
    #plot(Emesh,Aw,label='dos')
   

def fermi(T,E):
	"""
	 fermi function:when |x| in exp(x) is very large,it
	 cause overflow problem in python.So we use the following
	 method.
	"""
	if (E/T<-10.0):
	    return 1.0
	if (E/T>10.0):
	    return 0.0
	return 1.0/(exp(E/T)+1.0)

def sceqs(n,m,mu,T,U):
	nup=0.;ndn=0.;
	for ek in ekls:
	    Eupk=ek-U*m/2.0-mu+U*n/2.
	    Ednk=ek+U*m/2.0-mu+U*n/2.
	    nup+=fermi(T,Eupk)
	    ndn+=fermi(T,Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	mp=nup-ndn
	np=nup+ndn
	return array([np,mp])
	
	
if __name__=='__main__':
	"""
	This code is implementing the self-consistent solution
	of the possible ferromagnetic phase in the Hubbard model 
	using mean field.
	We take the case of 2D square lattice as an example,which
	has the band dispersion as ek=-2t( cos(kx)+cos(ky) ).We examine the temperature dependence of 
	the magnetization.
	Treatment of other dispersions and filling factor is straigtforward.
        Q.H 2015.09.24
	"""

	NT=40
	Nk=100
	t=1.0
	NU=1
	#Uls=linspace(0.,10.,NU)
	Uls=[6.0];
	dlt=0.01
	Ne=100
	Tmesh=linspace(1e-2,2.,NT)  #Temperature
	kmesh=linspace(-pi,pi,Nk)  #kmesh
	ekls=[]
	for kx in kmesh:
	    for ky in kmesh:
		ekls.append(-2.0*t*(cos(kx)+cos(ky)))   #Band dispersion
	ekls=array(ekls)
	#DOS()
	muls=arange(1.2,5.0,0.2);

	Dlt_n=1e-4;Dlt_m=1e-4;Nit=200;
	for U in Uls:
	  print "U=",U
	  for mu in muls:
	      print "---------------------"
	      print "mu=",mu
	      fmag=open("mag_mueq"+str(mu)+".txt",'w')
	      fcheck=open("OUT.dat",'w')
	      print>>fmag,"#chemical potential mu=",mu,"U=",U
	      print>>fmag,"#T,|m|,n,iteration times"
	      print>>fcheck,"#iter#,","n-nprime,","m-mprime"
	      n=0.;m=1.;
	      for T in Tmesh:
		  print "T=",T
	      	  print>>fcheck,"#chemical potential mu=",mu,"U=",U,"T=",T
		  it=0
		  while True:
		     (np,mp)=sceqs(n,m,mu,T,U)
		     it+=1
		     if (abs(np-n)<Dlt_n and abs(mp-m)<Dlt_m) or (it==Nit):
		        print>>fmag,T,mp,np,it
		        print>>fcheck,it,n-np,m-mp
			break
		     n=(n+np)/2.0
		     m=(m+mp)/2.0
