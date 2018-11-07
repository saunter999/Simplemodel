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

def sceqs(mu,m,n,T,U):
	nup=0.;ndn=0.;
	for ek in ekls:
	    Eupk=ek-U*m/2.0-mu+U*n/2.
	    Ednk=ek+U*m/2.0-mu+U*n/2.
	    nup+=fermi(T,Eupk)
	    ndn+=fermi(T,Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	mp=nup-ndn
	mup=optimize.brentq(func,-20.,20.,args=(mp,n,T,U))
	return array([mup,mp])
	
def func(mu,m,n,T,U):
	nup=0.;ndn=0.;
	for ek in ekls:
	    Eupk=ek-U*m/2.0-mu+U*n/2.
	    Ednk=ek+U*m/2.0-mu+U*n/2.
	    nup+=fermi(T,Eupk)
	    ndn+=fermi(T,Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	return nup+ndn-n
	
def DetermineTc(mu,Tc,n,U)
	nup=0.;ndn=0.;
	for ek in ekls:
	    Eupk=ek-mu+U*n/2.
	    Ednk=ek-mu+U*n/2.
	    nup+=fermi(Tc,Eupk)
	    ndn+=fermi(Tc,Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	mup=optimize.brentq(func,-20.,20.,args=(0.,n,Tc,U))
	Tcp=optimize.brentq(Tcfunc,0.,1.,args=(mup,n,U))
	return array([mup,Tcp]) 

def Tcfunc(Tc,mu,n,U)
	for ek in ekls:
	    obk=1./cosh( (ek-mu+U*n/2.) / (2.0*Tc) )**2 
	    ob+=obk
	ob=(ob/Nk**2)*U/(4.0*Tc)
	return ob-1.0
		
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

	NT=20
	Nk=40
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
	#nls=[0.65,0.7,0.75,1.25,1.3,1.35];
	nls=[1.0]

	Dlt_mu=1e-4;Dlt_m=1e-4;Nit=200;
	for U in Uls:
	  print "U=",U
	  for n in nls:
	      print "---------------------"
	      print "n=",n
	      fmag=open("mag_neq"+str(n)+".txt",'w')
	      fcheck=open("OUT.dat",'w')
	      print>>fmag,"#occupation n=",n,"U=",U
	      print>>fmag,"#T,|m|,mu,iteration times"
	      print>>fcheck,"#iter#,","mu-muprime,","m-mprime"
	      mu=0.;m=1.;
	      for T in Tmesh:
		  print "T=",T
	      	  print>>fcheck,"#occupation n=",n,"U=",U,"T=",T
		  it=0
		  while True:
		     (mup,mp)=sceqs(mu,m,n,T,U)
		     it+=1
		     if (abs(mup-mu)<Dlt_mu and abs(mp-m)<Dlt_m) or (it==Nit):
		        print>>fmag,T,mp,mup,it
		        print>>fcheck,it,mu-mup,m-mp
			break
		     mu=(mu+mup)/2.0
		     m=(m+mp)/2.0
		mu=0.0;Tc=0.01;Dlt_T=1e-2;
	        it=0;
		while True:
		     (mup,Tcp)=DetermineTc(mu,Tc,n,U)
		     it+=1
		     if (abs(mup-mu)<Dlt_mu and abs(Tcp-Tc)<Dlt_m) or (it==Nit):
			print Tc,it
			break
		     mu=(mu+mup)/2.0
		     Tc=(Tc+Tcp)/2.0
