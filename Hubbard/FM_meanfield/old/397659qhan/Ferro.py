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

def selffuns(x,n,T,U):
	"""
	x[0]--mu
	x[1]--m
	"""
	nup=0.;ndn=0.;
	for ek in ekls:
	    Eupk=ek-U*x[1]/2.0-x[0]+U*n/2.
	    Ednk=ek+U*x[1]/2.0-x[0]+U*n/2.
	    nup+=fermi(T,Eupk)
	    ndn+=fermi(T,Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	return [nup+ndn-n,nup-ndn-x[1]]
	
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

	NT=80
	Nk=100
	t=1.0
	NU=50
	Uls=linspace(0.,10.,NU)
	Uls=[5.0];
	dlt=0.01
	Ne=100
	Tmesh=linspace(1e-2,3.,NT)  #Temperature
	kmesh=linspace(-pi,pi,Nk)  #kmesh
	  #Hubbard U.It is strange that U has to be so high(as 10) to have a magnetization at low temperature.
	        #because the dos of paramagntic state at half filling is diverging(according to Mahan,Uc=0...)
	ekls=[]
	for kx in kmesh:
	    for ky in kmesh:
		ekls.append(-2.0*t*(cos(kx)+cos(ky)))   #Band dispersion
	ekls=array(ekls)
	#DOS()
	nls=[0.9,0.94,0.98,1.0,1.02,1.06,1.1];
#	for T in Tmesh:
#	    print T
#	    for n in nls:
#	      print "---------------------"
#	      fmag=open("mag_neq"+str(n)+".txt",'w')
#	      print>>fmag,"#occupation n=",n,"T=",T
#	      print>>fmag,"#U,|m|,mu,sucess"
#	      for U in Uls:
#		print "U=",U
#		sol = optimize.root(selffuns,[n*U/2., 10.0],args=(n,T,U),method="anderson") ##(mu,m)
#		print>>fmag,U,abs(sol.x[1]),sol.x[0],sol.success
#		if sol.success==0:break;
#		print sol.success
#	      fmag.close()
	for U in Uls:
	  print "U=",U
	  for n in nls:
	      print "---------------------"
	      print "n=",n
	      fmag=open("mag_neq"+str(n)+".txt",'w')
	      print>>fmag,"#occupation n=",n,"U=",U
	      print>>fmag,"#T,|m|,mu,sucess"
	      for T in Tmesh:
	        print "T=",T
		sol = optimize.root(selffuns,[n*U/2., 10.0],args=(n,T,U),method="anderson") ##(mu,m)
		print>>fmag,T,abs(sol.x[1]),sol.x[0],sol.success
		if sol.success==0:break;
		print sol.success
