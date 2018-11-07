#!/usr/bin/env python
from scipy import *
from pylab import *

def EigenE():
	print "--------EigenEnergy---------"
	clf();
	U=2.0;
	muls=linspace(-1.9*U,1.9*U,100)
	E0=[];E1=[];E3=[];
	for mu in muls:
	    E0.append(0.0)
	    E1.append(-mu)
	    E3.append(U-2.*mu)
	plot(muls,E0,label="holon")
	plot(muls,E1,label="eletron")
	plot(muls,E3,label="doublon")
	legend(loc=0)
	xlabel("mu",size=19)
	ylabel("E",size=19)
	grid()
	savefig("EigenE.pdf")

def occ():
	print "--------occupancy---------"
	clf();
	U=3.5;T=0.1;
	beta=1./T
	print "%8s %8s %8s" % ('U','T','beta')  
	print "%8.4f %8.4f %8.4f" % (U,T,beta) 
	muls=linspace(-1.9*U,1.9*U,100)
	nls=[]
	b=exp(-beta*U)
	for mu in muls:
	    if mu<0:
		a=exp(beta*mu)
	        nave=2.*(a+a**2*b)/(1.+2.*a+a**2*b) 
	    else:
		a=exp(-beta*mu)
		nave=2.*(a+b)/(a**2+2.*a+b)
	    nls.append(nave)
	plot(muls,nls,"*",c='g')
	xlabel("mu",size=19)
	ylabel("n_occ",size=19)
	axhline(y=1,ls='--',c='b')
	axvline(x=U/2.,ls='--',c='b')
	axvline(x=U/4.,ls='-',c='b')
	axvline(x=3.*U/4.,ls='-',c='b')
	savefig("occ.pdf")
	
	
def Energy():
	clf();
	U=0.5;T=0.1;
	beta=1./T
	print "--------Energy---------"
	print "%8s %8s %8s" % ('U','T','beta')  
	print "%8.4f %8.4f %8.4f" % (U,T,beta) 
	muls=linspace(-1.9*U,1.9*U,100)
	E=[]
	b=exp(-beta*U)
	for mu in muls:
	    if mu<0:
		a=exp(beta*mu)
	        Eave=( (U-2.*mu)*a**2*b-2.0*mu*a)/(1.+2.*a+a**2*b) 
	    else:
		a=exp(-beta*mu)
		Eave=((U-2.*mu)*b-2.0*mu*a)/(a**2+2.*a+b)
	    E.append(Eave)
	plot(muls,E,"*",c='g')
	xlabel("mu",size=19)
	ylabel("E_ave",size=19)
	grid()
	axvline(x=U/2.,ls='--',c='b')
	savefig("Energy.pdf")

def doubleocc():
	clf();
	U=0.5;T=1.;
	beta=1./T
	print "--------double occ---------"
	print "%8s %8s %8s" % ('U','T','beta')  
	print "%8.4f %8.4f %8.4f" % (U,T,beta) 
	muls=linspace(-3.9*U,3.9*U,100)
	dbocc=[]
	b=exp(-beta*U)
	for mu in muls:
	    if mu<0:
		a=exp(beta*mu)
	        db=(a**2*b)/(1.+2.*a+a**2*b) 
	    else:
		a=exp(-beta*mu)
		db=b/(a**2+2.*a+b)
	    dbocc.append(db)
	plot(muls,dbocc,"*",c='g')
	xlabel("mu",size=19)
	ylabel("Double_occ",size=19)
	grid()
	axvline(x=U/2.,ls='--',c='b')
	axhline(y=0.25,ls='--',c='b')
	savefig("double_occ.pdf")

def magnet_mom():
	clf();
	U=0.9;T=0.1;
	beta=1./T
	print "--------magnetic moment---------"
	print "%8s %8s %8s" % ('U','T','beta')  
	print "%8.4f %8.4f %8.4f" % (U,T,beta) 
	muls=linspace(-5.9*U,5.9*U,100)
	#muls=linspace(-2.,2.,100)
	m2ls=[]
	b=exp(-beta*U)
	for mu in muls:
	    if mu<0:
		a=exp(beta*mu)
	        m2=(2.*a)/(1.+2.*a+a**2*b) 
	    else:
		a=exp(-beta*mu)
		m2=2.*a/(a**2+2.*a+b)
	    m2ls.append(m2)
	plot(muls,m2ls,"*",c='g')
	xlabel("mu",size=19)
	ylabel("m^2",size=19)
	grid()
	axvline(x=U/2.,ls='--',c='b')
	#axhline(y=0.25,ls='--',c='b')
	savefig("magnet_mom.pdf")

if __name__=="__main__":

	EigenE();
	occ();
	Energy();
	doubleocc();
	magnet_mom();
	show()
