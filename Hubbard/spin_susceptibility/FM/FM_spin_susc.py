#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize

def dispersion(kls):
	Ek=[]
	for kx in kls:
	    for ky in kls:
	        Ek.append(-2.*t*(cos(kx)+cos(ky))-4.*tp*cos(kx)*cos(ky) )	
	Ek=array(Ek)
	return Ek

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


def scseqs(mu,m,n,T,U):
	nup=0.0;ndn=0.0
	for ek in Ekls:
	    Eupk=ek-mu+U*n/2.-U*m
	    Ednk=ek-mu+U*n/2.+U*m
	    nup+=fermi(T,Eupk)
	    ndn+=fermi(T,Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	mp=0.5*(nup-ndn)
	mup=optimize.brentq(mufunc,-20.,20.,args=(m,n,T,U))
	return array([mup,mp])
	    
def mufunc(mu,m,n,T,U):
	nup=0.;ndn=0.;
	for ek in Ekls:
	    Eupk=ek-mu+U*n/2.-U*m
	    Ednk=ek-mu+U*n/2.+U*m
	    nup+=fermi(T,Eupk)
	    ndn+=fermi(T,Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	return nup+ndn-n
	

def spin_susc(mu,m,n,T):
	chiup=0.0
	chidn=0.0
	for ek in Ekls:
	   Eupk=ek-mu+U*n/2.-U*m
	   Ednk=ek-mu+U*n/2.+U*m
	   chiup+=fermi_derivative(T,Eupk)
	   chidn+=fermi_derivative(T,Ednk)
	chiup=-1.0*chiup/Nk**2
	chidn=-1.0*chidn/Nk**2
	return [0.25*(chiup+chidn)/(1.0-0.5*U*(chiup+chidn) ),chiup*chidn/((chiup+chidn)-2.*U*chiup*chidn)]


def fermi_derivative(T,e):
	if abs(e/T)>10.0:return 0.0
	return -1.0/(4.*T* ( cosh(0.5*e/T) )**2)	

if __name__=="__main__":
	t=1.0;tp=-0.45;U=3;n=0.2
	Nk=80;NT=50
	kls=linspace(-pi,pi,Nk)
	Tls=linspace(1e-4,0.4,NT)
	Ekls=dispersion(kls)

	Dlt_mu=1e-4;Dlt_m=1e-4;Nit=200;
	fcheck=open("scsiterations.dat",'w')
	print>>fcheck,"#iter#,","mu-muprime,","m-mprime"
	fmag=open("mag_neq"+str(n)+".txt",'w')
	print>>fmag,"#occupation n=",n,"U=",U
	print>>fmag,"#T,|Sz|,mu,iteration times"
	fchi=open("chi"+str(n)+".txt",'w')
	print>>fchi,"#T,chi_para,chi_n"

	for T in Tls:
	    print "Temperature:",T
	    mu=0.;m=0.5
	    it=0
	    while True:
		(mup,mp)=scseqs(mu,m,n,T,U)	
		it+=1
		if (abs(mup-mu)<Dlt_mu and abs(mp-m)<Dlt_m) or (it==Nit):
		    print>>fmag,T,mp,mup,it
		    print>>fcheck,it,mu-mup,m-mp
		    break
		mu=(mu+mup)/2.0
		m=(m+mp)/2.0
	    chi_para,chi_n=spin_susc(mu,m,n,T)
	    print>>fchi,T,chi_para,chi_n
