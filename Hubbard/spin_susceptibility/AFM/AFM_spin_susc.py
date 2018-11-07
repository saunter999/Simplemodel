#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize

def magneticBz(kls):
	mbz=[]
	for kx in kls:
	    for ky in kls:
		if (abs(kx)+abs(ky))<=pi:
		    mbz.append([kx,ky])
	return mbz

def ksiekpm_eval(mbz): 
	ksiEkp=[];ksiEkm=[]
	for k in mbz:
	    kx=k[0];ky=k[1]
	    ek=-2.*t*(cos(kx)+cos(ky))-4.*tp*cos(kx)*cos(ky)
	    ekq=-2.*t*(cos(kx+pi)+cos(ky+pi))-4.*tp*cos(kx+pi)*cos(ky+pi)
	    ksiEkp.append(0.5*(ek+ekq))
	    ksiEkm.append(0.5*(ek-ekq))
	ksiEkp=array(ksiEkp)
	ksiEkm=array(ksiEkm)
	return (ksiEkp,ksiEkm)

def Dk_eval(E,mQ):
	 return sqrt( E**2+(U*mQ)**2 )
	   
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

def scseqs(mu,mQ,n,T,U):
	mQp=0.0
	for idx,e in enumerate(ksiEkp):
	    Dk=Dk_eval(ksiEkm[idx],mQ) 
	    Eupk=e+Dk-mu 
	    Ednk=e-Dk-mu 
	    mQp+=(fermi(T,Eupk)-fermi(T,Ednk))*(U*mQ)/Dk
	mQp=-1.0*mQp/Nk**2;
	mup=optimize.brentq(mufunc,-20.,20.,args=(mQ,n,T,U)) 
	return array([mup,mQp])

def mufunc(mu,mQ,n,T,U):
	nup=0.;ndn=0.;
	for idx,e in enumerate(ksiEkp):
	    Dk=Dk_eval(ksiEkm[idx],mQ) 
	    Eupk=e+Dk-mu 
	    Ednk=e-Dk-mu 
	    nup+=fermi(T,Eupk)
	    ndn+=fermi(T,Ednk)
	nup=2.0*nup/Nk**2;
	ndn=2.0*ndn/Nk**2;
	return nup+ndn-n

def spin_susc(mu,mQ,n,T):
	chibare=0.0
	a11=0.0
	a12=0.0
	for idx,e in enumerate(ksiEkp):
	    Dk=Dk_eval(ksiEkm[idx],mQ) 
	    Eupk=e+Dk-mu 
	    Ednk=e-Dk-mu 
	    Ak=(fermi(T,Eupk)-fermi(T,Ednk))*ksiEkm[idx]**2/Dk**3
	    fermi_deri_sum=fermi_derivative(T,Eupk)+fermi_derivative(T,Ednk)
	    fermi_deri_min=fermi_derivative(T,Eupk)-fermi_derivative(T,Ednk)
	    Bk=fermi_deri_sum*(U*mQ/Dk)**2
	    chibare+=Ak+Bk
	    a11+=fermi_deri_sum
	    a12+=fermi_deri_min*U*mQ/Dk
	chibare=-1.0*chibare/Nk**2
	a11=-2.0*a11/Nk**2
	a12=2.0*a12/Nk**2
	a21=0.5*a12
	a22=chibare
	chi_nbare=(a11*a22-a12*a21)/a11
	return (0.5*chibare/(1.0-U*chibare),0.5*chi_nbare/(1.0-U*chi_nbare))

def fermi_derivative(T,e):
	if abs(e/T)>10.0:return 0.0
	return -1.0/(4.*T* ( cosh(0.5*e/T) )**2)	

if __name__=="__main__":
	t=1.0;tp=-0.2;U=3;n=1.1
	Nk=100;NT=50
	kls=linspace(-pi,pi,Nk)
	Tls=linspace(1e-3,0.5,NT)
	mbz=magneticBz(kls)
	(ksiEkp,ksiEkm)=ksiekpm_eval(mbz)
	ksiEkpcp=ksiEkp
	Dlt_mu=1e-4;Dlt_m=1e-4;Nit=300;

	fcheck=open("scsiterations.dat",'w')
	print>>fcheck,"#iter#,","mu-muprime,","mQ-mQprime"
	fmag=open("mag_neq"+str(n)+".txt",'w')
	print>>fmag,"#occupation n=",n,"U=",U
	print>>fmag,"#T,|mQ|,mu,iteration times"
	fchi=open("chi"+str(n)+".txt",'w')
	print>>fchi,"#T,chi_para,chi_n"

	for T in Tls:
	    ksiEkp=ksiEkpcp
	    ksiEkp=ksiEkp+U*n/2.0
	    print "Temperature:",T
	    mu=0.;mQ=0.5
	    it=0
	    while True:
		(mup,mQp)=scseqs(mu,mQ,n,T,U)	
		it+=1
		if (abs(mup-mu)<Dlt_mu and abs(mQp-mQ)<Dlt_m) or (it==Nit):
		    print>>fmag,T,mQp,mup,it
		    print>>fcheck,it,mu-mup,mQ-mQp
		    break
		mu=(mu+mup)/2.0
		mQ=(mQ+mQp)/2.0
	    chi_para,chi_n=spin_susc(mu,mQ,n,T)
	    print>>fchi,T,chi_para,chi_n
