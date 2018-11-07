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
	   
def fermi_zeroT(E):
	"""
	 fermi function:when |x| in exp(x) is very large,it
	 cause overflow problem in python.So we use the following
	 method.
	"""
	if (E<0.0):
	    return 1.0
	if E==0.0:
	    return 0.5
	if (E>0.0):
	    return 0.0


def fmmufunc(mu,m,n,U):
	nup=0.;ndn=0.;
	for ek in Ekls:
	    Eupk=ek-mu+U*n/2.-U*m
	    Ednk=ek-mu+U*n/2.+U*m
	    nup+=fermi_zeroT(Eupk)
	    ndn+=fermi_zeroT(Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	return nup+ndn-n
	
def FM_freeeng(mu,m,n,U): 
	F1=-U*(n**2/4.0-m**2)
	F2=0.0
	for ek in Ekls:
	    Eupk=ek-mu+U*n/2.-U*m
	    Ednk=ek-mu+U*n/2.+U*m
	    logup=loge_zeroT(Eupk)
	    F2+=logup
	    logdn=loge_zeroT(Ednk)
	    F2+=logdn
	F2=-1.0*F2/Nk**2;
	F=F1+F2
	return F

def afmmufunc(mu,mQ,n,U):
	nup=0.;ndn=0.;
	for idx,e in enumerate(ksiEkp):
	    Dk=Dk_eval(ksiEkm[idx],mQ) 
	    Eupk=e+Dk-mu 
	    Ednk=e-Dk-mu 
	    nup+=fermi_zeroT(Eupk)
	    ndn+=fermi_zeroT(Ednk)
	nup=2.0*nup/Nk**2;
	ndn=2.0*ndn/Nk**2;
	return nup+ndn-n

def AFM_freeeng(mu,mQ,n,U): 
	F1=-U*(n**2/4.0-mQ**2)
	F2=0.0
	for idx,e in enumerate(ksiEkp):
	    Dk=Dk_eval(ksiEkm[idx],mQ) 
	    Eupk=e+Dk-mu 
	    Ednk=e-Dk-mu 
	    logup=loge_zeroT(Eupk)
	    F2+=logup
	    logdn=loge_zeroT(Ednk)
	    F2+=logdn
	F2=-2.0*F2/Nk**2;
	F=F1+F2
	return F

def loge_zeroT(e):
	if e<0:return -e
	if e>=0:return 0.0

if __name__=="__main__":
	t=1.0;tp=-0.0;n=0.9
	Nk=59; ###it turns out it is more accurate to determine chemical potentail at half filling for odd Nk because a bunch of k points are on the fermi surface for even k
	Uls=[1.,3.,6.,9.,12.]
	kls=linspace(-pi,pi-0.001,Nk)
	Ekls=dispersion(kls)
	mbz=magneticBz(kls)
	(ksiEkp,ksiEkm)=ksiekpm_eval(mbz)
	ksiEkpcp=ksiEkp
	Dlt_mu=1e-4;Dlt_m=1e-4;


	mls=linspace(0,0.5*n-0.001,20)
	for U in Uls:
	    print "U",U
	    ksiEkp=ksiEkpcp
	    ksiEkp=ksiEkp+U*n/2.0
	    fmeng=open("fmeng_neq"+str(n)+"Ueq"+str(U)+".txt",'w')
	    afmeng=open("afmeng_neq"+str(n)+"Ueq"+str(U)+".txt",'w')
	    print>>fmeng,"#occupation n=",n,"U=",U
	    print>>fmeng,"#|m|,F_fm"
	    print>>afmeng,"#occupation n=",n,"U=",U
	    print>>afmeng,"#|mQ|,F_afm"
	    for m in mls:
	        mu=optimize.brentq(fmmufunc,-30.,30.,args=(m,n,U))
		print "FM---","m,mu=",m,mu
		F_fm=FM_freeeng(mu,m,n,U)
		print>>fmeng,m,F_fm
	        mu=optimize.brentq(afmmufunc,-30.,30.,args=(m,n,U))
		print "AFM---","m,mu=",m,mu
		F_afm=AFM_freeeng(mu,m,n,U)
		print>>afmeng,m,F_afm
