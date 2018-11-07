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

def afm_sol(mu,mQ,n,U):
	it=0
	while True:
	    (mup,mQp)=afmscseqs(mu,mQ,n,U)	
	    it+=1
	    if (abs(mup-mu)<Dlt_mu and abs(mQp-mQ)<Dlt_m) or (it==Nit):
		print>>afmmag,U,mQp,mup,it
		print>>afmcheck,it,mu-mup,mQ-mQp
		break
	    mu=(mu+mup)/2.0
	    mQ=(mQ+mQp)/2.0
	return (mu,mQ)

def afmscseqs(mu,mQ,n,U):
	mQp=0.0
	for idx,e in enumerate(ksiEkp):
	    Dk=Dk_eval(ksiEkm[idx],mQ) 
	    Eupk=e+Dk-mu 
	    Ednk=e-Dk-mu 
	    mQp+=(fermi_zeroT(Eupk)-fermi_zeroT(Ednk))*(U*mQ)/Dk
	mQp=-1.0*mQp/Nk**2;
	mup=optimize.brentq(afmmufunc,-20.,20.,args=(mQ,n,U)) 
	return array([mup,mQp])

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


	  
def fm_sol(mu,m,n,U):
	it=0
	while True:
	    (mup,mp)=fmscseqs(mu,m,n,U)	
	    it+=1
	    if (abs(mup-mu)<Dlt_mu and abs(mp-m)<Dlt_m) or (it==Nit):
		print>>fmmag,U,mp,mup,it
		print>>fmcheck,it,mu-mup,m-mp
		break
	    mu=(mu+mup)/2.0
	    m=(m+mp)/2.0
	return (mu,m)

def fmscseqs(mu,m,n,U):
	nup=0.0;ndn=0.0
	for ek in Ekls:
	    Eupk=ek-mu+U*n/2.-U*m
	    Ednk=ek-mu+U*n/2.+U*m
	    nup+=fermi_zeroT(Eupk)
	    ndn+=fermi_zeroT(Ednk)
	nup=nup/Nk**2;
	ndn=ndn/Nk**2;
	mp=0.5*(nup-ndn)
	mup=optimize.brentq(fmmufunc,-30.,30.,args=(m,n,U))
	return array([mup,mp])

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
	t=1.0;tp=-0.0;n=1.0
	Nk=59;###it turns out it is more accurate to determine chemical potentail at half filling for odd Nk because a bunch of k points are on the fermi surface for even k
	Uls=[1.,3.,6.,9.,12.]
	kls=linspace(-pi,pi-0.001,Nk)
	Ekls=dispersion(kls)
	mbz=magneticBz(kls)
	(ksiEkp,ksiEkm)=ksiekpm_eval(mbz)
	ksiEkpcp=ksiEkp
	Dlt_mu=1e-4;Dlt_m=1e-4;Nit=1000;

	fmcheck=open("fmscsiterations.txt",'w')
	print>>fmcheck,"#iter#,","mu-muprime,","m-mprime"
	fmmag=open("fmmag_neq"+str(n)+".txt",'w')
	print>>fmmag,"#occupation n=",n
	print>>fmmag,"#U,|m|,mu,iteration times"
	afmcheck=open("afmscsiterations.txt",'w')
	print>>afmcheck,"#iter#,","mu-muprime,","mQ-mQprime"
	afmmag=open("afmmag_neq"+str(n)+".txt",'w')
	print>>afmmag,"#occupation n=",n
	print>>afmmag,"#U,|mQ|,mu,iteration times"
	feng=open("eng_U.txt","w")
	print>>feng,"#U,F_fm,F_afm,F_para"

	for U in Uls:
	    print "U",U
	    mu_para=optimize.brentq(fmmufunc,-30.,30.,args=(0.0,n,U))
	    ksiEkp=ksiEkpcp
	    ksiEkp=ksiEkp+U*n/2.0
	    mu=0.;m=0.5
	    mu,m=fm_sol(mu,m,n,U) 
	    print "mupara,mu_mag,m:",mu_para,mu,m
	    F_fm=FM_freeeng(mu,m,n,U)
	    F_para=FM_freeeng(mu_para,0.0,n,U)
	    mu=0.;mQ=0.5
	    mu,mQ=afm_sol(mu,mQ,n,U)
	    F_afm=AFM_freeeng(mu,mQ,n,U)    
	    print U,F_fm,F_afm,F_para
	    print>>feng,U,F_fm,F_afm,F_para
