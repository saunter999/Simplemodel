#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize

def fermi(E,T):
    if (E/T<-10.0):
	return 1.0;
    if (E/T>10.0):
	return 0.0;
    return 1.0/(exp(E/T)+1.0)

def ocup(mu,T,h):
    """
    return the occupation of electrons with spin up and down
    nelec[0]->the ocp of electrons with spin up
    nelec[1]->the ocp of electrons with spin dn
    mu:chemical potential
    T:temperature
    h:external magnetic field
    """
    nelec=zeros(2)
    for E in Ek:
	Eup=E+0.5*g*h-mu
	Edn=E-0.5*g*h-mu
	nelec[0]+=fermi(Eup,T)
	nelec[1]+=fermi(Edn,T)
    nelec=nelec/Nk**2
    return nelec

def mu_scfeq(mu,T,h):
    """
    Defining the function as the equation of mu by requiring 
    nup+ndn=n(filling factor)
    """
    nelec=zeros(2)
    nelec=ocup(mu,T,h)
    return (nelec[0]+nelec[1]-n)

def Dos(omg):
	Gloc=0.0j
	for ek in Ek:
	    Gloc+=1./(omg-ek+1.0j*dlt)
	Gloc=Gloc/Nk**2
	return -(Gloc.imag)/pi

if __name__=='__main__':
    """
	This script is an implementation of calculating the magnetic response of one-band 
        tight-binding model in 2d(or U=0 hubbard model) to an external uniform magnetic
	field h.We implement the caluation of Mz and Chiz as a function (T,h).
	As an benchmark,we compare the numerical result to the analytic result in the limit
	of zero T and h,which are in good match.
	Q.H. 2015.10.18
    """
    Nk=200;NT=20;Nh=20;Tmax=0.5;hmax=0.3 
    Ts=1e-2;hs=1e-8
    # To evaluate chemical potentail mu accurately,Nk needs to be greater than 50. 
    t=1.0;dlt=0.01;n=1.1        #n->> filling factor
    g=2.0              #(g->>lande factor,$\mu_b$==1 here.)
    kmesh=linspace(-pi,pi,Nk)
    ##   numerical test for Ts and hs ##
    ##   We tend to compare the analytical result of Mz and Chiz for small T and h with 
    ##   the numerial implementation.The test shows that the value Ts is tricky.For very 
    ##   small Ts,the relative error is big,for Ts around 1e-2,the analytic and numerial
    ##   results are in good match.For hs,even if it is very small,they match well.
    Tmesh=linspace(Ts,Tmax,NT)
    hmesh=linspace(hs,hmax,Nh) 
    dh=hmesh[1]-hmesh[0]
    Mz=zeros((NT,Nh))  #Magnetization 
    chiz=zeros((NT,Nh)) #Susceptibility
    nelec=zeros(2)
    Ek=[]
    for kx in kmesh:
        for ky in kmesh:
	    Ek.append(-2.*t*(cos(kx)+cos(ky)))  #dispersion of the unperturbed system.
    Ek=array(Ek)
    #Determination of the chemical potential of the unperturbed system.
    mu0=optimize.brentq(mu_scfeq,-8.0*t,8.0*t,args=(1e-8,0.0)) 
    for iT,T in enumerate(Tmesh):
	for ih,h in enumerate(hmesh):
	  mu=optimize.brentq(mu_scfeq,-8.0*t,8.0*t,args=(T,h))
#	  print "for filling n=",n,"with T=",T,",h=",h,"...mu=",mu
	  nelec=ocup(mu,T,h)
	  Mz[iT,ih]=-0.5*g*(nelec[0]-nelec[1])
    for i in arange(NT):
    	for j in arange(Nh):
	    if j==0:
		#when h is very small,we evaulte the susc against (Mz=0,h=0) point.
	    	chiz[i,j]=(Mz[i,j]-0.0)/hmesh[0]
	    else:
	    	chiz[i,j]=(Mz[i,j]-Mz[i,j-1])/dh
     #comparing to the analytical result for T=0 and h=0
    Mzana=0.5*g**2*Dos(mu0)*hmesh[0]
    print "Comparing Mz(T,h very small) evaluated by numericals and exact result"
    print "Mz_Ana=",Mzana,"Mz_num=",Mz[0,0],"rela_err=",abs(Mzana-Mz[0,0])/Mz[0,0]
    chizana=Mzana/hmesh[0]
    print "Comparing chiz(T,h very small) evaluated by numericals and exact result"
    print "chiz_Ana=",chizana,"chiz_num=",chiz[0,0],"rela_err=",abs(chizana-chiz[0,0])/chiz[0,0]

    fmag=open("magn.txt","w")
    print >>fmag,"# Temp,hmesh,Mz"
    fsus=open("suscept.txt","w")
    print >>fsus,"# Temp,hmesh,chiz"
    for i in arange(NT):
	for j in arange(Nh):
	    print>>fmag,Tmesh[i],hmesh[j],Mz[i,j]
	    print>>fsus,Tmesh[i],hmesh[j],chiz[i,j]


