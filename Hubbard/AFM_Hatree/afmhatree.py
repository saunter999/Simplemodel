#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize

def Banddisp(kx,ky,Dta):
	return sqrt( (2.*t*(cos(kx)+cos(ky)))**2+Dta**2)

def Deltascfeq(Delta,U):
	kp=linspace(-pi,pi,Nk)
	scfsumk=0.0
	for kx in kp:
	    for ky in kp:
		scfsumk+=1./Banddisp(kx,ky,Delta)
	scfsumk=scfsumk/Nk**2
	return scfsumk-1./U

if __name__=="__main__":
     ## Ploting SDW band structure along high symmetry lines
     ## in the Magnetic Brillion zone.
     t=1.0
     Dta=0.5
     Nk=200
     plotflag=0
     Gamma=array([0.0,0.0]);M=array([pi,0.0]);Mpm=array([0.0,pi])
     symklist=array([Gamma,M,Mpm])
     Dir=array([M-Gamma,Mpm-M,Gamma-Mpm])
     Ekup=[]
     Ekdn=[]
     for idk,symk in enumerate(symklist):
	 for i in arange(Nk):
	     k=symk+Dir[idk]*float(i)/Nk
	     Ekup.append(Banddisp(k[0],k[1],Dta))
	     Ekdn.append(-Banddisp(k[0],k[1],Dta))
     Ekdn=array(Ekdn)
     Ekup=array(Ekup)
     knum=arange(Ekup.shape[0])
     if(plotflag==1):
       figure(1)
       title(r"$\Delta_{sdw}=$"+str(Dta))
       ypos=-4.9*t
       axvline(x=0*Ekup.shape[0]/3.0,c='r')
       text(0*Ekup.shape[0]/3.0-1, ypos, r"$\Gamma$", fontsize=20)
       axvline(x=1*Ekup.shape[0]/3.0,c='r')
       text(1*Ekup.shape[0]/3.0-1, ypos, r"$M$", fontsize=20)
       axvline(x=2*Ekup.shape[0]/3.0,c='r')
       text(2*Ekup.shape[0]/3.0-1, ypos, r"$M^{\prime}$", fontsize=20)
       axvline(x=3*Ekup.shape[0]/3.0,c='r')
       text(3*Ekup.shape[0]/3.0-1, ypos, r"$\Gamma$", fontsize=20)
       plot(knum,Ekup,label="up")
       plot(knum,Ekdn,label="dn")
       xlabel(r"$k$",size=20)
       ylabel(r"$E(k)$",size=20)
       grid()
       legend(loc=0)
       savefig("Bandstru_SDW.png")
     
     ## Determine Delta_sdw self-consistently by minimizing the ground state energy.##
     Nk=1200
     Ulist=linspace(0.15,0.9,30)
     Delta=[]
     fUdel=open("UDeltaSDW.dat","w")
     print>>fUdel,"#U","Delta_sdw"
     for idU,U in enumerate(Ulist):
	 print str(idU)+" iteration of Ulist"
	 Dlt=optimize.brentq(Deltascfeq,1.0e-15,5.0,args=(U,))
         Delta.append(Dlt)
	 print>>fUdel,U,Dlt  
     Delta=array(Delta)
     if(plotflag==1):
       figure(2)
       plot(Ulist,Delta,"r*")
       grid()
       xlabel(r"$U$",size=20)
       ylabel(r"$\Delta_{sdw}$",size=20)
       savefig("U_DeltaSDW.png")
       show()


