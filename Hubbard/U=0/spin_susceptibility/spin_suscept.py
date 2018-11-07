#!/usr/bin/env python
from scipy import *
from pylab import *
def Banddisp(kx,ky):
       return -2.*t*( cos(kx)+cos(ky) )

def fermi(E,T):
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

def Dos(emesh):
	dlt=0.01
	Ek=[]
	for kx in kp:
	    for ky in kp:
		Ek.append(-2.*t*(cos(kx)+cos(ky))) 
	Ek=array(Ek)
	dos=[]
	for omg in emesh:
	  Gloc=0.0j
	  for ek in Ek:
	      Gloc+=1./(omg-ek+1.0j*dlt)
	  Gloc=Gloc/Nk**2
	  dos.append( -(Gloc.imag)/pi)
	return dos

if __name__=="__main__":
     """
	We evaluate the static spin susceptibitilty along (0,0)->(pi,pi) and (0,0)->(pi,0)  
        for hubbard model with U=0 at half filling
     """
     t=1.0
     Nk=300
     kp=linspace(-pi,pi,Nk)
     emesh=linspace(-5,5,100)
     #dos=Dos(emesh)
     #figure(1)
     #plot(emesh,dos)
     #savefig("dos.png")

     ## q path 1 :(0,0)->(pi,pi)  
     figure(2)
     Nq=10
     onedq=[1e-5]
     qls=[[1e-5,1e-5]]
     for i in range(Nq):
	 qp=float(i+0.99)/Nq*pi
	 onedq.append(qp)
	 qls.append([qp,qp])
	
     qls=array(qls) 
     onedq=array(onedq)
     mu=0.0;
     Tmesh=linspace(1e-3,0.1,3)
     for T in Tmesh:
       chispin=[]	
       for q in qls:
	   print "Evaluating static spin susceptibility for q=",q
	   sumchi=0.0	
	   for kx in kp:
	       for ky in kp:
		   ek=Banddisp(kx,ky)
		   ekq=Banddisp(kx+q[0],ky+q[1])
		   sumchi+=- ( fermi(ekq-mu,T)-fermi(ek-mu,T) )/ (ekq-ek)
	   sumchi=sumchi/Nk**2
	   chispin.append(sumchi)
       plot(onedq,chispin,'o-',label='T='+str(T))
     xlabel("q((0,0)->(pi,pi))",size=19)
     ylabel(r"$\chi_0$"+"(q,"+r"$\Omega=0$"+")",size=19)
     legend(loc=0.0)
     savefig("spinsusc1.png")

     ## q path 2 :(0,0)->(pi,0)  
     figure(3)
     Nq=10
     onedq=[1e-5]
     qls=[[1e-5,0]]
     for i in range(Nq):
	 qp=float(i+0.99)/Nq*pi
	 onedq.append(qp)
	 qls.append([qp,0])
	
     qls=array(qls) 
     onedq=array(onedq)
     mu=0.0;
     Tmesh=linspace(1e-3,0.1,3)
     for T in Tmesh:
       chispin=[]	
       for q in qls:
	   print "Evaluating static spin susceptibility for q=",q
	   sumchi=0.0	
	   for kx in kp:
	       for ky in kp:
		   ek=Banddisp(kx,ky)
		   ekq=Banddisp(kx+q[0],ky+q[1])
		   sumchi+=- ( fermi(ekq-mu,T)-fermi(ek-mu,T) )/ (ekq-ek)
	   sumchi=sumchi/Nk**2
	   chispin.append(sumchi)
       plot(onedq,chispin,'o-',label='T='+str(T))
     xlabel("q((0,0)->(pi,0))",size=19)
     ylabel(r"$\chi_0$"+"(q,"+r"$\Omega=0$"+")",size=19)
     legend(loc=0.0)
     savefig("spinsusc2.png")
	
     ##Contour plot
    
     figure(4)
     qp=linspace(-pi,pi,Nq)
     T=0.01
     chispin=zeros((Nq,Nq))
     for idx,qx in enumerate(qp):
	 for idy,qy in enumerate(qp):
	   print "iqx,iqy=",idx,idy
	   sumchi=0.0	
	   for kx in kp:
	       for ky in kp:
		   ek=Banddisp(kx,ky)
		   ekq=Banddisp(kx+qx,ky+qy)
		   sumchi+=- ( fermi(ekq-mu,T)-fermi(ek-mu,T) )/ (ekq-ek)
	   sumchi=sumchi/Nk**2
	   chispin[idx,idy]=sumchi
     imshow(chispin,origin='lower',extent=[-pi, pi, -pi, pi],interpolation='bilinear',aspect='equal')
     colorbar()
     xlabel("qx",size=19)
     ylabel("qy",size=19)
     title(r"$\chi_0$"+"(q,"+r"$\Omega=0$"+")"+"T="+str(T),size=19)
     savefig("contour_spinsusc.png")

     show()


