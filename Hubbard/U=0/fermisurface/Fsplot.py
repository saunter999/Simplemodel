#!/usr/bin/env python
from scipy import *
from pylab import *
def Banddisp():
       bands=zeros((Nk,Nk))
       for idx,kx in enumerate(kp):
	 for idy,ky in enumerate(kp):
	     bands[idx,idy]=-2.*t*( cos(kx)+cos(ky) )
       return bands

if __name__=="__main__":
     """
	plotting the Fermi surface for several mu
        for hubbard model with U=0
     """
     t=1.0
     Nk=800
     kp=linspace(-pi,pi,Nk)
     muls=[-1.0,0.0,1.0]
     eta=5e-2
     bands=Banddisp()
     for idx,mu in enumerate(muls):
       subplot(1,3,idx+1)
       Fsw=zeros((Nk,Nk))
       for ix in range(Nk):
         for iy in range(Nk):
	   if ( abs(bands[ix,iy]-mu) <eta):
		Fsw[ix,iy]=10.0
	   else:
	        Fsw[ix,iy]=0.0
#	   Fsw[idx,idy]=exp(-alpha*(Banddisp(kx,ky)-mu)**2)
       imshow(Fsw,origin='lower',extent=[-pi, pi, -pi, pi],interpolation='bilinear',aspect='equal')
       title(r"$\mu$="+str(mu),size=25)
       xlabel(r"$k_x$",size=20)
       ylabel(r"$k_y$",size=20)
     savefig("Fermisurface_mu.png")
     show()


