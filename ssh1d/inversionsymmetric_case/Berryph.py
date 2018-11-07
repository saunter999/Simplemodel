#!/usr/bin/env python
from scipy import *
from pylab import *
import numpy as np
from numpy import linalg as LA
from scipy.integrate import quad


def berrycnt(k,a,b):
	"""
	 a=t;b=dlt
	"""
	norm2=4*(a*cos(k)**2+b*sin(k)**2)**2+(b-a)**2*sin(2*k)**2+4*(a**2*cos(k)**2+b**2*sin(k)**2)
	return 2*(b-a)*(a*(1.+cos(2*k))+b*(cos(2*k)-1.))/norm2
if __name__=='__main__':
	"""
	  Evaulation of the Berry phase of ssh model using the Bloch wavefunction of the Bloch Hamiltonian(or the Kernal of Hamiltonian).
	  For dlt<0,the berry phase is pi;dlt>0,the berry phase is zero.
	  so the former is topological nontrival while the latter is trival.
	  This is in accord with the edge state calculation.
	"""
	t=1
	Nk=100
	kl=linspace(0,pi,Nk)	
	dk=pi/Nk
	Hamk=zeros((2,2),dtype=complex)
	uk_ocp=zeros((2,Nk),dtype=complex)  #used to store the eigenvector of the Bloch Hamiltonian at each k point. 
	uk_unocp=zeros((2,Nk),dtype=complex)
	ukext_ocp=zeros((2,Nk),dtype=complex)
	ukext_unocp=zeros((2,Nk),dtype=complex)
	bryph1_ocp=[]
	bryph1_unocp=[]
	bryph2_ocp=[]
	bryph2_unocp=[]
	bryph3_ocp=[]
	dltls=linspace(-t,t,200)
	for dlt in dltls:
	  for idk,k in enumerate(kl):
	     #set the Bloch Hamiltonian.
	     Hamk[0,1]=2.*t*cos(k)**2+2.*dlt*sin(k)**2+1j*(dlt-t)*sin(2.0*k)
	     Hamk[1,0]=conjugate(Hamk[0,1])
	     Eig,v=LA.eig(Hamk)
	     idx=Eig.argsort() 
	     Eig=Eig[idx]
	     v=v[:,idx]
	     uk_ocp[:,idk]=v[:,0]
	     uk_unocp[:,idk]=v[:,1]

	  #method 1):Using the definition of berry phase to calculate it,which turns out to be impractical.
	  brp_ocp=0.0
	  for i in range(1,Nk):
	     udif_k=(uk_ocp[:,i]-uk_ocp[:,i-1])/dk
	     Ak=-1j*np.vdot(uk_ocp[:,i],udif_k)
	     brp_ocp+=Ak*dk
	  #print brp_ocp/pi
	
	  #method 2):Using the analytical formulua for eigenvector 
	  for idk,k in enumerate(kl):
	     ukext_ocp[0,idk]=2*(t*cos(k)**2+dlt*sin(k)**2)+1j*(dlt-t)*sin(2*k) 
	     ukext_ocp[1,idk]=-2*sqrt(t**2*cos(k)**2+dlt**2*sin(k)**2)
	     ukext_unocp[0,idk]=2*(t*cos(k)**2+dlt*sin(k)**2)+1j*(dlt-t)*sin(2*k) 
	     ukext_unocp[1,idk]=2*sqrt(t**2*cos(k)**2+dlt**2*sin(k)**2)
	     norm=sqrt(ukext_ocp[0,idk]*conjugate(ukext_ocp[0,idk])+ukext_ocp[1,idk]**2)
	     ukext_ocp[:,idk]=ukext_ocp[:,idk]/norm     #normalized eigenvector of the lower band
	     ukext_unocp[:,idk]=ukext_unocp[:,idk]/norm #normalized eigenvector of the top band
	  brp_ocp=0.0
	  brp_unocp=0.0
	  for i in range(1,Nk):
	     udif_k=(ukext_ocp[:,i]-ukext_ocp[:,i-1])/dk
	     Ak=-1j*np.vdot(ukext_ocp[:,i],udif_k)
	     brp_ocp+=Ak*dk
	     udif_k=(ukext_unocp[:,i]-ukext_unocp[:,i-1])/dk
	     Ak=-1j*np.vdot(ukext_unocp[:,i],udif_k)
	     brp_unocp+=Ak*dk
	  bryph1_ocp.append(brp_ocp.real) 
	  bryph1_unocp.append(brp_unocp.real) 
	  

	  #method 3):Using the equivalent form introduced by D.Vanderbilt to evaluate it,which is numerically stable.
	  umul1=1.0
	  umul2=1.0
	  for i in range(1,Nk):
	      uovlap1=np.vdot(uk_ocp[:,i-1],uk_ocp[:,i])
	      uovlap2=np.vdot(uk_unocp[:,i-1],uk_unocp[:,i])
	      umul1=umul1*uovlap1	
	      umul2=umul2*uovlap2	
	  bryph2_ocp.append((np.log(umul1)).imag)   
	  bryph2_unocp.append((np.log(umul2)).imag)   
	  
	  #method 4)analytical formula by doing the calculation.
	  bry=quad(berrycnt,0,pi,args=(t,dlt))		
	  bryph3_ocp.append(bry[0])

	bryph1_ocp=array(bryph1_ocp)
	bryph1_unocp=array(bryph1_unocp)
	bryph2_ocp=array(bryph2_ocp)
	bryph2_unocp=array(bryph2_unocp)
	bryph3_ocp=array(bryph3_ocp)
	figure(1)
	bryph1_ocp=(np.around(bryph1_ocp/pi))%2
	bryph1_unocp=(np.around(bryph1_unocp/pi))%2
	plot(dltls,bryph1_ocp,'r-', linewidth=4,label='Bph_dnband') 
	plot(dltls,bryph1_unocp,'g--', linewidth=3,label='Bph_upband') 
	plot(dltls,(bryph1_ocp+bryph1_unocp)%2,'b--', linewidth=2,label='Bph_sum') 
	plot(dltls,bryph3_ocp/pi,'k--', linewidth=1,label='Bph_dnband:exact_result') 
	grid()
	xlabel(r'$\frac{\delta t}{t}$',fontsize=23)
	ylabel(r'$\frac{\phi_{Berry}}{\pi}$',fontsize=23,rotation=0)
	legend(loc='upper right')
        savefig("berry_phase_ssh_method2.png")
	figure(2)
	bryph2_ocp=(np.around(bryph2_ocp/pi))%2
	bryph2_unocp=(np.around(bryph2_unocp/pi))%2
	plot(dltls,bryph2_ocp,'r-',linewidth=3,label='Bph_dnband') 
	plot(dltls,bryph2_unocp,'g--',linewidth=2,label='Bph_upband') 
	plot(dltls,(bryph2_ocp+bryph2_unocp)%2,'b--',linewidth=1,label='Bph_sum') 
	grid()
	xlabel(r'$\frac{\delta t}{t}$',fontsize=23)
	ylabel(r'$\frac{\phi_{Berry}}{\pi}$',fontsize=23,rotation=0)
	legend(loc=0)
        savefig("berry_phase_ssh_method3.png")
	show()
