#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
import numpy as np

def dx(kx,ky):
	return sin(ky)

def dy(kx,ky):
	return -sin(kx)

def dz(kx,ky,es):
	return c*(2.0-cos(kx)-cos(ky)-es)

def Hin(ky):
	H=zeros((2,2),dtype=complex)
	H=sin(ky)*sigmax+c*(2.0-cos(ky)-es)*sigmaz
	return H

def HacR():
	H=zeros((2,2),dtype=complex)
	H=0.5j*sigmay-c/2.*sigmaz
	return H
	
def HacL():
	H=zeros((2,2),dtype=complex)
	H=-0.5j*sigmay-c/2.*sigmaz
	return H

def set_ham(ky,Nj):
	ham=zeros((2*Nj,2*Nj),dtype=complex)
	for j in range(Nj):
	    HinM=Hin(ky)	
	    ham[2*j,2*j]=HinM[0,0]
	    ham[2*j,2*j+1]=HinM[0,1]
	    ham[2*j+1,2*j]=HinM[1,0]
	    ham[2*j+1,2*j+1]=HinM[1,1]
	    if (j !=0) and (j !=(Nj-1)):
		HacLM=HacL()	
		HacRM=HacR()	
		ham[2*j,2*j-2]=HacLM[0,0]
		ham[2*j,2*j-1]=HacLM[0,1]
		ham[2*j+1,2*j-2]=HacLM[1,0]
		ham[2*j+1,2*j-1]=HacLM[1,1]
		ham[2*j,2*j+2]=HacRM[0,0]
		ham[2*j,2*j+3]=HacRM[0,1]
		ham[2*j+1,2*j+2]=HacRM[1,0]
		ham[2*j+1,2*j+3]=HacRM[1,1]
	    if j==0:
		HacRM=HacR()	
		ham[0,2]=HacRM[0,0]
		ham[0,3]=HacRM[0,1]
		ham[1,2]=HacRM[1,0]
		ham[1,3]=HacRM[1,1]
	    if j==Nj-1:
		HacLM=HacL()	
		ham[2*j,2*j-2]=HacLM[0,0]	
		ham[2*j,2*j-1]=HacLM[0,1]	
		ham[2*j+1,2*j-2]=HacLM[1,0]	
		ham[2*j+1,2*j-1]=HacLM[1,1]	
	return ham
      
def Ham_solver(H):
	(Eig,v)=LA.eig(H)
	idx=Eig.argsort() #This one is needed as the eigenvaules returned by LA.eig aren't in the order of -/+ levels.
	Eig=Eig[idx]
	v=v[:,idx]
	return (Eig,v)

if __name__=='__main__':
	"""
	  Implementation of intepreation of hall conductivity in terms of winding number...
	  The model considered here is taken from eq(7) in the paper PRB 74,085308(2006) XL.Qi
	"""
	sigmax=array([[0.,1.],[1.,0.]]);sigmay=array([[0.,-1.0j],[1.0j,0.]]);sigmaz=array([[1.0,0.],[0.,-1.]])
	### Model parameter ###
	c=1.0
        ##########################
        Ne=66;elist=linspace(-6.,6.,Ne)
        f=open("taskcont_TImodel1.dat","r")
        lines=f.readlines()
        l=lines[0]
        taskcont=l.split()
        itask1=taskcont[0]
        itask2=taskcont[1]
	###1).winding number of the d vector....
	if itask1=="1":
            print "--------------------------------------------------"
	    print "task1:calculating the winding number of the d vector:running" 
	    Nk=600;kp=linspace(-pi,pi,Nk)
	    f=open("TI_windnum.dat","w")
	    JM=zeros((3,3)) ###Jacobi matrix
	    WNe=[]
	    for es in elist:
	      WN=0.0
	      for ky in kp:
		for kx in kp:
		    d1=dx(kx,ky);d2=dy(kx,ky);d3=dz(kx,ky,es)
		    JM[0,0]=d1;JM[0,1]=d2;JM[0,2]=d3;
		    JM[1,0]=0.;JM[1,1]=-cos(kx);JM[1,2]=c*sin(kx)
		    JM[2,0]=cos(ky);JM[2,1]=0.;JM[2,2]=c*sin(ky)
		    WN+=1./(sqrt(d1**2+d2**2+d3**2))**3*LA.det(JM)
	      WN=pi*WN/Nk**2
	      print>>f,es,WN
	      WNe.append(WN)
        else:
            print "--------------------------------------------------"
	    print "task1:calculating the winding number of the d vector:NOT running" 

	###2).Calculating the edge states.
	es=2.0;TIind=1
	if itask2=="1":
            print "--------------------------------------------------"
	    print "task2:calculating the edge state:running"
	    Nk=100;kp=linspace(-pi,pi,Nk)
	    Nj=100;  #number of chains in x direction
	    print "model parameter:"
	    print "c=",c,"es=",es
	    print "Nk=",Nk
	    print "Nj=",Nj
	    Eigky=zeros((2*Nj,Nk))
	    for indk,ky in enumerate(kp):
		ham=zeros((2*Nj,2*Nj),dtype=complex)
		ham=set_ham(ky,Nj)
		print "indk=",indk,",","ky=",ky
		Eig,v=Ham_solver(ham)
		###we need to transpose v to make that v[i,:] is the eigenvector of eigvalue w[i]. 
		v=v.transpose()
		for j in range(2*Nj):
		    Eigky[j,indk]=Eig[j].real
		if TIind==1:      #####The edge state energy is the (Nj-1)'th and Nj'th eigenvaules of the Hamiltonian for a given ky,we need to switch Nj-1 and Nj for half of ky point when it is TI
		  if ky>0:
		    E_edge=Eigky[Nj-1,indk]
		    Eigky[Nj-1,indk]=Eigky[Nj,indk]
		    Eigky[Nj,indk]=E_edge
		  if indk==Nk/2-3:
		    print "Edge ky=",kp[indk]
		    figure(1)
		    den_Nj=zeros(Nj)
		    den_Njminus1=zeros(Nj)
		    for j in arange(Nj):
			den_Nj[j]=abs(v[Nj,2*j])+abs(v[Nj,2*j+1])   ## density distribution of the edge state for a ky,here we sum over the interal degree of freedom,i.e the peseudospin space.
			den_Njminus1[j]=abs(v[Nj-1,2*j])+abs(v[Nj-1,2*j+1])
		    plot(arange(Nj),den_Nj,label="E(Nj),k<0") ##0end
		    plot(arange(Nj),den_Njminus1,label="E(Nj-1),k<0")  ##Nj end
		  if indk==Nk/2+3:
		    v1=np.copy(v[Nj-1])
		    v0=np.copy(v[Nj])
		    v[Nj-1]=v0
		    v[Nj]=v1
		    print "Edge ky=",kp[indk]
		    for j in arange(Nj):
			den_Nj[j]=abs(v[Nj,2*j])+abs(v[Nj,2*j+1])
			den_Njminus1[j]=abs(v[Nj-1,2*j])+abs(v[Nj-1,2*j+1])
		    plot(arange(Nj),den_Nj,label="E(Nj),k>0",ls="--",marker='*',markersize=6) ##Nj end without switching
		    plot(arange(Nj),den_Njminus1,label="E(Nj-1),k>0",ls="--",marker="8",markersize=6)  ##0end without switching
		    legend(loc=0)
	    	    xlabel(r"$x_a$",size=20)
	            ylabel(r"$|\psi(x_a)|^2$",size=20)
		    grid()
	    	    savefig("Edgestate_density.png")
	    figure(2)		   
	    scatter(kp,Eigky[Nj-1,:],c="g")
	    scatter(kp,Eigky[Nj,:],c="b")
	    for j in range(Nj-1):
	        plot(kp,Eigky[j,:],"r")
	    for j in range(Nj+1,2*Nj):
	        plot(kp,Eigky[j,:],"r")
	    grid()
	    title("strip calculation with "+"es="+str(es),size=20)
	    xlabel(r"$k_y$",size=20)
	    ylabel(r"$E$",size=20)
	    savefig("strip_cal.png")
	else:
            print "--------------------------------------------------"
	    print "task2:calculating the edge state:NOT running"
	show()
