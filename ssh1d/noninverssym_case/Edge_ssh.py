#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA




if __name__=='__main__':
	"""
	Edge states in 1d ssh model including the staggered mass term in open boundary condition:
	1)m=0,dlt<0,it is topological nontrival,there is a zero mode,which is localized at the two boundaries
	(by inspecting the wavefucntion distribution);
	2)m=0,dlt>0,it is topological trival,no zero mode.
	3)m/=0: no zero mode.  
	"""
	N=50            #number of unit cell.
	t=1.0 
	Npara=50          #Choose Npara to be even to make sure m=0&dlt=0 is included.  
	mls=linspace(-t,t,Npara+1)       
	dltls=linspace(-t,t,Npara+1)       
	ham=zeros((2*N,2*N))
	x=zeros(2*N)
	f=open("zeromode.dat",'w')
	print>>f,"#m","dlt","zero mode--:1->yes;-1->no"
	for m in mls:
	  for dlt in dltls:
	    for i in range(2*N):
		    if(i==0):
		       ham[i,i+1]=t+dlt
		       ham[i,i]=m	
		    if(i==2*N-1):
		       ham[i,i-1]=t+dlt
		       ham[i,i]=-m	
		    if(i>0) and (i<2*N-1):	
		      if i%2==1:
			 ham[i,i]=-m	
			 ham[i,i-1]=t+dlt
			 ham[i,i+1]=t-dlt
		      if i%2==0:
			 ham[i,i]=m	
			 ham[i,i-1]=t-dlt
			 ham[i,i+1]=t+dlt
	    Eig,v=LA.eig(ham)
	    idx=Eig.argsort()
	    Eig=Eig[idx]
	    v=v[:,idx]
	    for idx,E in enumerate(Eig):
		if abs(E)<abs(0.0001):
		  print >>f,m,dlt,1
		  break
		else:
		  if(idx==2*N-1):print >>f,m,dlt,-1
	f.close()
	path='./zeromode.dat'
	dat=loadtxt(path).transpose()
	zerom=dat[2]
	dim=int(sqrt(len(zerom)))
	mapzerom=zeros((dim,dim))
	for j in range(dim):
	    for i in range(dim):
		mapzerom[i,j]=zerom[i+j*dim]
	imshow(mapzerom,extent=[-t,t,-t,t],cmap=cm.jet,origin='lower')
	colorbar(shrink=1.0)
	xlabel(r'$\frac{m}{t}$',size=15)
	ylabel(r'$\frac{\delta t}{t}$',size=15,rotation=0)
	title('map of ' r'$zeromode:-1:NO;1:Yes$')
	savefig("zeromode.png")
	show()
