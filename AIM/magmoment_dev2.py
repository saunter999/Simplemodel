#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D

def mueq(mu,nd,Edup,Eddn):
	"""
	scf equation determining mu
	"""
	return (1./pi*(arctan(mu-Edup)+arctan(mu-Eddn))-(nd-1.))

def md_eva(nd,U,md): 
	###For given(trial) md we could determine the chemical potential mu by scf equation of nd.###
	Edup=U*(nd-md)/2.
	Eddn=U*(nd+md)/2.
	mu=optimize.brentq(mueq,-10.*Umax,10.*Umax,args=(nd,Edup,Eddn))
	### Using the mu evaluated above,we calculate md 
	mdeva=1./pi*(arctan(mu-Edup)-arctan(mu-Eddn))
	return (mdeva,mu)

def md_scfeq(md_cmp,mdmesh):
	for idx in arange(Nm-1):
	    md_ldif=md_cmp[idx]-mdmesh[idx]
	    md_rdif=md_cmp[idx+1]-mdmesh[idx+1]
	    ###We use the fact that around the scf md point,the difference (md_cmp-md_trial) changes its sign,
	    ###which is used to locate scf md and find the correspounding U and nd for this magnetic phase. 
	    if (md_ldif*md_rdif<0.0):
		return md_cmp[idx] 
	    if (idx==Nm-2) and (md_ldif*md_rdif>0.0):
		return "Noscfmd"
def Ucrext(nd):	
	return pi* (1+(tan(pi/2.0*(nd-1.0)))**2 ) 

def ddos(w,U,nd,md):
	Edup=U*(nd-md)/2.
	Eddn=U*(nd+md)/2.
	Awup=1./(pi*(1.+(w-Edup)**2))
	Awdn=1./(pi*(1.+(w-Eddn)**2))
	return (Awup,Awdn)

if __name__=='__main__':
	"""
	Implementation of study on the formation of magnetic moment as a function of
	nd(occupation of d electron) and U.
	(we set impurity energy level ed=0.0 and hybridization Delta=1) within Hartree
	approximation. 
	The output is the phase diagram of the model as a function of U and nd
	with two phases:magnetic one(square points) and nonmagnetic one
	Q.H. 2015/10/30
	"""
	Nd=50
	Nu=Nd
	Nm=50
	devhf=0.8
	###Generating mesh for nd,U and md 
	ndmesh=linspace(1.-devhf,1.+devhf,Nd)
	Umax=Ucrext(1.+devhf)      ##Diverging as devhf is 1.0. 
	Umesh=linspace(0.,Umax,Nu)
	mdmesh=linspace(0.0,1.5,Nm)
	ndmag=[]
	Umag=[]
	mdmag=zeros((Nd,Nu))
	fmag=open("magn_moment.txt",'w')
	print >>fmag,"#","U","nd","md"
	for id1,nd in enumerate(ndmesh):
	    for id2,U in enumerate(Umesh):
		md_cmp=[]    ###This stores the computed magnetic moment md using the trial md in the mdmesh.
		for md in mdmesh:
			(mdeva,muscf)=md_eva(nd,U,md)
			md_cmp.append(mdeva)	
		md_cmp=array(md_cmp)
		mdscf=md_scfeq(md_cmp,mdmesh)
		if(mdscf!='Noscfmd'):
		    ndmag.append(nd)
		    Umag.append(U)
		    mdmag[id2,id1]=mdscf
		    print>>fmag,U,nd,mdscf
		else:
		    mdmag[id2,id1]=0.0
	fmag.close()
	ndmag=array(ndmag)
	Umag=array(Umag)
	Ucr=[]
	for nd in ndmesh:
		Ucr.append( Ucrext(nd))
	Ucr=array(Ucr)
	figure(1)	
	plot(Umag,ndmag,"s",label="numerics")
	plot(Ucr,ndmesh,label="U_c(nd)(Exact)")
	xlabel("U",size=20)
	ylabel("nd",size=20)
	legend(loc=0)
	grid()
	savefig("magphasedgm.png")
	
	figure(2)
	ndmesh,Umesh=meshgrid(ndmesh,Umesh)
	ax = plt.axes(projection='3d')
	surf=ax.plot_surface(ndmesh,Umesh,mdmag,rstride=1,cstride=1,linewidth=0,cmap=cm.coolwarm, antialiased=False)
	ax.set_xlabel(r'$n_d$',size=20)
#	ax.set_xlim3d(0, 10)
	ax.set_ylabel('U',size=20)
#	ax.set_ylim3d(0, Umax)
	ax.set_zlabel(r'$m_d$',size=20)
#	ax.set_zlim3d(0, 1)
	colorbar(surf, shrink=0.5, aspect=10)
	savefig("magplot.png")

	figure(3)
	md_cmp=[];nd=1.3;U=8.;
	for md in mdmesh:
	    (mdeva,mu)=md_eva(nd,U,md)    
	    md_cmp.append(mdeva)	
	md_cmp=array(md_cmp)
	mdscf=md_scfeq(md_cmp,mdmesh)
	(mdeva,muscf)=md_eva(nd,U,mdscf)
	print mdscf,muscf
	wmesh=linspace(-20,20,200)
	Delta=1e-1
	Awup=[];Awdn=[];Awatom=[];	
	for w in wmesh:
	    (Aw1,Aw2)=ddos(w,U,nd,mdscf)
	    Awup.append(Aw1)
	    Awdn.append(Aw2)
	    Awatom.append(Delta/(pi*(w**2+Delta**2)))
	Awup=array(Awup)
	Awdn=array(Awdn)
	Awatom=array(Awatom)
	plot(wmesh,Awup,label="Awup")
	plot(wmesh,Awdn,label="Awdn")
	plot(wmesh,Awatom,label="Awatom"+r'$U=\Delta=0$')
	axvline(lw=2,x=muscf,c='k',ls='--',label=r"$\mu$")
	legend(loc=0)
	xlabel(r"$\omega$",size=20)
	ylabel("DOS",size=20)
	grid()
	savefig("Aw_dorb.png")
	show()
				
			
 
