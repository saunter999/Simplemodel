#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize

def mueq(mu,nd,Edup,Eddn):
	"""
	scf equation determining mu
	"""
	return (1./pi*(arctan(mu-Edup)+arctan(mu-Eddn))-(nd-1.))

if __name__=='__main__':
	"""
	Implementation of study on the formation of magnetic moment as a function of
	nd(occupation of d electron) and U.
	(we set impurity energy level ed=0.0 and hybridization Delta=1) within Hartree
	approximation. 
	The output is the phase diagram of the model as a function of U and nd
	with two phases:magnetic one(square points) and nonmagnetic one
	Q.H. 2015/10/28
	"""
	Nu=100
	Nm=50
	devhf=0.8
	###Generating mesh for nd,U and md 
	ndmesh=linspace(1.-devhf,1.+devhf,50)
	Umax=pi* (1+(tan(pi/2.0*devhf))**2 ) ##Diverging as devhf is 1.0. 
	Umesh=linspace(0.,Umax,Nu)
	mdmesh=linspace(0.0,1.5,Nm)
	ndmag=[]
	Umag=[]
	fmag=open("magn_moment.txt",'w')
	print >>fmag,"#","U","nd","md"
	for nd in ndmesh:
	    for U in Umesh:
		md_cmp=[]    ###This stores the computed magnetic moment md using the trial md in the mdmesh.
		for md in mdmesh:
			###For given(trial) md we could determine the chemical potential mu by scf equation of nd.###
			Edup=U*(nd-md)/2.
			Eddn=U*(nd+md)/2.
			mu=optimize.brentq(mueq,-10.*Umax,10.*Umax,args=(nd,Edup,Eddn))
			### Using the mu evaluated above,we calculate md 
			mdeva=1./pi*(arctan(mu-Edup)-arctan(mu-Eddn))
			md_cmp.append(mdeva)	
		md_cmp=array(md_cmp)
	
		for idx in arange(Nm-1):
		    md_ldif=md_cmp[idx]-mdmesh[idx]
		    md_rdif=md_cmp[idx+1]-mdmesh[idx+1]
		    ###We use the fact that around the scf md point,the difference (md_cmp-md_trial) changes its sign,
		    ###which is used to locate scf md and find the correspounding U and nd for this magnetic phase. 
		    if (md_ldif*md_rdif<0.0):
			    ndmag.append(nd)
			    Umag.append(U)
			    print>>fmag,U,nd,mdeva
			    break
	fmag.close()
	ndmag=array(ndmag)
	Umag=array(Umag)
	Ucr=[]
	for nd in ndmesh:
		Ucr.append(pi* (1+(tan(pi/2.0*(nd-1.0)))**2 ) )
	Ucr=array(Ucr)
	
	plot(Umag,ndmag,"s",label="numerics")
	plot(Ucr,ndmesh,label="U_c(nd)(Exact)")
	xlabel("U",size=20)
	ylabel("nd",size=20)
	legend(loc=0)
	grid()
	savefig("magmoment_phasediagram.png")
	show()
				
			
 
