#!/usr/bin/env python
from scipy import *
from pylab import *



if __name__=='__main__':
	"""	
	Eigenvaules of Two-site Hubbard model in the sector Sz=0 for half-filled case..
	"""
	t=1.;NU=20
	Uls=linspace(0,20.,NU)
	E0=zeros(NU)
	E1=Uls
	E2=0.5*( Uls-sqrt(16*t**2+Uls**2) )
	E3=0.5*( Uls+sqrt(16*t**2+Uls**2) )
	plot(Uls,E0,label="E0")
	plot(Uls,E1,label="E1")
	plot(Uls,E2,label="E2")
	plot(Uls,E3,label="E3")
	xlabel("U/t",size=19)
	ylabel("E/t",size=19)
	grid()
	legend(loc=0)
	savefig("Hub_twosite.pdf")
	show()
	
