#!/usr/bin/env python
from scipy import *
import glob
from pylab import *



if __name__=="__main__":
	f1=glob.glob("mag_neq*txt")
	f2=glob.glob("chi*txt")
	mag=loadtxt(f1[0]).transpose()
	plot(mag[0],mag[1]*30,"*-",label='m')
	chi=loadtxt(f2[0]).transpose()
	plot(chi[0],chi[1],"o-",label=r"$\chi_{para}$")
	plot(chi[0],chi[2],"^-",label=r"$\chi_n$")
	legend(loc=0)
        xlabel("T/ev",size=19)
        ylabel(r"$\chi$"+'or m*30',size=19)
	axvline(x=0.187,c='k',ls='--')
        savefig("mag_chiT.png")
	show()
