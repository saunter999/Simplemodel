#!/usr/bin/env python
from scipy import *
import glob
from pylab import *



if __name__=="__main__":
	f=glob.glob("mag_mueq*")
	for fname in f:
	    data=loadtxt(fname).transpose()
	    plot(data[0],data[1],"*-",label=fname)
	    legend(loc=0)
	    if fname==f[-1]: 
	       xlabel("T/ev")
	       ylabel("m")
	       grid()
	       savefig("mag_mu.png")
