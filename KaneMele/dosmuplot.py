#!/usr/bin/env python
from scipy import * 
from pylab import *

data0=loadtxt("doskanemele.dat").transpose()
data1=loadtxt("occ_mu.dat").transpose()
subplot(2,1,1)
title(r"$\lambda=0.13$")
plot(data0[0],data0[3],label=r"$t=1.0$")
xlabel(r"$\omega$",size=20)
ylabel("DOS",size=20)
grid()
legend(loc=0)

subplot(2,1,2)
plot(data1[0],data1[1],'r')
xlabel(r"$\mu$",size=20)
ylabel(r"$n_e$",size=20)
grid()
xlim([-4.0,4.0])
savefig("dosmukanemele.png")

show()
