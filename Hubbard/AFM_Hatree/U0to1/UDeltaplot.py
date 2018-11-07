#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import optimize

def func(x,a,b,c):
	return a+b*x**2+c*x**4

data=loadtxt("UDeltaSDW.dat").transpose()
xdata=data[0]
ydata=data[1]
fitparm=optimize.curve_fit(func,xdata,ydata)[0]
print fitparm
fity=[]
extpx=linspace(0.0,max(xdata),20)
for x in extpx:
	fity.append(func(x,fitparm[0],fitparm[1],fitparm[2]))
fity=array(fity)	
plot(xdata,ydata,"r*",label='raw data')
plot(extpx,fity,"go",label="fitting data")
legend(loc=0)
grid()
show()

