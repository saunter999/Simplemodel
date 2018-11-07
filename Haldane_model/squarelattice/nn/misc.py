#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
#import linalg as LA
#from linalg import *

print exp(-1j*pi)
print -2.0*exp(-1j*pi)
print conjugate(1j)
(w,v)=LA.eig(diag((1,2,3)))
print w,type(w),w[0]
print v
B=range(5)
for i in reversed(B):
	print i
