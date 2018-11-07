#!/usr/bin/env python
from scipy import *
import numpy as np
import numpy.linalg as linalg

#testing that python stores the eigenvector in columns not in rows.
B=array([[1,1j],[-2j,1]])
eig,eigv=linalg.eig(B)
print eig
print eigv

#sorting the eigensystems in the order of increasing magnitude of the eigenvaules.
A = np.random.random((3,3))
eigenValues,eigenVectors = linalg.eig(A)
print eigenValues,eigenVectors
idx = eigenValues.argsort()   
eigenValues = eigenValues[idx]
eigenVectors = eigenVectors[:,idx]
print eigenValues,eigenVectors
