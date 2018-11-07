#!/usr/bin/env python
from scipy import *

####python uses pass by reference for array and pass by value for number.
def func(a):
	a[0]=3;a[1]=4;
	return a



if __name__=="__main__":
	b=array([1,2]);
	c=func(b);
	print b,c
