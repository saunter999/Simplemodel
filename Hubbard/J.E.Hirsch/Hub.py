#!/usr/bin/env python
"""This Benchmarking codes for the well-known one_band Hubbard model serves as a (priliminary) test of the codes of dp_model 
in the determination of mu for a given filling factor,the fermi surface and the susceptibility.
QH 03/07/2015 
Reference: PhysRevB.31.4403  J.E.Hirsch.
"""
from scipy import * 
import scipy.weave as weave
import io
import os,sys,subprocess
import band 
##################################
#==    one band Hubbard Model parameters      ==#
d=2
beta=200
t=1.
occ=0.
Nk=1000
Emax=5.
Ne=500
dlt=1.0e-2
Uc=0 # Uc is set to be zero for multi-bands
##################################

#### Set MPI prefix ##############
fileMPI = 'mpi_prefix.dat'
if (os.path.exists(fileMPI)):
    print "create mpi_prefix."
    mpi_prefix = ' '.join(loadtxt(fileMPI, dtype = str))
else :
    print "mpi_prefix not found"
    mpi_prefix=''
#################################
#########0.Settingband parameters######
band.band.input_params(beta,d,Nk,Emax,Ne,dlt,occ,t)

########### 1.kmesh####################    
print "1.generating kmesh in first BZ"
band.band.kmesh()
band.band.qmesh()
band.band.ek_bz()

#sys.exit(0)
###########2.mu_fix####################

f=open("occ_mu.out",'w') 
ocup=arange(1.0,1.21,0.25)
for i,occ in enumerate(ocup):
    print i,occ
    band.band.input_params(beta,d,Nk,Ne,dlt,occ,t)
    print "2.calculating mu for nonint bands for a given occup"
    murng=[-5,5]  #Guess for the right bound for the mu(chemical potential).
    mu=band.band.mu_fix(murng)
    print>>f,occ,mu
    ########### 3.dos####################    
    print "3.calcuating dos,fermi surface and band structure along the high symmetry lines  for nonint bands"
    band.band.dos_mu(mu,i)
    band.band.band_struc(mu)
    band.band.fermi_surface(mu,i)
    band.band.suscpt_screen(mu,i)
