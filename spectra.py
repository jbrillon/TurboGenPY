# -*- coding: utf-8 -*-
"""
Created on Tuesday May  29 13:08:01 2018

@author: Tony Saad
"""
# !/usr/bin/env python
from scipy import interpolate
import numpy as np
from numpy import pi

class file_spectrum:
  def __init__(self,specfname):
    filespec = np.loadtxt(specfname)
    kfile=filespec[:,0]
    efile=filespec[:,1]
    self.especf = interpolate.interp1d(kfile, efile,'cubic')
    self.kmin = kfile[0]
    self.kmax = kfile[len(kfile) - 1]
  def evaluate(self,k):
    return self.especf(k)


class cbc_spectrum:
  def __init__(self):
    cbcspec = np.loadtxt('cbc_spectrum.txt')
    kcbc=cbcspec[:,0]*100 # [1/m]
    ecbc=cbcspec[:,1]*1e-6 # [m3/s/s]

    self.especf = interpolate.interp1d(kcbc, ecbc,'cubic')
    self.kmin = kcbc[0]
    self.kmin_paper = kcbc[1] # minimum wavenumber comes from table 3 of Comte-Bellot and Corrsin's 1971 paper, corresponds to 15 1/m
    self.kmax = kcbc[len(kcbc) - 1]
  def evaluate(self,k):
    return self.especf(k)

class ml_spectrum:
  def __init__(self):
    ''' 
      Implements the Misra and Lund (1996) spectrum,
      a nondimensionalized CBC spectrum. 
      Reference: https://ntrs.nasa.gov/citations/19970014674
    '''
    cbcspec = np.loadtxt('cbc_spectrum.txt')
    # Misra and Lund non-dimensionalization:
    M = 5.08 # [cm] (mesh size from experiment)
    u_rms = 22.2 # [cm/s] (rms velocity from experiment)
    L_ref = 11.0*M/(2.0*np.pi) # cm
    U_ref = np.sqrt(3.0/2.0)*u_rms # cm/s
    energy_ref = U_ref*U_ref*L_ref # cm3/s2
    # -- nondimensionalize experiment values
    kcbc=cbcspec[:,0]*L_ref
    ecbc=cbcspec[:,1]/energy_ref

    self.especf = interpolate.interp1d(kcbc, ecbc,'cubic')
    self.kmin = kcbc[0]
    self.kmin_paper = kcbc[1] # minimum wavenumber comes from table 3 of Comte-Bellot and Corrsin's 1971 paper
    self.kmax = kcbc[len(kcbc) - 1]
  def evaluate(self,k):
    return self.especf(k)

class vkp_spectrum:
  def __init__(self,ke=40.0,nu=1.0e-5,urms=0.25):
    self.nu = nu
    self.urms = urms
    self.ke = ke
    self.kmin = 0
    self.kmax = 1e6
  def evaluate(self, k):
    ke = self.ke
    Nu = self.nu
    urms = self.urms    
    # computed from input to satisfy homogeneous turbulence properties
    Kappae = np.sqrt(5.0/12.0)*ke
    L = 0.746834/Kappae #integral length scale
    Alpha = 1.452762113;
    #  L = 0.05 # integral length scale
    #  Kappae = 0.746834/L
    Epsilon = urms*urms*urms/L;
    KappaEta = pow(Epsilon,0.25)*pow(Nu,-3.0/4.0);
    r1 = k/Kappae
    r2 = k/KappaEta
    espec = Alpha*urms*urms/ Kappae * pow(r1,4)/pow(1.0 + r1*r1,17.0/6.0)*np.exp(-2.0*r2*r2)
    return espec

class kcm_spectrum:
  def __init__(self,station=0):
    self.station_ = station
    self.epst_ = [22.8, 9.13, 4.72, 3.41]
    self.lt_  = [0.25, 0.288, 0.321, 0.332]
    self.etat_= [0.11e-3, 0.14e-3, 0.16e-3, 0.18e-3]
    self.ck_ = 1.613;
    self.alfa_ = [0.39, 1.2, 4.0, 2.1, 0.522, 10.0, 12.58]
    self.eps_ = self.epst_[station]
    self.ckEpsTwo3rd_ = self.ck_*pow(self.eps_, 0.666666666667)
    self.l_ = self.lt_[station]
    self.eta_ = self.etat_[station]
  def evaluate(self, k):
    kl = k*self.l_;
    keta = k*self.eta_;
    term1 = kl / pow( (pow(kl, self.alfa_[1]) + self.alfa_[0]), 0.83333333333333)
    term1 = pow(term1,5.66666666666667)
    term2 = 0.5 + 0.31830988618379*np.arctan(self.alfa_[5]*np.log10(keta) + self.alfa_[6])
    espec = self.ckEpsTwo3rd_*pow(k,-1.66666666666667)*term1 * np.exp(-self.alfa_[3]*keta) * (1.0 + self.alfa_[4]*term2)
    return espec
    
class pq_spectrum:
  # Implements the Passot-Pouquet Spectrum
  # http://www.wseas.us/e-library/conferences/2011/Corfu/ASM/ASM-19.pdf
  def __init__(self,ke=506, uavg=1.5):
    self.uavg = uavg
    self.ke = ke
    self.kmin = 0
    self.kmax = 1e6
  def evaluate(self, k):
    ke = self.ke
    uavg = self.uavg    
    kke = k/ke
    espec = 16.0*uavg*uavg/ke * np.sqrt(2.0/np.pi) * pow(kke,4) * np.exp(-2.0*(kke)*(kke))
    return espec
