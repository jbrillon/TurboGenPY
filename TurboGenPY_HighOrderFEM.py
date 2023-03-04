# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 17:12:30 2023

@author: Julien Brillon
"""
# !/usr/bin/env python
from scipy import interpolate
from scipy import integrate
import numpy as np
from numpy import pi
import time
import scipy.io
from tkespec import compute_tke_spectrum
import isoturb

import spectra

# ----------------------------------------------------------------------------------------------
#  
#  
#    _  _ ___ ___ _  _    ___  ___ ___  ___ ___   ___ ___ __  __   _____ _   _ ___ ___  ___  
#   | || |_ _/ __| || |  / _ \| _ \   \| __| _ \ | __| __|  \/  | |_   _| | | | _ \ _ )/ _ \ 
#   | __ || | (_ | __ | | (_) |   / |) | _||   / | _|| _|| |\/| |   | | | |_| |   / _ \ (_) |
#   |_||_|___\___|_||_|  \___/|_|_\___/|___|_|_\ |_| |___|_|  |_|   |_|  \___/|_|_\___/\___/ 
#                                                                                            
#  
#  
# ----------------------------------------------------------------------------------------------

def generate_isotropic_turbulence_high_order_fem(
    number_of_elements_per_direction,
    poly_degree,
    output_filename="velocity_equidistant_nodes.fld",
    number_of_modes=5000, # suggested in TurboGenPY paper
    spectra_name='ml'):

    number_of_unique_points_per_direction = number_of_elements_per_direction*(poly_degree+1) - (number_of_elements_per_direction-1)

    nx = 1*number_of_unique_points_per_direction
    ny = 1*number_of_unique_points_per_direction
    nz = 1*number_of_unique_points_per_direction

    # specify which spectrum you want to use
    # inputspec = 'cbc'
    inputspec = str(spectra_name)

    # Default values for domain size in the x, y, and z directions. This value is typically
    # based on the largest length scale that your data has. For the cbc data,
    # the largest length scale corresponds to a wave number of 15, hence, the
    # domain size is L = 2pi/15.
    if(inputspec == 'cbc'):
        lx = 9 * 2.0 * pi / 100.0 # [m]
        ly = 9 * 2.0 * pi / 100.0 # [m]
        lz = 9 * 2.0 * pi / 100.0 # [m]
    elif(inputspec == 'ml'):
        lx = 2.0 * pi
        ly = 2.0 * pi
        lz = 2.0 * pi

    # number of modes
    nmodes = 1.0*number_of_modes

    # specify the spectrum name to append to all output filenames
    fileappend = inputspec + '_' + str(nx) + '.' + str(ny) + '.' + str(nz) + '_' + str(nmodes) + '_modes'

    print('input spec', inputspec)
    if inputspec != 'cbc' and inputspec != 'ml' and inputspec != 'vkp' and inputspec != 'kcm' and inputspec != 'pq':
        print('Error: ', inputspec, ' is not a supported spectrum. Supported spectra are: cbc, ml, vkp, kcm, and pq. Please revise your input.')
        exit()
    inputspec += '_spectrum'
    # now given a string name of the spectrum, find the corresponding function with the same name. use locals() because spectrum functions are defined in this module.
    # whichspec = locals()[inputspec]
    # whichspec = spectra.cbc_spectrum().evaluate
    whichspec = getattr(spectra, inputspec)().evaluate

    # smallest wavenumber that can be represented by this grid
    wn1_grid = min(2.0*pi/lx, min(2.0*pi/ly, 2.0*pi/lz))
    # enter the smallest wavenumber represented by this spectrum
    if(inputspec=='cbc_spectrum' or inputspec=='ml_spectrum'):
        # NOTE: this 15 [1/m] minimum wavenumber comes from table 3 of comte's original (CBC) paper
        # determined here from cbc spectrum properties (suggested in TurboGenPY paper)
        wn1 = getattr(spectra, inputspec)().kmin_paper
    else:
        # default to grid wavenumber
        wn1 = 1.0*wn1_grid

    # summarize user input
    print('-----------------------------------')
    print('SUMMARY OF USER INPUT (TurboGenPY):')
    print('Domain size:', lx, ly, lz)
    print('Grid resolution:', nx, ny, nz)
    print('Fourier accuracy (modes):', nmodes)
    print('Smallest wavenumber represented by this spectrum: %.3f' % wn1)
    print('Smallest wavenumber represented by this grid: %.3f' % wn1_grid)
    print('Which spectra: %s' % inputspec)
    print('Poly degree: ', poly_degree)
    print('Number of elements per direction: ', number_of_elements_per_direction)
    print('DOFs per dim: ', ((poly_degree+1.0)*number_of_elements_per_direction))
    print('Output filename: %s' % output_filename)

    # ------------------------------------------------------------------------------
    # END USER INPUT
    # ------------------------------------------------------------------------------

    t0 = time.time()
    # generate isotropic turbulence and write to file
    u, v, w = isoturb.generate_isotropic_turbulence(lx, ly, lz, nx, ny, nz, nmodes, wn1, whichspec, 
        write_fem_grid_files=True,
        number_of_elements_per_direction=number_of_elements_per_direction,
        poly_degree=poly_degree,
        velocity_field_filename=output_filename)
    t1 = time.time()
    elapsed_time = t1 - t0
    print('it took me ', elapsed_time, 's to generate the isotropic turbulence.')
    print('done.')
    print('-----------------------------------')
    return
