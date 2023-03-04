#
#  isoturb.py
#
#  The MIT License (MIT)
#
#  Copyright (c) 2015, Tony Saad. All rights reserved.
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
#

# -*- coding: utf-8 -*-
"""
Created on Mon May 12 09:31:54 2014

@author: tsaad
"""
import numpy as np
from numpy import sin, cos, sqrt, ones, zeros, pi, arange


def generate_isotropic_turbulence(lx, ly, lz, nx, ny, nz, nmodes, wn1, especf, 
  write_fem_grid_files=False,
  number_of_elements_per_direction=4,
  poly_degree=5,
  velocity_field_filename="velocity_equidistant_nodes.fld"):
    """
    Given an energy spectrum, this function computes a discrete, staggered, three
    dimensional velocity field in a box whose energy spectrum corresponds to the input energy
    spectrum up to the Nyquist limit dictated by the grid

    This function returns u, v, w as the axial, transverse, and azimuthal velocities.

    Parameters:
    -----------
    lx: float
      The domain size in the x-direction.
    ly: float
      The domain size in the y-direction.
    lz: float
      The domain size in the z-direction.
    nx: integer
      The number of grid points in the x-direction.
    ny: integer
      The number of grid points in the y-direction.
    nz: integer
      The number of grid points in the z-direction.
    wn1: float
      Smallest wavenumber. Typically dictated by spectrum or domain size.
    espec: functor
      A callback function representing the energy spectrum.
    """

    if(write_fem_grid_files==True):
      dx = lx / (nx-1)
      dy = ly / (ny-1)
      dz = lz / (nz-1)
    else:
      # generate cell centered x-grid
      dx = lx / nx
      dy = ly / ny
      dz = lz / nz

    # START THE FUN!

    # compute random angles
    phi = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes)
    nu = np.random.uniform(0.0, 1.0, nmodes)
    theta = np.arccos(2.0 * nu - 1.0)
    psi = np.random.uniform(-pi / 2.0, pi / 2.0, nmodes)

    # highest wave number that can be represented on this grid (nyquist limit)
    wnn = max(np.pi / dx, max(np.pi / dy, np.pi / dz))
    print('I will generate data up to wave number: ', wnn)

    # wavenumber step
    dk = (wnn - wn1) / nmodes

    if(write_fem_grid_files==True):
      # wavenumber as described in TurboGenPY paper (step 4)
      wn = wn1 + arange(0, nmodes) * dk
    else:
      # wavenumber at cell centers
      wn = wn1 + 0.5 * dk + arange(0, nmodes) * dk

    dkn = ones(nmodes) * dk

    #   wavenumber vector from random angles
    kx = sin(theta) * cos(phi) * wn
    ky = sin(theta) * sin(phi) * wn
    kz = cos(theta) * wn

    # create divergence vector -- pretty sure this comes from eq.(7) without the 2/delta factor at the front
    ktx = np.sin(kx * dx / 2.0) / dx
    kty = np.sin(ky * dy / 2.0) / dy
    ktz = np.sin(kz * dz / 2.0) / dz

    # Enforce Mass Conservation
    phi1 = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes)
    nu1 = np.random.uniform(0.0, 1.0, nmodes)
    theta1 = np.arccos(2.0 * nu1 - 1.0)
    zetax = sin(theta1) * cos(phi1)
    zetay = sin(theta1) * sin(phi1)
    zetaz = cos(theta1)
    sxm = zetay * ktz - zetaz * kty
    sym = -(zetax * ktz - zetaz * ktx)
    szm = zetax * kty - zetay * ktx
    smag = sqrt(sxm * sxm + sym * sym + szm * szm)
    sxm = sxm / smag
    sym = sym / smag
    szm = szm / smag

    # verify that the wave vector and sigma are perpendicular
    kk = np.sum(ktx * sxm + kty * sym + ktz * szm)
    print('Orthogonality of k and sigma (divergence in wave space):', kk)

    # get the modes
    km = wn

    espec = especf(km)
    espec = espec.clip(0.0)

    
    um = sqrt(espec * dkn)
    u_ = zeros([nx, ny, nz])
    v_ = zeros([nx, ny, nz])
    w_ = zeros([nx, ny, nz])

    # generate turbulence at cell centers
    xc = dx / 2.0 + arange(0, nx) * dx
    yc = dy / 2.0 + arange(0, ny) * dy
    zc = dz / 2.0 + arange(0, nz) * dz

    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # for every grid point (i,j,k) do the fourier summation
                arg = kx * xc[i] + ky * yc[j] + kz * zc[k] - psi
                bmx = 2.0 * um * cos(arg - kx * dx / 2.0)
                bmy = 2.0 * um * cos(arg - ky * dy / 2.0)
                bmz = 2.0 * um * cos(arg - kz * dz / 2.0)
                u_[i, j, k] = np.sum(bmx * sxm)
                v_[i, j, k] = np.sum(bmy * sym)
                w_[i, j, k] = np.sum(bmz * szm)

    if(write_fem_grid_files==True):
      # generate files for PHiLiP
      # IDEA: I could compute the velocity field directly at the GLL nodes instead of having to convert it after

      number_of_points_per_direction = number_of_elements_per_direction*(poly_degree+1)
      file = open(velocity_field_filename,"w") # for testing
      # - write total DOFs
      wstr = "%i\n" % (number_of_points_per_direction*number_of_points_per_direction*number_of_points_per_direction)
      file.write(wstr)

      nElements_per_direction = 1*number_of_elements_per_direction
      nQuadPoints_per_element = poly_degree+1
      x_element_faces = np.linspace(0,lx,nElements_per_direction+1)
      y_element_faces = np.linspace(0,ly,nElements_per_direction+1)
      z_element_faces = np.linspace(0,lz,nElements_per_direction+1)
      
      i = 0
      for ez in range(0,nElements_per_direction):
        for qz in range(0,nQuadPoints_per_element):
          z_val = z_element_faces[ez] +  np.float64(qz) * dz # z-coordinate
          for ey in range(0,nElements_per_direction):
            for qy in range(0,nQuadPoints_per_element):
              y_val = y_element_faces[ey] +  np.float64(qy) * dy # y-coordinate
              for ex in range(0,nElements_per_direction):
                for qx in range(0,nQuadPoints_per_element):
                  x_val = x_element_faces[ex] +  np.float64(qx) * dx # x-coordinate

                  # for every grid point (i,j,k) do the fourier summation
                  arg = kx * x_val + ky * y_val + kz * z_val - psi
                  bmx = 2.0 * um * cos(arg - kx * dx / 2.0)
                  bmy = 2.0 * um * cos(arg - ky * dy / 2.0)
                  bmz = 2.0 * um * cos(arg - kz * dz / 2.0)
                  u_val = np.sum(bmx * sxm)
                  v_val = np.sum(bmy * sym)
                  w_val = np.sum(bmz * szm)

                  wstr = " %18.16e %18.16e %18.16e %18.16e %18.16e %18.16e\n" % (x_val, y_val, z_val, u_val, v_val, w_val)
                  file.write(wstr)
                  i += 1

      file.close()

    print('done generating turbulence.')
    return u_, v_, w_


def generate_scalar_isotropic_turbulence(lx, ly, lz, nx, ny, nz, nmodes, wn1, especf):
    """
    Given an energy spectrum, this function computes a discrete, staggered, three
    dimensional velocity field in a box whose energy spectrum corresponds to the input energy
    spectrum up to the Nyquist limit dictated by the grid

    This function returns u, v, w as the axial, transverse, and azimuthal velocities.

    Parameters:
    -----------
    lx: float
      The domain size in the x-direction.
    ly: float
      The domain size in the y-direction.
    lz: float
      The domain size in the z-direction.
    nx: integer
      The number of grid points in the x-direction.
    ny: integer
      The number of grid points in the y-direction.
    nz: integer
      The number of grid points in the z-direction.
    wn1: float
      Smallest wavenumber. Typically dictated by spectrum or domain size.
    espec: functor
      A callback function representing the energy spectrum.
    """

    # generate cell centered x-grid
    dx = lx / nx
    dy = ly / ny
    dz = lz / nz

    # START THE FUN!

    # compute random angles
    phi = 2.0 * pi * np.random.uniform(0.0, 1.0, nmodes)
    nu = np.random.uniform(0.0, 1.0, nmodes);
    theta = np.arccos(2.0 * nu - 1.0);
    psi = np.random.uniform(-pi / 2.0, pi / 2.0, nmodes)

    # highest wave number that can be represented on this grid (nyquist limit)
    wnn = max(np.pi / dx, max(np.pi / dy, np.pi / dz))
    print('I will generate data up to wave number: ', wnn)

    # wavenumber step
    dk = (wnn - wn1) / nmodes

    # wavenumber at cell centers
    wn = wn1 + 0.5 * dk + arange(0, nmodes) * dk

    dkn = ones(nmodes) * dk

    #   wavenumber vector from random angles
    kx = sin(theta) * cos(phi) * wn
    ky = sin(theta) * sin(phi) * wn
    kz = cos(theta) * wn

    # get the modes
    km = wn

    espec = especf(km)
    espec = espec.clip(0.0)

    # generate turbulence at cell centers
    um = sqrt(espec * dkn)
    scalar_ = zeros([nx, ny, nz])

    xc = dx / 2.0 + arange(0, nx) * dx
    yc = dy / 2.0 + arange(0, ny) * dy
    zc = dz / 2.0 + arange(0, nz) * dz

    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # for every grid point (i,j,k) do the fourier summation
                arg = kx * xc[i] + ky * yc[j] + kz * zc[k] - psi
                bm = 2.0 * um * cos(arg)
                scalar_[i, j, k] = np.sum(bm)

                print('done. I am awesome!')
    return scalar_

def generate_isotropic_turbulence_2d(lx, ly, nx, ny, nmodes, wn1, especf):
    """
    This is the 2D version of the isotropic turbulence generator.
    Given an energy spectrum, this function computes a discrete, staggered, three
    dimensional velocity field in a box whose energy spectrum corresponds to the input energy
    spectrum up to the Nyquist limit dictated by the grid

    This function returns u, v as the axial and transverse velocities.

    Parameters:
    -----------
    lx: float
      The domain size in the x-direction.
    ly: float
      The domain size in the y-direction.
    nx: integer
      The number of grid points in the x-direction.
    ny: integer
      The number of grid points in the y-direction.
    wn1: float
      Smallest wavenumber. Typically dictated by spectrum or domain size.
    espec: functor
      A callback function representing the energy spectrum.
    """

    # generate cell centered x-grid
    dx = lx / nx
    dy = ly / ny

    # START THE FUN!

    # compute random angles
    psi = np.random.uniform(-pi / 2.0, pi / 2.0, nmodes)

    # highest wave number that can be represented on this grid (nyquist limit)
    wnn = max(np.pi / dx, np.pi / dy)
    print('I will generate data up to wave number: ', wnn)

    # wavenumber step
    dk = (wnn - wn1) / nmodes

    # wavenumber at cell centers
    wn = wn1 + 0.5 * dk + arange(0, nmodes) * dk

    dkn = ones(nmodes) * dk

    #   wavenumber vector from random angles
    theta = np.random.uniform(0.0,2.0*np.pi,nmodes)
    kx = cos(theta) * wn
    ky = sin(theta) * wn

    # create divergence vector
    ktx = np.sin(kx * dx / 2.0) / dx
    kty = np.sin(ky * dy / 2.0) / dy

    # Enforce Mass Conservation
    sxm = -kty
    sym = ktx
    
    smag = sqrt(sxm * sxm + sym * sym)
    sxm = sxm / smag
    sym = sym / smag

    # verify that the wave vector and sigma are perpendicular
    kk = np.sum(ktx * sxm + kty * sym)
    print('Orthogonality of k and sigma (divergence in wave space):', kk)

    # get the modes
    km = wn

    espec = especf(km)
    espec = espec.clip(0.0)

    # generate turbulence at cell centers
    um = sqrt(espec * dkn)
    u_ = zeros([nx, ny])
    v_ = zeros([nx, ny])

    xc = dx / 2.0 + arange(0, nx) * dx
    yc = dy / 2.0 + arange(0, ny) * dy
    for j in range(0, ny):
    	for i in range(0, nx):
    		# for every grid point (i,j,k) do the fourier summation
    		arg = kx * xc[i] + ky * yc[j]  - psi
    		bmx = 2.0 * um * cos(arg - kx * dx / 2.0)
    		bmy = 2.0 * um * cos(arg - ky * dy / 2.0)
    		u_[i, j] = np.sum(bmx * sxm)
    		v_[i, j] = np.sum(bmy * sym)

    print('done generating turbulence.')
    return u_, v_