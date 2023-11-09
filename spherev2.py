#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 18:44:03 2023

@author: rodrigo
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
#import os

def spherev2(radius,jmolout = False, ddaout = False):
    '''
    Function that generates points forming a spherical particle in a square 
    computational grid, specifically in the positive octant.

    Parameters
    ----------
    radius : integer
        Radius of the sphere in the computational grid. 
        In principle it can accept any positive integer, but values > 5 are preferred
        
    jmolout : Boolean, optional
        To write a .xyz file for jmol. The default is False.
        
    ddaout : Boolean, optional
        To write a .tgt file for OpenDDA. The default is False.

    Returns
    -------
    None.

    '''
        
    radius = int(radius)
    input_radius = radius
    Nangles = 1000
    
    theta_original = np.linspace(0, pi, Nangles, endpoint = True)
    phi_original = np.linspace(0, 2*pi, Nangles)
    radius_original = np.linspace(0,input_radius,input_radius + 1, endpoint = True)
    
    [radius,theta,phi] = np.meshgrid(radius_original,theta_original,phi_original)
    
    # Make all arrays one dimensional
    radius = radius.flatten()
    theta = theta.flatten()
    phi = phi.flatten()
    
    # Parametric equations for a sphere
    xval = radius * np.sin(theta) * np.cos(phi)
    yval = radius * np.sin(theta) * np.sin(phi)
    zval = radius * np.cos(theta)
    
    # Shift the points towards the positive octant
    xval = xval - min(xval) + 1
    yval = yval - min(yval) + 1
    zval = zval - min(zval) + 1
    
    Npoints = len(zval)
    print(f'Number of points {Npoints}')
    
    # Convert all arrays to integers
    xval = xval.astype(int)
    yval = yval.astype(int)
    zval = zval.astype(int)
    
    xmax = max(xval)
    ymax = max(yval)
    zmax = max(zval)
    
    occupied_points = np.zeros((xmax + 1,ymax + 1,zmax + 1))
    
    for k in range(0,Npoints):
        occupied_points[xval[k],yval[k],zval[k]] = 1
    
    xval = yval = zval = []
    
    # Retrieve xyz coordinates that are occupied
    x_new, y_new, z_new = np.where(occupied_points == 1)
    
    occupied_points = []
    
    # Make a scatter plot
    fig = plt.figure(1)
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(x_new,y_new,z_new)
    
    Npoints = len(z_new)
    print(f'Number of unique points {Npoints}')
       
    # Write files with locations
    script_name = 'spherev2.py'
    
    if jmolout:
        
        filename = f'sphere_R{input_radius}.xyz'
        f = open(filename,'w')
        f.write(f'{Npoints}\n')
        f.write(f'Sphere generated with {script_name}\n')
        
        for k in range(0,Npoints):
            f.write('Si {} {} {}\n'.format(x_new[k],y_new[k],z_new[k]))
        
        f.close()
        
        print(f'{filename} written')
    
    if ddaout:
        
        filename = f'sphere_R{input_radius}.tgt'
        f = open(filename,'w')
        f.write(f'# Sphere generated with {script_name}\n')
        f.write(f'# Radius {input_radius}\n')
        f.write(f'# Number of dipoles {Npoints}\n')
        
        for k in range(0,Npoints):
            f.write('{},{},{},0,0,0\n'.format(x_new[k],y_new[k],z_new[k]))
        
        f.close()
        
        print(f'{filename} written')
        