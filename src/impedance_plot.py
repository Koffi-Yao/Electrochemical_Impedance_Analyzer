# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 19:16:46 2019

@author: Koffi Yao
"""
import numpy as np
from matplotlib import pyplot as plt
from cmath import *
import math
font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}
plt.rc('font', **font)  # pass in the font dict as kwargs



def ZR(R, w):
    '''
    Resistor impedance
    '''
    return R + 0j

def ZC(C, w):
    '''
    impedance of capacitor
    '''
    return 0 - (1/(w*C))*1j

def warburg(Y0, w):
    '''
    '''
    return 0 + (1/Y0)*(w*1j)**(-1/2)


def series(Z1, Z2):
    '''
    Two impedances in series add up
    '''
    return Z1 + Z2

def parallel(Z1, Z2):
    '''
    Two impedances in parallel
    '''
    return 1 / (1/Z1 + 1/Z2)


if __name__ == '__main__':
    
    R1 = 20
    R2 = 100 # Ohms
    C = 1E-6 # Farad
    Y0 = 1E-6
    
    #Ang_freqs = np.linspace(0.01, 2E6, 10000)
    Ang_freqs = np.logspace(-4, 6, 10000)
    
    Z_vs_w = []
    Zprime = []
    Zprime2 = []
    Z_mag = []
    phi = []
    for w in Ang_freqs:
        '''
        Z1 = series(warburg(Y0, w), R2)
        Z2 = parallel(Z1, ZC(C, w))
        Z = series(ZR(R1, w), Z2)
        '''
        Z = warburg(Y0, w)
        
        Zprime.append(Z.real)
        Zprime2.append(-Z.imag)
        Z_mag.append(np.sqrt((Z.real)**2 + (Z.imag)**2))
        phi.append(math.atan2(-Z.imag, Z.real))
        Z_vs_w.append(Z)
    
    
    # NYQUIST PLOT
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(1,1,1)
    ax.plot(Zprime, Zprime2, 'o-', linewidth = 5)
    ax.axis('equal')
    ax.set_xlabel('Z\' (Ohm)')
    ax.set_ylabel('Z\'\' (Ohm)')
    
    # BODE PLOT
    fig1 = plt.figure(figsize = (10,10))
    ax1 = fig1.add_subplot(2,1,1)
    ax2 = fig1.add_subplot(2,1,2)
    ax1.plot(np.log10(Ang_freqs), np.log10(Z_mag), linewidth = 5)
    ax1.set_ylabel('log(|Z|)')
    ax2.plot(np.log10(Ang_freqs), np.array(phi)*(180/np.pi), linewidth = 5)
    ax2.set_ylabel('phase angle (deg)')
    ax2.set_xlabel('log(w)')
        
        
        
        

