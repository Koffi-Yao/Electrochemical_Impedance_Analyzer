
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 19:16:46 2019

@author: Koffi Yao
"""
import numpy as np
from matplotlib import pyplot as plt
from cmath import *
import math
from scipy.optimize import minimize
import os

from EIS_analyzer import *

if __name__ == "__main__":
    
    folder = os.path.dirname(os.getcwd())
    folder = os.path.join(folder, 'resources')
    fname = r'DEMO EIS data.csv'
    fpath = os.path.join(folder, fname)
    data = np.genfromtxt(fpath, delimiter = ',', skip_header = 1)
    
    
    Relt = Z(name = 'Relt', value = 100, elementKind = 'resistor')
    Rct = Z(name = 'Rct', value = 400, elementKind = 'resistor')
    Celd = Z(name = 'Celd', value = 1E-4, elementKind = 'capacitor')
    Wd = Z(name = 'Wd', value = 1E-2, elementKind = 'warburg')
    CPEeld = Z(name = 'CPEeld', value = (1E-2,0.7), elementKind = 'CPE')
    
    # BUILD THE CIRCUIT FROM SOME ELEMENTS
    #Ze = Relt + ((Rct + Wd) // Celd)
    Ze = Relt + ((Rct + Wd) // CPEeld)
    
    solution = Ze.fit(data[:,0], data[:,1], data[:,2])
    Ze.plot_fit()
    
    print(f'Solution = {solution}')
    
    





# =============================================================================
# if __name__ == "__main__":
#     #Ang_freqs = np.logspace(-4, 6, 100)
#     folder = r'G:\.shortcut-targets-by-id\1BnxvHTeLpN2mitgFkKSqAAZhrBfKHgNB\UD\CLASSES TAUGHT\MEEG 467-667 DIAGNOSTICS OF BATTERIES\Fall 2019\Class demo codes\EIS code\resources'
#     fname = r'DEMO EIS data.csv'
#     fpath = os.path.join(folder, fname)
#     data = np.genfromtxt(fpath, delimiter = ',', skip_header = 1)
#     
#     dataComplexZ = {'freq': data[:,0].tolist(), 'Z': []}
#     for i in range(len(dataComplexZ['freq'])):
#         dataComplexZ['Z'].append(data[i,1] - 1j*data[i,2])
#     
#     Relt = Z(name = 'Relt', value = 100, elementKind = 'resistor')
#     Rct = Z(name = 'Rct', value = 400, elementKind = 'resistor')
#     Celd = Z(name = 'Celd', value = 1E-4, elementKind = 'capacitor')
#     Wd = Z(name = 'Wd', value = 1E-2, elementKind = 'warburg')
#     #Cb = Z(name = 'Cb', value = 1E-2, elementKind = 'capacitor')
#     CPEeld = Z(name = 'CPEeld', value = (1E-2,0.7), elementKind = 'CPE')
#     
#     #Ze = Relt + ((Rct + Wd) // Celd)
#     Ze = Relt + ((Rct + Wd) // CPEeld)
#     
#     initGuess = {z.name: z.value for z in Ze.elements}
#     solution, x = fit_model(Ze, dataComplexZ, initGuess, 'Nelder-Mead')
#     
#     W = dataComplexZ['freq']
#     func = Ze.function(W)
#     modelZpred = func(x)
#     P = np.array([[z.real,-z.imag] for z in modelZpred])
#     
#     fig = plt.figure(figsize = (10,10))
#     ax = fig.add_subplot(111)
#     
#     ax.plot(data[:,1], data[:,2], 'ko')
#     ax.plot(P[:,0], P[:,1], '-', label = Ze.name)
#     ax.set_xlabel('Z\'[Ohms]')
#     ax.set_ylabel('Z\"[Ohms]')
#     ax.set_xlim(0,2000)
#     ax.set_ylim(0, 1000)
#     ax.set_aspect('equal', 'box')
# 
#     plt.legend()
#     
#     print(f'Solution = {solution}')
# =============================================================================
