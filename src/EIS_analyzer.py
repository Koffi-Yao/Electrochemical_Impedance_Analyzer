
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from cmath import *
import math
from scipy.optimize import minimize
import os
from copy import deepcopy


class circuit():

    def __init__(self):
        self.elements = []
        self.name = '' # The circuit layout
        self.formula = None
        self.params = []
        
        self.__fit_data_input = []
        self.ordered_params_values = []

    def __call__(self, w):

        return None

    def function(self, w):

        def f(x):
            local = {'w': w}
            for i, n in enumerate(self.params):
                local[n] = x[i] 
                
            if not hasattr(w, '__iter__'):
                return eval(self.formula, None, local)
            else:
                ans = []
                for h in w:
                    local['w'] = h
                    ans.append(eval(self.formula, None, local))
                return ans
                
        return f
    
    def __add__(self, other):
        '''
        Two impedances in series add up
        Access using +
        '''
        result = circuit()
        result.elements = self.elements + other.elements
        result.formula = f'({self.formula} + {other.formula})'
        result.name = f'({self.name} + {other.name})'
        
        for z in result.elements:
            result.params.extend(z.params)

        return result
            
    def __floordiv__(self, other):
        '''
        Two impedances in parallel add their inverses
        Access using //
        '''
        result = circuit()
        result.elements = self.elements + other.elements
        result.formula = f'((1/{self.formula} + 1/{other.formula})**(-1))'
        result.name = f'({self.name} || {other.name})'
        for z in result.elements:
            result.params.extend(z.params)

        return result
    
    def fit(self, data_frequency, data_Z_real, data_Z_imaginary, initial_guess = None, optimizer_method = 'Nelder-Mead'):
        '''
        docs:
            data_frequency: A 1D iterable list of the frequency Omega in Hertz unit
            data_Z_real: The real part of the impedance (positive value)
            data_Z_imaginary: The imaginary part of the impedance (must be already -Zim form)
            initial_guess: the initial guess for the optimizer
            optimizer_method: The method to pass to the scipy minimizer. 
                See https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
            
            
            Output: solution

        
        '''
        
        # Convert inputs to complex data to feed the optimizer
        dataComplexZ = {'freq': [v for v in data_frequency], 'Z': [(re - 1j*im) for (re, im) in zip(data_Z_real, data_Z_imaginary)]}
        
        if initial_guess == None: 
            initial_guess = {z.name: z.value for z in self.elements}
        
        solution, x = fit_model(self, dataComplexZ, initial_guess, optimizer_method) 
        
        self.ordered_params_values = x
        self.__fit_data_input = (data_frequency, data_Z_real, data_Z_imaginary)
            
        return solution
    
    def fit_Z_values_at_frequencies(self, frequencies):
        func = self.function(frequencies)
        modelZpred = func(self.ordered_params_values)
        return modelZpred

    
    def plot_fit(self):
        
        Zpred = self.fit_Z_values_at_frequencies(self.__fit_data_input[0])
        P = np.array([[z.real,-z.imag] for z in Zpred])
        
        fig = plt.figure(figsize = (10,10))
        ax = fig.add_subplot(111)
        
        ax.plot(self.__fit_data_input[1], self.__fit_data_input[2], 'ko')
        ax.plot(P[:,0], P[:,1], '-', label = self.name)
        ax.set_xlabel('Z\'[Ohms]')
        ax.set_ylabel('Z\"[Ohms]')

        ax.set_aspect('equal', 'box')

        plt.legend()
        
        return fig
    
class Z(circuit):

    def __init__(self, **kwargs):
        '''
        kwargs: elementKind, name, value
            elementKind: 'resistor', 'capacitor', 'warburg'
            name: R1, R_eltd, Cs1, etc... A name associated with the element
            value: The value of Resistance in Ohms, Capacitance in Farad, warburg in warburg units
                !!!Note: For a warburg element, the correponding CPE Y0 is expected
                !!! This class does not consider imperfect CPE with alpha not 0, 1, or 1/2
        '''
        self.name = kwargs.get('name')
        self.value = kwargs.get('value') # CPEs take a tuple (Y0, alpha) as value
        self.elementKind = kwargs.get('elementKind')
        self.params = []
        self.Zformula()
        
    
    def Zformula(self):
        '''
        Constant phase element calculator
        '''
        self.params.append(self.name)
        self.elements = [self]

        if self.elementKind == 'resistor':
            self.formula = f'({self.name})'
        elif self.elementKind == 'capacitor':
            self.formula = f'((w*{self.name}*1j)**(-1))'
        elif self.elementKind == 'warburg':
            self.formula = f'({self.name}**(-1)*(w*1j)**(-0.5))'
        elif self.elementKind == 'CPE': # General CPE
            self.formula = f'({self.name}**(-1)*(w*1j)**(-a{self.name}))'
            self.params.append(f'a{self.name}')
        else:
            raise Exception('Unrecognized circuit element kind')

    def Zvalue(self, w):
        local = {'w': w}
        local[self.name] = self.value
        # Handling the second variable of a CPE
        if f'a{self.name}' in self.params:
            local[self.name] = self.value[0]
            local[f'a{self.name}'] = self.value[1]
        return eval(self.formula, None, local)




# Free standing functions
def targetConstructor(Zequiv, dataComplexZ):
    '''
    Zeq: The Zeq model built from the Zeq class
    dataComplexZ: The dataset as a dictionary of frequencies (freq) and and complex impedance values (Z)
    '''
    W = dataComplexZ['freq']
    func = Zequiv.function(W)
    def target(x):
        error = 0
        modelZ = func(x)
        for i in range(len(W)):
            diff = modelZ[i] - dataComplexZ['Z'][i]
            error += (diff.real)**2 + (diff.imag)**2
        return (error**(1/2))/len(W)
    return target


def fit_model(Zequiv, dataComplexZ, initGuess, method):
    '''
    initGuess: initial guess dictionary of the parameters
    '''
    # Create target function
    params = Zequiv.params
    X0 = []
    for p in params:
        if 'CPE' not in p:
            X0.append(initGuess[p])
        else:
            if 'aCPE' in p:
                X0.append(initGuess[p[1:]][1])
            else: 
                X0.append(initGuess[p][0])
        
    
    solution = {p: X0[i] for i, p in enumerate(params)}
    solution['success'] = False
    target = targetConstructor(Zequiv, dataComplexZ)
    result = minimize(target, X0, method = method, options = {'maxiter': 20000, 'disp': True})
    x = X0
    #print(result.message)
    if result.success:
        solution['success'] = True
        x = result.x
        for i, p in enumerate(params):
            solution[p] = x[i]
    return solution, x
        
def unpackComplex(z):
    
    return z.real, z.imag



    