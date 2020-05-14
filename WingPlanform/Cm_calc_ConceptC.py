# -*- coding: utf-8 -*-
"""
Created on Thu May 14 12:22:43 2020

@author: Ihab Benyahia
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt("ConceptC_T1-28_0 m_s-LLT.txt",skiprows=1)

def linear(x,a,b):
    return a*x+b
alpha = data[:,0]
Cm = data[:,8]

coeff,covar = curve_fit(linear,alpha,Cm)
Cmalpha,Cm0 = coeff

plt.plot(alpha,Cm,alpha,linear(alpha,Cmalpha,Cm0))

print(f'Cmalpha = {Cmalpha} [deg^-1]')
print(f'Cm0 = {Cm0} [-]')