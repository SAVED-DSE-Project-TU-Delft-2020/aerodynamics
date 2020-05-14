# -*- coding: utf-8 -*-
"""
Created on Wed May 13 11:26:35 2020

@author: halms
"""


import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

AR = 9/1.3

data_14 = np.loadtxt("ConceptA_D_T1-14_0 m_s-LLT.txt",skiprows=1)

alpha_14 = data_14[:,0]
CL_14 = data_14[:,2]
CD_14 = data_14[:,5]
CDi_14 = data_14[:,3]

e_14 = CL_14**2/(np.pi*AR*CDi_14)

def linear(x,a,b):
    return a*x + b

a,b = curve_fit(linear,alpha_14,CL_14)[0]


data_28 = np.loadtxt("ConceptA_D_T1-28_0 m_s-LLT.txt",skiprows=1)

alpha_28 = data_28[:,0]
CL_28 = data_28[:,2]
CD_28 = data_28[:,5]
CDi_28 = data_28[:,3]

def polarf(x,a,c):
    return a*x**2 + c

coeff_28 = curve_fit(polarf,CL_28,CD_28)[0]


CD0_28 = coeff_28[1]
e_28 = CL_28**2/(np.pi*AR*CDi_28)

#%%


