# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:51:15 2020

@author: Ihab Benyahia
"""


import numpy as np
import matplotlib.pyplot as plt

Fx_alpha = np.loadtxt("Fx_alpha.txt",skiprows = 1)
Fz_alpha = np.loadtxt("Fz_alpha.txt",skiprows = 1)
Cm_alpha = np.loadtxt("Cm_alpha.txt",skiprows = 1)
Fx_V     = np.loadtxt("Fx_V.txt",skiprows = 1)
Fz_V     = np.loadtxt("Fz_V.txt",skiprows = 1)
Cm_V     = np.loadtxt("Cm_V.txt",skiprows = 1)
L_beta   = np.loadtxt("L_beta.txt",skiprows = 1)
N_beta   = np.loadtxt("N_beta.txt",skiprows = 1)


beta  = L_beta[:,0]
L     = L_beta[:,1]
N     = N_beta[:,1]


alpha    = Fx_alpha[:,0]
X_alpha  = Fx_alpha[:,1]
Z_alpha  = Fz_alpha[:,1]
Cm_alpha = Cm_alpha[:,1]

V        = Fx_V[:,0]
X_V      = Fx_V[:,1]
Z_V      = Fz_V[:,1]
Cm_V     = Cm_V[:,1]


def central(x,y):
    res = np.empty(len(x))
    res[0] = (y[1] - y[0])/(x[1] - x[0])
    res[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])
    for i in range(1,len(x)-1):
        res[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
    return res



# Xu_xflr = -0.57771
# Xa_xflr = 

# Zu_xflr = 0
# Za_xflr = 0

# Cmu_xflr = 0
# Cma_xflr = 0