# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:22:56 2020

@author: halms
"""

import DATCOM_functions as dc
import numpy as np
import matplotlib.pyplot as plt
import isacalculator as isa


# ------------------------------- Parameters ------------------------------- #

#===========Constant Parameters===============================================

g   = 9.80665                                          # [m/s^2]
h   = 500                                              # [m]

#===========Varying Parameters================================================

rho   = isa.compute_isa(h)[1]                          # [kg/m^3]
MTOM  = 17.536                                         # [kg]
MTOW  = MTOM * g                                       # [N]
taper = 0.35                                           # [-]
cr    = 0.696                                          # [m]
ct    = cr * taper                                     # [m]
MAC   = 2/3 * cr * (1 + taper + taper**2)/(1 + taper)  # [m]
b     = 3                                              # [m]
S     = (cr + ct) * b/2                                # [m^2]
AR    = b**2/S                                         # [-]
V_cr  = 28                                             # [m/s]
V_tr  = 14                                             # [m/s]
LE_sw = 23.1                                           # [deg]
clmax = 1.37                                           # [-]
QC_sw = dc.compute_sweep(LE_sw,taper,0.25,cr,b)        # [deg]
cldes = dc.compute_cldes(MTOM,rho,V_cr,QC_sw,S)        # [-]



#===========DATCOM Parameters=================================================

C1    = dc.compute_C1(taper)                           # [-]
cond  = dc.compute_ARcondition(C1, LE_sw, AR)          # [-]

###########################         HIGH AR         ###########################


if cond == 'High AR':
    HC_sw   = dc.compute_sweep(LE_sw,taper,0.5,cr,b)   # [deg]
    CLalpha = dc.compute_CLa(AR,HC_sw)                 # [deg^-1]
    
    
    