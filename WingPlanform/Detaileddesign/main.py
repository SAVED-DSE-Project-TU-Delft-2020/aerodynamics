# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:22:56 2020

@author: halms
"""

import DATCOM_functions as dc
import numpy as np
import matplotlib.pyplot as plt


# ------------------------------- Parameters ------------------------------- #

#===========Constant Parameters===============================================

g = 9.80665                                            # [m/s^2]
rho = 1.225                                            # [kg/m^3]

#===========Varying Parameters================================================

MTOM  = 17.5                                           # [kg]
MTOW  = MTOM * g                                       # [N]
taper = 0.35                                           # [-]
cr    = 0.671                                          # [m]
ct    = cr * taper                                     # [m]
MAC   = 2/3 * cr * (1 + taper + taper**2)/(1 + taper)  # [m]
b     = 3                                              # [m]
S     = (cr + ct) * b/2                                # [m^2]
AR    = b**2/S                                         # [-]
V_cr  = 28                                             # [m/s]
LE_sw = 23.1                                           # [deg]

#===========DATCOM Parameters=================================================

C1    = dc.compute_C1(taper)                           # [-]
cond  = dc.compute_ARcondition(C1, LE_sw, AR)          # [-]

###########################         HIGH AR         ###########################



if cond == 'High AR':
    HC_sw   = dc.compute_sweep(LE_sw,taper,0.5,cr,b)   # [deg]
    CLalpha = dc.compute_CLa(AR,HC_sw)                 # [deg^-1]
    