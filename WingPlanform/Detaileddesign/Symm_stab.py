# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:51:36 2020

@author: Ihab Benyahia
"""

import numpy as np
import main as m
import DATCOM_functions as dc
import matplotlib.pyplot as plt


# Cruise parameters

CL    = m.CL
CD    = m.CD
CLa   = m.CLalpha
ind   = np.argmax(CL/CD)
a0    = m.aoa[ind]
Vcr   = m.V_cr
Vcr_u = m.V_cr_update
rho   = m.rho
S     = m.S
MTOW  = m.MTOW
AR    = m.AR
e     = m.e
MAC   = m.MAC
xcg   = 0.3043
xcgw  = 0.447

# Symmetric stability derivatives method I

CX0   = 0
CZ0   = -CL[ind]
Cm0   = 0

CXu   =
CZu   = -2*CL[ind]
Cmu   = 0

CXa   = CL[ind]*(1-CLa/(np.pi*AR*e))
CZa   = -CLa
Cma   = (xcg-xcgw)*CLa/MAC

CXq   = 0
CZq   = 0
Cmq   = 0

CXa_d = 0
CZa_d = 0
Cma_d = 0