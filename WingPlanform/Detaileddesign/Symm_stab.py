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
Vcr   = m.V_cr
Vcr_u = m.V_cr_update
rho   = m.rho
S     = m.S
MTOW  = m.MTOW
AR    = m.AR
e     = m.e

# Symmetric stability derivatives method I

CX0   = -CD[ind]
CZ0   = -CL[ind]
Cm0   = 0

CXu   = 2*CL[ind]*np.tan(np.arcsin(CD[ind]*0.5*rho*Vcr_u**2*S/MTOW))
CZu   = -2*CL[ind]
Cmu   = 0

CXa   = CL*(1-CLa/(np.pi*AR*e))
CZa   = -CLa

