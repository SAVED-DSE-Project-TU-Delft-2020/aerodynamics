# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:22:56 2020

@author: halms
"""

import DATCOM_functions as dc
import numpy as np
import matplotlib.pyplot as plt
import isacalculator as isa
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


# ------------------------------- Parameters ------------------------------- #

#===========Constant Parameters===============================================

g       = 9.80665                                        # [m/s^2]
h       = 500                                            # [m]
twist   = False                                          # [-]
#===========Varying Parameters================================================

rho     = isa.compute_isa(h)[1]                          # [kg/m^3]
MTOM    = 17.536                                         # [kg]
MTOW    = MTOM * g                                       # [N]
taper   = 0.35                                           # [-]
cr      = 0.696                                          # [m]
ct      = cr * taper                                     # [m]
MAC     = 2/3 * cr * (1 + taper + taper**2)/(1 + taper)  # [m]
b       = 3                                              # [m]
S       = (cr + ct) * b/2                                # [m^2]
AR      = b**2/S                                         # [-]
V_cr    = 28                                             # [m/s]
V_tr    = 14                                             # [m/s]
LE_sw   = 23.1                                           # [deg]
clmax   = 1.37                                           # [-]
QC_sw   = dc.compute_sweep(LE_sw,taper,0.25,cr,b)        # [deg]
cldes   = dc.compute_cldes(MTOM,rho,V_cr,QC_sw,S)        # [-]
clalpha = 0.11067579                                     # [deg^-1]
cl0     = 0.0074437                                      # [-]

#===========DATCOM Parameters=================================================

C1      = dc.compute_C1(taper)                           # [-]
C2      = 1                                              # [-]
cond    = dc.compute_ARcondition(C1, LE_sw, AR)          # [-]

###########################         HIGH AR         ###########################


if cond == 'High AR':
    
    HC_sw = dc.compute_sweep(LE_sw,taper,0.5,cr,b)                   # [deg]
    CLalpha = dc.compute_CLa(AR,HC_sw)                               # [rad^-1]
    
    if twist == False:
        aoa_des = (cldes - cl0)/clalpha                              # [deg]
        aoa_0 = dc.compute_AOA_0_lift(clalpha,cldes,aoa_des)         # [deg]
    
    else:
        aoa_zero = -cl0 / clalpha                                    # [deg]
        theta =  ___                                                 # [deg]
        aoa0_theta = ___                                             # [-]
        aoa_0 = dc.compute_AOA_0_lift_tw(aoa0_theta,theta,aoa_zero)  # [deg]
        
    delta_y = 3.72                                                   # [-]
    CL_cl = 0.83                                                     # [-]
    delta_CLmax = 0                                                  # [-]
    CLmax = dc.compute_CLmax_high(CL_cl,clmax,delta_CLmax)           # [-]
    
    delta_aoa = 2.4                                                  # [deg]
    aoa_stall = dc.compute_aoa_stall_high(CLmax,CLalpha,aoa_0,delta_aoa) # [deg]
    
    CN_prime_max = dc.compute_CN_prime_CLmax(CLmax,aoa_stall)        # [-]
    J = dc.compute_Jpar(C1,C2,LE_sw,AR)                              # [-]
    aoa = np.linspace(aoa_0,aoa_stall,100)                           # [deg]
    tan_ratio = np.tan(aoa*np.pi/180) / np.tan(aoa_stall*np.pi/180)  # [-]
    delta_CNaa = np.empty(len(aoa))                                  # [deg^-1]
    
    for i in range(len(aoa)):
        if tan_ratio[i] <= .6:
            delta_CNaa[i] = 2.8
        else:
            delta_CNaa[i] = dc.compute_linear(0.6,1,2.8,0,tan_ratio[i])
            
    CNaa_ref = dc.compute_CNaa_ref(CN_prime_max,CLalpha,aoa_stall)   # [-]
    CNaa_below = dc.compute_CNaa_below(CNaa_ref,delta_CNaa)          # [-]
    CN_prime = dc.compute_CN_prime(CLalpha, aoa, CNaa_below)         # [-]
    CL_below = dc.compute_CL(CN_prime,aoa)                           # [-]
    
#     daoa = aoa[1] - aoa[0]                                           # [deg]
#     aoa_aft = np.arange(aoa_stall+daoa,30+daoa,daoa)               # [deg]
#     CNaa_90 = 1.3                                                    # [-]
#     tan_ratio_rev = np.tan(aoa_stall*np.pi/180) / np.tan(aoa_aft*np.pi/180) # [-]
    
    tc_avg = 0.08                                                  # [-]
    x_tcmax = 0.2494                                                 # [-]
    Sref = S                                                         # [m^2]
    Swet = 2*S                                                       # [m^2]
    Cf = 0.0046                                                      # [-]
    TMAX_sw = dc.compute_sweep(LE_sw,taper,x_tcmax,cr,b)             # [-]
    Rls = 1.06                                                       # [-]
    CD0_wing = dc.compute_CD0_wing(Cf,tc_avg,x_tcmax,Swet,Sref,Rls)  # [-]
    
    
            
# def ra(x,a,b,c,d):
#     return a*np.sin(b*x+c) + d

# def ra2(x,a,b,c):
#     return a*x**2 + b*x + c

# xdata = np.array([1,.87,.8,.65,.6])
# ydata = [0,-0.5,-1,-1.5,-1.6]
# a,b,c,d = curve_fit(ra,xdata,ydata)[0]
# x = np.linspace(.6,1,100)
# y = ra(x,a,b,c,d)

# udata = np.array([.6,.4,.2,0])
# vdata = [-1.6,-1.25,-0.7,0]
# e,f,g = curve_fit(ra2,udata,vdata)[0]
# u = np.linspace(0,.6,100)
# v = ra2(u,e,f,g)

# plt.plot(x,y)
# plt.plot(u,v)


# --------------------------- Extra Computations --------------------------- #
# This is done to determine the average thickness

data = np.loadtxt("CAL4014L.dat")
res = [data[:,1][i] for i in range(len(data[:,0]))]
neg = []
pos = []
for i in res:
    if i<0:
        neg.append(i)
    else:
        pos.append(i)

neg = np.array(neg)
pos = np.array(pos)

lists = []
ilist = []
maxlist = []
while len(pos) >= 1:
    for i in range(len(pos)):
        lists.append(np.max(pos[i] - neg))
    maxlist.append(np.max(lists))
    pos = np.delete(pos,np.argmax(lists))
    lists = []
    