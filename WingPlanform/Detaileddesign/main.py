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
import scipy.integrate as sc


# Import data

xflr = np.loadtxt("T1-28_0 m_s-VLM2.txt",skiprows=8)
CL_xflr = xflr[:,2]
CD_xflr = xflr[:,5]
CD0_wing_xflr = dc.compute_linear(CL_xflr[20],CL_xflr[21],CD_xflr[20],CD_xflr[21],0)


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
taper_t = 0.2                                            # [-]
cr_t    = 0.15                                           # [m]
ct_t    = cr_t * taper_t                                 # [m]
MAC_t   = 2/3*cr_t*(1+taper_t+taper_t**2)/(1+taper_t)    # [m]
b       = 3                                              # [m]
S       = (cr + ct) * b/2                                # [m^2]
AR      = b**2/S                                         # [-]
b_t     = 0.882                                          # [m]
S_t     = (cr_t + ct_t) * b_t/2                          # [m^2]
AR_t    = b_t**2/S_t                                     # [-]
V_cr    = 28                                             # [m/s]
V_tr    = 14                                             # [m/s]
LE_sw   = 23.1                                           # [deg]
LE_sw_t = 39.1                                           # [deg]
clmax   = 1.37                                           # [-]
QC_sw   = dc.compute_sweep(LE_sw,taper,0.25,cr,b)        # [deg]
QC_sw_t = dc.compute_sweep(LE_sw_t,taper_t,0.25,cr_t,b_t)# [deg]
cldes   = dc.compute_cldes(MTOM,rho,V_cr,QC_sw,S)        # [-]
clalpha = 0.11067579                                     # [deg^-1]
cl0     = 0.0074437                                      # [-]
Re_m    = rho*V_cr/0.0000179579                          # [m^-1]
M_cr    = 0.084                                          # [-]

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
        aoa_des = (cldes - cl0)/clalpha                              # [deg]
        aoa_zero = dc.compute_AOA_0_lift(clalpha,cldes,aoa_des)      # [deg]
        theta = -4                                                   # [deg]
        aoa0_theta = -0.39                                           # [-]
        aoa_0 = dc.compute_AOA_0_lift_tw(aoa0_theta,theta,aoa_zero)  # [deg]
    
    # Maximum lift coefficient parameters    
    delta_y = 3.72     # Airfoil sharpness parameter                 # [-]
    CL_cl = 0.83                                                     # [-]
    delta_CLmax = 0                                                  # [-]
    CLmax = dc.compute_CLmax_high(CL_cl,clmax,delta_CLmax)           # [-]
    
    # Angle of attack at maximum lift
    delta_aoa = 2.4                                                  # [deg]
    aoa_stall = dc.compute_aoa_stall_high(CLmax,CLalpha,aoa_0,delta_aoa) # [deg]
    
    # Lift coefficients below stall
    CN_prime_max = dc.compute_CN_prime_CLmax(CLmax,aoa_stall)        # [-]
    J = dc.compute_Jpar(C1,C2,LE_sw,AR)                              # [-]
    aoa = np.linspace(-5,aoa_stall,10000)                            # [deg]
    # aoa_e = dc.compute_aoa_camber(aoa, aoa_0, aoa_stall)
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
    CL = dc.compute_CL(CN_prime,aoa)                           # [-]
    
    
    # Wing parameters to calculate zero-lift drag coefficient of the wing
    tc_avg = 0.1                  # Thickness of CAL4014L            # [-]
    x_tcmax = 0.2494                                                 # [-]
    Sref = S                                                         # [m^2]
    Swet = 2*S                                                       # [m^2]
    Cf = dc.compute_Cf(M_cr,Re_m,taper,cr)                           # [-]
    TMAX_sw = dc.compute_sweep(LE_sw,taper,x_tcmax,cr,b)             # [-]
    Rls = dc.compute_Rls(M_cr,QC_sw)                                 # [-]
    # CD0_wing = dc.compute_CD0_wing(Cf,tc_avg,x_tcmax,Swet,Sref,1)    # [-]
    CD0_wing = CD0_wing_xflr

    # Tail parameters to calculate zero-lift drag coefficient of the tail
    NACA0009 = np.loadtxt("NACA0009.txt")                            # [-]
    tc_avg_t = 0.09                                                  # [-]
    x_tcmax_t = 0.2903                                               # [-]
    Swet_t = 2*S_t                                                   # [m^2]
    Cf_t = dc.compute_Cf(M_cr,Re_m,taper_t,cr_t)                     # [-]
    TMAX_sw_t = dc.compute_sweep(LE_sw_t,taper_t,x_tcmax_t,cr_t,b_t) # [-]
    Rls_t = dc.compute_Rls(M_cr,QC_sw_t)                             # [-]
    CD0_tail = dc.compute_CD0_tail(Cf_t,tc_avg_t,x_tcmax_t,Swet_t,Sref,Rls_t)  # [-]
    
    
    # Body parameters to calculate zero-lift drag coefficient of the body
    lb = 0.743                                                       # [m]
    d = 0.379    # Max diameter                                      # [m]
    h = 0.167    # Max height                                        # [m]
    Cf_b = dc.compute_Cf_b(M_cr,Re_m,lb)                             # [-]
    Sb = 2*((0.488-0.135)*d + (0.743-0.488)*d/2 + np.pi*0.125*d**2) # [m^2]    
    Sb_Sref = Sb/Sref                                                # [-]
    CD0_body = dc.compute_CD0_body(Cf_b,lb,Sb_Sref,Sb,h)             # [-]
    
    # Total zero-lift drag coefficient
    CD0 = CD0_wing + CD0_tail + CD0_body                             # [-]

    # Lift induced drag and Oswald efficiency factor
    R = 0.94
    if twist == False:
        CDi_wing,e = dc.compute_CD_ind_wing(CLalpha,AR,clalpha,0,0,0,R,CL) # [-]
            
    else:
        CDi_wing,e = dc.compute_CD_ind_wing(CLalpha,AR,clalpha,0.000625,0.0019,theta,R,CL) # [-]
    
    # Total drag
    CD = dc.compute_CD(CD0,CDi_wing)                                 # [-]
    
    # Updated cruise speed
    V_cr_update = np.sqrt(MTOW*2/(rho*S*CL[np.argmax(CL/CD)])) # [m/s]
    
    
    
print("CD0 =", CD0)
print("Max. L/D =",np.max(CL/CD))
print("Cruise AoA =", aoa[np.argmax(CL/CD)])
print("Oswald efficiency factor e =", e)
print("Max CL =", CLmax)
print("Lift-slope CLalpha =",CLalpha)
print("Updated cruise speed =",V_cr_update)
print("Cruise CL =", CL[np.argmax(CL/CD)])
print("Climb CD =", CD[np.argmax(CL)])

V = np.sqrt(MTOW*2/(rho*S*CL[CL>0]))

# --------------------------- Extra Computations --------------------------- #

data = np.loadtxt("CAL4014L.dat")

# Cf = dc.compute_Cf(0.084,rho*V_cr/0.0000179579,taper,cr)
# Cf_b = dc.compute_Cf_b(0.084,rho*V_cr/0.0000179579,lb)
# Rls = dc.compute_Rls(0.084, QC_sw)
# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiFvvr06-rpAhWPzKQKHSQyAJQQFjABegQIAhAB&url=https%3A%2F%2Fwww.mdpi.com%2F2076-3417%2F9%2F15%2F3043%2Fpdf&usg=AOvVaw1BMA1rrBr4OXRl5WVuMs9s
# dc.compute_CD0_wing(test_cf,0.098,0.276,2*1.07,1.07,1.06)