# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 09:22:23 2020

@author: Casper

This program models the propeller wing interaction in order to compute a more accurate representation of the aircraft's lifting line. 


"""
# =============================================================================
# Importing relevant modules 
import numpy as np 
from isacalculator import compute_isa
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

# =============================================================================
# Parameters 
V_cruise = 28                                               # [m/s]. Cruise free stream velocity          
V_stall = 14                                                # [m/s]. Stall free stream velocity
h_cruise = 500                                              # [m]. Cruise altitude 
h_stall = 20                                                # [m]. Transition altitude 
p_cruise,rho_cruise,T_cruise = compute_isa(h_cruise)        # [Pa],[kg/m^3],[K]. Cruise air pressure, density, and temperature
p_stall, rho_stall, T_stall = compute_isa(h_stall)          # [Pa],[kg/m^3],[K]. Stall air pressure, density, and temperature
M_cruise = V_cruise/np.sqrt(1.4*287*T_cruise)               # [-]. Cruise Mach number
b = 3                                                       # [m]. Wing span
r0 = 1                                                      # [m^2/s]. Starting value of circulation at the origin.    
alpha_l0 = -0.06 /180*np.pi                                 # [rad]. Airfoil zero-lift angle of attack. This will become a function of theta if wing twist is introduced.
W_S = 122.23                                                # [N/m^2]. Wing loading 
m = 17.53                                                   # [kg]. MTOM
S = m*9.81/W_S                                              # [m^2]. Wing surface area 

# TEMPORARY DUMMY VALUES
w_p = 10                                                    # [m/s]. Propeller induced downwash velocity
u_p = 10                                                    # [m/s]. Propeller induced axial velocity          

# =============================================================================
# Functions 
def c(theta):
    span = 3
    y = -span/2*np.cos(theta)
    c_root = 0.696
    c_tip = 0.244
    dydx = (c_root - c_tip)/(span/2)
    b = c_root 
    if y >= 0:
        c = b - y*dydx
    else: 
        c = b + y*dydx
    
    return c
    
# =============================================================================
# MAIN. Creating a propeller adapted LL method. 

N = 200                         # [-]. Number of spanwise stations. Change to increase accuracy 
alpha_geo = 0                   # [rad]. Geometric angle of attack for analysis
theta = np.linspace(0,np.pi,N)  # [rad]. Theta locations for spanwise evaluation.

# Creating a system of equations according to the fundamental equation of Prandtl's LL theory expressed in Fourrier coefficients. 
system_matrix = np.zeros((N,N))

for i in range(len(theta)):
    for n in range(1,N+1):
        if n*theta[i] == theta[i]:
            k = 2*b / (np.pi*c(theta[i])) * np.sin(n*theta[i]) + n*1
        else:
            k = 2*b / (np.pi*c(theta[i])) * np.sin(n*theta[i]) + n*np.sin(n*theta[i])/np.sin(theta[i])
        system_matrix[i,n-1] = k
        
        
# Creating the result vector for the system matrix to be equated to. 
result = np.zeros((N,1))

y =  alpha_geo - alpha_l0 - w_p/(V_cruise + u_p)
for i in range(len(result)):
    result[i,0] = y 
    
A_coefs = np.linalg.solve(system_matrix,result)

# Obtaining a circulation distribution from the Fourier coefficients. 
R = []
for i in range(len(theta)): 
    total = 0 
    for j in range(len(A_coefs)):
        total = total + A_coefs[j]*np.sin(j*theta[i])
    R_i = 2*b*V_cruise*total 
    R.append(R_i)
    
    
y = -b/2*np.cos(theta)

plt.plot(y,R)

    
        
        






    

        
        
    












