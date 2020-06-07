# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 09:22:23 2020

@author: Casper

This program models the propeller wing interaction in order to compute a more accurate representation of the aircraft's lifting line. 
The propeller-wing interaction is modelled according to propeller actuator disk theory. The lifting line is model according to an adapted version of 
Prandtl's lifting line by E.K. Epema. 

Verification of Actuator disk model: - ? 
Verification of the LL model: XFLR5


"""
# =============================================================================
# SWITCHES 
""" 
Switch the propellers on or off:
    PROPELLERS = True : Propellers are turned on 
    PROPELLERS = False: Propellers are turned off
"""

PROPELLERS = True   

# =============================================================================
# Importing relevant modules 
import numpy as np 
from isacalculator import compute_isa
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from numpy.linalg import lstsq
import seaborn as sns
import statsmodels.api as sm 
sns.set()

# =============================================================================
# Parameters 
V_cruise = 28                                               # [m/s]. Cruise free stream velocity          
V_stall = 14                                                # [m/s]. Stall free stream velocity
V_VTOL = 0                                                  # [m/s]. Start of VTOL free stream velocity 
h_cruise = 500                                              # [m]. Cruise altitude 
h_stall = 20                                                # [m]. Transition altitude 
h_VTOL = 0                                                  # [m]. Start of VTOL altitude 
p_cruise,rho_cruise,Temp_cruise = compute_isa(h_cruise)     # [Pa],[kg/m^3],[K]. Cruise air pressure, density, and temperature
p_stall, rho_stall, T_stall = compute_isa(h_stall)          # [Pa],[kg/m^3],[K]. Stall air pressure, density, and temperature
p_VTOL,rho_VTOL,T_VTOL = compute_isa(h_VTOL)                # [Pa],[kg/m^3],[K]. VTOl air pressure, density, and temperature 
M_cruise = V_cruise/np.sqrt(1.4*287*Temp_cruise)            # [-]. Cruise Mach number
b = 3                                                       # [m]. Wing span
W_S = 122.23                                                # [N/m^2]. Wing loading 
m = 17.53                                                   # [kg]. MTOM
S = m*9.81/W_S                                              # [m^2]. Wing surface area 
twist = 0                                                   # [deg]. Wing tip twist angle (linear distribution)
AR = b**2/S                                                 # [-]. Wing aspect ratio 
taper = 0.35                                                # [-]. Wing taper ratio 
t_c = 0.10                                                  # [-]. Wing airfoil maximum thickness over chord ratio.
c_root = 2*S/(b*(1+taper))                                  # [m]. Wing root chord                
c_tip = c_root*taper                                        # [m]. Wing tip chord

# =============================================================================
# Propeller actuator disk model.
D_p = 0.3                                                   # [m]. Propeller disk diameter 
R_p = D_p/2                                                 # [m]. Propeller disk radius 
S_p = np.pi*R_p**2                                          # [m^2]. Propeller disk area 
T_cruise = 20                                               # [N]. Propeller thrust during cruise 
T_VTOL = 50                                                 # [N]. Propeller thrust during VTOL 
N_p = 4                                                     # [-]. Number of propellers                                                  
omega_cruise = 3                                            # [rad/s]. Propeller angular velocity during cruise            
omega_VTOL = 3                                              # [rad/s]. Propeller angular velocity during VTOL 

# Wing discretization
N = 1000                                                    # [-]. Number of spanwise wing stations. 
K = 100                                                     # [-]. Number of Fourier modes used. N>K for a solvable system
alpha_fly = 2                                               # [deg]. Geometric angle of attack. 
dy = b/N                                                    # [m]. Width of panels 
y = np.linspace(-b/2+dy/2,b/2-dy/2,N)                       # [-]. Control point locations on mid-panel 
theta = np.arccos(-2*y/b)                                   # [-]. Coordinate transformation 

# Defining and computing the axial velocity induced by the propeller. 
def v_axial_propeller(V_0,T,rho,S_p):
    v_a = 0.5*(-V_0 + np.sqrt(V_0**2 + 2*T/(rho*S_p))) #Equation A.29 PhD Veldhuis
    return v_a 

# Defining and computing the radial/swirl velocity induced by the propellers
def v_swirl_propeller(V_0,v_a,omega,R_p):
    v_swirl = (2*V_0*v_a)/(omega*R_p) #Ferarri 1957
    return v_swirl 

V_a_cruise = v_axial_propeller(V_cruise,T_cruise,rho_cruise,S_p)
V_swirl_cruise = v_swirl_propeller(V_cruise,V_a_cruise,omega_cruise,R_p)

# Creating a matrix of induced propeller velocities according to the propeller placements along the span. 
inducedVelocity = np.zeros((N,3))

# Propeller positions. Engines are placed symmetrically (laterally) at (b/2) 0.35 and (b/2)*0.7
y_inner = 0.35
y_outer = 0.7

y_lim_inner_Pinner = b/2 * y_inner - R_p
y_lim_outer_Pinner = b/2 * y_inner + R_p 
y_lim_inner_Pouter = b/2 * y_outer - R_p 
y_lim_outer_Pouter = b/2 * y_outer + R_p 

for i in range(N): 
    if -y_lim_outer_Pouter <= y[i] <= - y_lim_inner_Pouter:
        inducedVelocity[i,0] = V_a_cruise
        inducedVelocity[i,2] = V_swirl_cruise
    elif -y_lim_outer_Pinner <= y[i] <= -y_lim_inner_Pinner:
        inducedVelocity[i,0] = V_a_cruise
        inducedVelocity[i,2] = V_swirl_cruise
    elif y_lim_inner_Pinner <= y[i] <= y_lim_outer_Pinner:
        inducedVelocity[i,0] = V_a_cruise
        inducedVelocity[i,2] = V_swirl_cruise
    elif y_lim_inner_Pouter <= y[i] <= y_lim_outer_Pouter:
        inducedVelocity[i,0] = V_a_cruise
        inducedVelocity[i,2] = V_swirl_cruise
        
if PROPELLERS == True: 
    inducedVelocity = inducedVelocity
else: 
    inducedVelocity = np.zeros((N,3))

# =============================================================================
# Prandtl's Adapted lifting line model with induced propeller velocities. Computing the circulation and induced lift distribution 

# Finding airfoil's clalpha and alpha 0 lift. 
reynolds_data_transposed = np.genfromtxt("reynoldsdistribution.txt")
reynolds_data = np.transpose(reynolds_data_transposed)
Re_min = min(reynolds_data[1])
Re_max = max(reynolds_data[1])
airfoil_data_transposed = np.genfromtxt("cal4014lcla.txt")
airfoil_data = np.transpose(airfoil_data_transposed)
airfoil_alpha = airfoil_data[0][0:25]
airfoil_cl = airfoil_data[1][0:25]
plt.plot(airfoil_alpha,airfoil_cl)
x = sm.add_constant(airfoil_alpha)
results = sm.OLS(airfoil_cl,x).fit()
clalpha = np.rad2deg(0.1154)

alpha_asfunctionof_cl = interp1d(airfoil_cl,airfoil_alpha)
alpha_l0 = np.deg2rad(float(alpha_asfunctionof_cl(0)))

# Linear twist distribution
twist_distribution_function = interp1d([-b/2,0,b/2],[twist,0,twist])
twist_distribution = []
for i in range(N): 
    twist_distribution.append(twist_distribution_function(y[i]))
twist_distribution = np.array(twist_distribution)
    
# Calculating geometric angle of attack [rad]
alpha_geo = (alpha_fly + twist_distribution)*np.pi/180

# Linear chord distribution 
chord_distribution_function = interp1d([-b/2,0,b/2],[c_tip,c_root,c_tip])
chord_distribution = []
for i in range(N): 
    chord_distribution.append(chord_distribution_function(y[i]))
chord_distribution = np.array(chord_distribution)

# Generating the Fourier coefficient matrix to solve the system. 
M = np.zeros((N,K))
for i in range(N):
    for j in range(K): 
        M[i,j] = np.sin((j+1)*theta[i])*(np.sin(theta[i]) + (j+1)*(chord_distribution[i]*clalpha/(4*b)))
        
# Inducing the propeller velocities onto the total velocity vector. 
V_tot = np.zeros((N,3))
for i in range(N): 
    V_tot[i,0] = V_cruise + inducedVelocity[i,0]
    V_tot[i,2] = inducedVelocity[i,2] 
    
# Creating a normalized vector for each N.
V_norm = np.zeros((N,1))
for i in range(N):
    V_norm_i = np.sqrt(V_tot[i][0]**2 +V_tot[i][1]**2 + V_tot[i][2]**2)
    V_norm[i][0] = V_norm_i

# Generating the solution vector (b_vector)
B = np.zeros((N,1))
for i in range(N):
    B[i,0] = (clalpha*chord_distribution[i]/(4*b))*(V_cruise/V_norm[i][0])*((alpha_geo[i] - alpha_l0 - twist_distribution[i])+V_tot[i][2]/V_norm[i][0])*np.sin(theta[i])

# Solving the system and calculating the Fourier coefficients
A = lstsq(M,B,rcond=None)[0]

# Creating a circulation distribution matrix 
G = []
for i in range(N):
    total = 0
    for j in range(K): 
        total = total + A[j]*np.sin(theta[i]*(j+1))
    total = total * 2 * b * V_cruise 
    G.append(total)

# Computing the induced angle of attack distribution from the circulation matrix
alpha_induced = []
for i in range(N):
    total = 0 
    for j in range(K): 
        total = total + (j+1)*A[j]*(np.sin((j+1)*theta[i])/np.sin(theta[i]))
    alpha_induced.append(total)
    
# =============================================================================
# Post processing. Computing all aerodynamic coefficients based off the adapted LL model. 
 
# 3D lift coefficient at this angle of attack. 
C_L = float(A[0][0]*np.pi*AR)
        
# Lift induced drag coefficient at this angle of attack. 
delta = 0
for i in range(1,K):
    delta = delta + (i+1)*(A[i][0]/A[0][0])**2
C_Di = float(np.pi*AR*A[0][0]**2 *(1+delta))

# Calculating the span efficiency factor. 
e1 = C_L**2/(np.pi*C_Di*AR)
e2 = 1/(1+delta)
# Assert whether AR checks out
AR_check = C_L**2/(np.pi*C_Di*e2)

"""VERIFICIATION STEP"""
assert np.isclose(AR,float(AR_check),rtol=0.01), "Verification failed: Lifting line method not correctly implemented"
    
# =============================================================================
# Propellers induction

# =============================================================================
# # Calculating the downwash at the wing [m/s]. 
# downwash = np.zeros((N,1))
# for i in range(N):
#     downwash[i,0] = -V_norm[i][0]*alpha_induced[i][0] 
# 
# for i in range(N):
#     V_tot[i][2] = downwash[i][0]
# =============================================================================
    
# Computing the Clprandtl. This is a vector containing cl values normalized to the local V so not Vinfinity. 
clprandtl = np.zeros((N,1))
for i in range(N): 
    clprandtl[i,0] = 2*G[i][0] / (V_cruise * chord_distribution[i])  

# =============================================================================
# Generating the lift distribution 

# Kutta-Jouwkowski theorem: Applying it with the freestream velocity instead of propeller induced velocity according to Dr. Sinnige. 
L = np.zeros((N,1))
for i in range(N):
    L[i,0] = rho_cruise*V_cruise*G[i][0]
    
c_l = np.zeros((N,1))
for i in range(N): 
    c_l[i,0] = 2*G[i][0]/(V_cruise*chord_distribution[i])
    
# =============================================================================
# Importing the XFLR5 LL for verification purposes. 
data_xflr_transposed = np.genfromtxt("LLxflr5SAVEDunswept.txt")
data_xflr = np.transpose(data_xflr_transposed)
y_xflr = data_xflr[0]
cl_xflr = data_xflr[1]
    
# Plotting and formatting the data 
fig = plt.figure(figsize = (10,5),dpi=250)
plt.plot(y,L, label = "Adapted LL model")
plt.scatter(y_xflr,cl_xflr, label = "XFLR5 LL model",color = "indianred")
plt.xlabel("b [m]")
plt.ylabel("$C_L$ [-]")
plt.legend()
plt.savefig("liftdistribution.png")

lift_file = open("liftdistribution.txt","w")
lift_file.close()

lift_file = open("liftdistribution.txt","w")
lines = []
for i in range(N): 
    lines.append(str(L[i][0])+"\n")
    print(i,str(L[i][0]))
lift_file.writelines(lines)

# =============================================================================













    
        
        






    

        
        
    












