# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 09:22:23 2020

@author: Casper Kanaar

This program validates the adapted LL model against experimental tractor propeller configuration data by Dr. Ir. Sinnige. 


"""

# =============================================================================
# Importing relevant modules 
import numpy as np 
from isacalculator import compute_isa
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from numpy.linalg import lstsq
import pandas as pd
import seaborn as sns
sns.set()
# =============================================================================
# Functions 

def compute_n(J,V,D): 
    return V/(J*D)

def compute_T(rho,n,D,CT):
    return rho* n**2 * D**4 * CT

def v_axial_propeller(V_0,T,rho,S_p):
    v_a = 0.5*(-V_0 + np.sqrt(V_0**2 + 2*T/(rho*S_p))) #Equation A.29 PhD Veldhuis
    return v_a 

def v_swirl_propeller(V_0,v_a,omega,r):
    if omega == 0:
        v_swirl = 0
    else:
        v_swirl = (2*V_0*v_a)/(omega*r) #Ferarri 1957
    return v_swirl 

def compute_Re(rho,V,c,mu):
    return rho*V*c/mu

# =============================================================================
# Creating a function of thrust coefficient CT vs. advance ratio J. 
propdata = pd.read_csv("propdata.csv")
CTdata = propdata["CT=T/rho*n2*D4"]
Jdata = propdata["J=Vinf/nD"]
CT_asfunctionof_J = interp1d(Jdata,CTdata)

# =============================================================================
# Parameters 
alpha_fly = 0
y_prop = 0.327
V_infty = 40
mu_infty = 17.89e-6
h = 0
p_infty,rho_infty,T_infty = compute_isa(h)
c = 0.240
b = 0.292
airfoil = "naca64A015"
flap = 0.25
D_p = 0.237
R_p = D_p/2
S_p = R_p**2 *np.pi
N_p = 1
J = [1.0,0.9,0.8,0.7]
n = np.array([compute_n(J[0],V_infty,D_p),compute_n(J[1],V_infty,D_p),compute_n(J[2],V_infty,D_p),compute_n(J[3],V_infty,D_p)])
n = n
CT = [CT_asfunctionof_J(J[0]),CT_asfunctionof_J(J[1]),CT_asfunctionof_J(J[2]),CT_asfunctionof_J(J[3])]
T = [compute_T(rho_infty,n[0],D_p,CT[0]),compute_T(rho_infty,n[1],D_p,CT[1]),compute_T(rho_infty,n[2],D_p,CT[2]),compute_T(rho_infty,n[3],D_p,CT[3])]
alpha_twist = 0 
rotation = ["IU","OU"]

# =============================================================================
# Propeller velocities

# Wing discretization
N = 1000                                                                # [-]. Number of spanwise wing stations. 
K = 100                                                                 # [-]. Number of Fourier modes used. N >= K for a solvable system
 
dy = b/N
y = np.linspace(-b/2+dy/2,b/2-dy/2,N)
theta = np.arccos(-2*y/b)

y_prop = b/2 + (y_prop - b)
y_prop_in = y_prop - R_p
y_prop_out = y_prop + R_p

# Creating the propeller induced velocities.
N_J = len(J)
inducedVelocity_IU = []
inducedVelocity_OU = []
for i in range(N_J): 
    inducedVelocity_i_IU = np.zeros((N,3))
    inducedVelocity_i_OU = np.zeros((N,3))
    Vaxial = v_axial_propeller(V_infty,T[i],rho_infty,S_p)
    
    for j in range(N): 
        if y_prop_in <= y[j] <= y_prop_out:
            inducedVelocity_i_IU[j,0] = Vaxial
            inducedVelocity_i_OU[j,0] = Vaxial
            if y_prop_in <= y[j] <= y_prop:
                inducedVelocity_i_IU[j,2] = v_swirl_propeller(V_infty,Vaxial,n[i],np.abs(y[i] - y_prop))
                inducedVelocity_i_OU[j,2] = -v_swirl_propeller(V_infty,Vaxial,n[i],np.abs(y[i] - y_prop))
            elif y_prop <= y[j] <= y_prop_out: 
                inducedVelocity_i_IU[j,2] = -v_swirl_propeller(V_infty,Vaxial,n[i],np.abs(y[i] - y_prop))
                inducedVelocity_i_OU[j,2] = v_swirl_propeller(V_infty,Vaxial,n[i],np.abs(y[i] - y_prop))
    
    inducedVelocity_IU.append(inducedVelocity_i_IU)
    inducedVelocity_OU.append(inducedVelocity_i_OU)

inducedVelocity = [inducedVelocity_IU,inducedVelocity_OU]

# =============================================================================
# Prandtl's Adapted lifting line model with induced propeller velocities. Computing the circulation and induced lift distribution 

# Finding airfoil's clalpha and alpha 0 lift. 
Re = compute_Re(rho_infty,V_infty,c,mu_infty)

airfoil_data = np.transpose(np.genfromtxt("naca64A015.txt"))
airfoil_alpha = airfoil_data[0]
airfoil_cl = airfoil_data[1]
airfoil_alpha_fit = airfoil_alpha[2:10]
airfoil_cl_fit = airfoil_cl[2:10]
clalpha = (airfoil_cl_fit[-1]-airfoil_cl_fit[0])/(airfoil_alpha_fit[-1]-airfoil_alpha_fit[0])*180/np.pi

alpha_asfunctionof_cl = interp1d(airfoil_cl,airfoil_alpha)
alpha_l0 = np.deg2rad(float(alpha_asfunctionof_cl(0)))

# Linear twist distribution
twist_distribution = []
twist_distribution_function = interp1d([-b/2,0,b/2],[alpha_twist,0,alpha_twist])
for i in range(N): 
    twist_distribution.append(float(twist_distribution_function(y[i])))
twist_distribution = np.array(twist_distribution)

# Calculating geometric angle of attack [rad]
alpha_geo = np.deg2rad(alpha_fly + twist_distribution)

# Linear chord distribution 
chord_distribution = np.full(shape = N,fill_value = c,dtype = np.float)

filenames_cl_distribution = ["clvsbvalitioniu.png","clvsbvalitionou.png"]
IUdata = pd.read_csv("IUdata.csv")
OUdata = pd.read_csv("OUdata.csv")
data = [IUdata,OUdata]

cl = []

# Creating a system of equations to solve for the Fourier coefficients for each of the load cases. 
for i in range(len(rotation)):
    
    
    for j in range(N_J):
        print(" Evaluating "+rotation[i]+" at J = "+str(J[j]))
        
        # Generating the Fourier coefficient matrix to solve the system. 
        M = np.zeros((N,K))
        for k in range(N):
            for l in range(K): 
                M[k,l] = np.sin((l+1)*theta[k])*(np.sin(theta[k]) + (l+1)*(chord_distribution[k]*clalpha/(4*b)))
      
        # Inducing the propeller velocities onto the total velocity vector. 
        V_tot = np.zeros((N,3))
        for k in range(N): 
            V_tot[k,0] = V_infty + inducedVelocity[i][j][k][0]
            V_tot[k,2] = inducedVelocity[i][j][k][2]
        
        # Creating a normalized vector for each N.
        V_norm = np.zeros((N,1))
        for k in range(N):
            V_norm_k = np.sqrt(V_tot[k][0]**2 +V_tot[k][1]**2 + V_tot[k][2]**2)
            V_norm[k][0] = V_norm_k
     
        # Generating the solution vector (b_vector)
        B = np.zeros((N,1))
        for k in range(N):
            B[k,0] = (clalpha*chord_distribution[k]/(4*b))*(V_infty/V_norm[k][0])*((alpha_geo[k] - alpha_l0)+V_tot[k][2]/V_norm[k][0])*np.sin(theta[k])
        
        # Solving the system and calculating the Fourier coefficients. Using a least squares solution to solve the system as N =! K. 
        A = lstsq(M,B,rcond=None)[0]

        # Creating a circulation distribution matrix 
        G = []
        for k in range(N):
            total = 0
            for l in range(K): 
                total = total + A[l]*np.sin(theta[k]*(l+1))
            gamma = total * 2 * b * V_infty
            G.append(gamma)

        # Computing the Clprandtl. This is a vector containing cl values normalized to the local V so not Vinfinity. 
        clprandtl = np.zeros((N,1))
        for k in range(N): 
            clprandtl[k,0] = 2*G[k] / (V_infty *chord_distribution[k])
        
        cl.append(clprandtl)
        
# =============================================================================
# Plotting the cl distributions against the validation data. 
colors = ["red","green","orange","purple"]
markers = ["o","v","^","D","p"]

for i in range(len(rotation)): 
    fig = plt.figure(dpi = 250)
    
    cl_data = data[i]["cl"]
    polar_data = data[i]["polar"]
    y_data = data[i]["Y/s"]
    
    if i == 0:
        
        for j in range(0,4):
            cl_j  = []
            y_j = []
            for k in range(len(cl_data)): 
                if polar_data[k] - 2 == j: 
                    cl_j.append(cl_data[k])
                    y_j.append(y_data[k])
            cl_j = np.array(cl_j)
            y_j = np.array(y_j)
            y_j =  y_j*b - b/2 
            y_j = np.fliplr([y_j])[0]
            cl_j = np.fliplr([cl_j])[0]
            plt.scatter(y_j,cl_j,color = colors[j],marker = markers[j],label = "Experimental data. J = "+str(J[j]))
            plt.plot(y,cl[j],color=colors[j],label = "Adapted LL model. J = "+str(J[j]))
        
        plt.xlabel("b [m]")
        plt.ylabel("$C_l$ [-]")
        plt.legend()
        plt.savefig(filenames_cl_distribution[i])
    
    elif i == 1:
        fig = plt.figure(dpi=250)
        for j in range(4,8):
            cl_j  = []
            y_j = []
            for k in range(len(cl_data)): 
                if polar_data[k] - 2 == j-4: 
                    cl_j.append(cl_data[k])
                    y_j.append(y_data[k])
            cl_j = np.array(cl_j)
            y_j = np.array(y_j)
            y_j =  y_j*b - b/2 
            y_j = np.fliplr([y_j])[0]
            cl_j = np.fliplr([cl_j])[0]
            plt.scatter(y_j,cl_j,color = colors[j-4],marker = markers[j-4],label = "Experimental data. J = "+str(J[j-4]))
            plt.plot(y,cl[j],color=colors[j-4],label = "Adapted LL model. J = "+str(J[j-4]))
        
        plt.xlabel("b [m]")
        plt.ylabel("$C_l$ [-]")
        plt.legend()
        fig.savefig(filenames_cl_distribution[i])
        
# =============================================================================
# End 
print("Propeller-wing interaction validation finished")        
                
     
        
    
    

 

    
        
        






    

        
        
    












