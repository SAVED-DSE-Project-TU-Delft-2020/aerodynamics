# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 09:22:23 2020

@author: Casper Kanaar

This program models the propeller wing interaction in order to compute a more accurate representation of the aircraft's lifting line. 
The propeller-wing interaction is modelled according to propeller actuator disk theory. The lifting line is model according to an adapted version of 
Prandtl's lifting line by E.K. Epema and R. Nederhof. 

The model features 3 load cases.
    - SAVED in cruise with propellers switched on (LC1).
    - SAVED in VTOL (hovering mode) with propellers switched on (LC2).
    - Arbitrary, low AR wing for verification with XFLR5's LL model (LC3). 
    
[L1,LC2,LC3]

Verification of Actuator disk model: - ? 
Verification of the LL model: XFLR5

"""
# =============================================================================
# SWITCHES AND INPUTS  
""" 
Switch the propellers on or off:
    PROPELLERS = True : Propellers are turned on 
    PROPELLERS = False: Propellers are turned off
"""

PROPELLERS = True   

""" 
Angle of attack [deg] the aircraft is flying at for each of the load cases. 
V_cruise [m/s] is the cruise speed of each of the cases.  
""" 
alpha_fly = [5.32,0,2]
V_cruise = [22.5,0.0001,28]                                                                                                   # [m/s]. Stall free stream velocity


"""
Switch the results of XFLR5 on or off in the plot for verification purposes.
    XFLR5 = True: Scatter is turned on 
    XFLR5 = False: Scatter is turned off

Note: 
    - The verification lifting line will always have XFLR5's scatter turned on. 
    - XFLR5 does not include the propeller-wing interaction. 
""" 
XFLR5 = False  

# =============================================================================
# Importing relevant modules 
import numpy as np 
from isacalculator import compute_isa
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from numpy.linalg import lstsq
import seaborn as sns
sns.set()

# =============================================================================
# Parameters [Load Case 1, Load Case 2, Load Case 3]
N_LC = 3                                                                 # [-]. Number of load cases 
h_cruise = [500,0,500]                                                   # [m]. Cruise altitude 
p_cruise = []                                                            # [Pa]. Cruise air pressure 
rho_cruise = []                                                          # [kg/m^3]. Cruise air density 
Temp_cruise = []                                                         # [K]. Cruise air temperature 

for i in range(N_LC):
    p,rho,Temp = compute_isa(h_cruise[i])
    p_cruise.append(p)
    rho_cruise.append(rho)
    Temp_cruise.append(Temp)

M_cruise = []                                                            # [-]. Cruise Mach number

for i in range(N_LC): 
    M_i = V_cruise[i]/np.sqrt(1.4*287*Temp_cruise[i])
    M_cruise.append(M_i)
            
sweep = np.deg2rad(np.array([0,0,0]))                                   # [rad]. Main wing sweep
b = [3,3,2.5]                                                            # [m]. Wing span
W_S = 130.034                                                              # [N/m^2]. Wing loading 
m = 17.4754                                                            # [kg]. MTOM
S = [m*9.81/W_S,m*9.81/W_S, 1.25]                                        # [m^2]. Wing surface area 
twist = [0,0,0]                                                          # [deg]. Wing tip twist angle (linear distribution)

AR = []
for i in range(N_LC):
    AR.append(b[i]**2/S[i])

taper = [0.4,0.4,1]                                                    # [-]. Wing taper ratio 
t_c = [0.10,0.10,0.12]                                                   # [-]. Wing airfoil maximum thickness over chord ratio.

c_root = []                                                              # [m]. Wing root chord
for i in range(N_LC-1):
    c_root.append(2*S[i]/(b[i]*(1+taper[i])))
c_root.append(0.5)

c_tip = []                                                               # [m]. Wing tip chord
for i in range(N_LC): 
    c_tip.append(c_root[i]*taper[i])
    
CD = 0.019728716372218168                                                # [m]. Cruise drag coefficient 
# =============================================================================
# Propeller actuator disk model.
D_p = 15.5*0.0254                                                       # [m]. Propeller disk diameter 
R_p = D_p/2                                                             # [m]. Propeller disk radius 
S_p = np.pi*R_p**2                                                      # [m^2]. Propeller disk area 
N_p = 4                                                                 # [-]. Number of propellers                                                  

T_cruise = [0.5*rho_cruise[0]*V_cruise[0]**2 * S[0]*CD,1.2*m*9.81,0]    # [N]. Propeller thrust 

# Computing the propeller angular velocity based on a linear regression between thrust and RPM.
def omega(T):
    omega_rpm = 106*T + 2831
    omega_rads = 0.104719755*(omega_rpm)
    
    return omega_rads

omega_cruise = [omega(T_cruise[0]),omega(T_cruise[1]),omega(T_cruise[2])]                                                  # [rad/s]. Propeller angular velocity during cruis

# Wing discretization
N = 1000                                                                # [-]. Number of spanwise wing stations. 
K = 10                                                                 # [-]. Number of Fourier modes used. N >= K for a solvable system

dy = []
for i in range(N_LC):
    dy.append(b[i]/N)
    
y = []                                                                  # [-]. Control point locations on mid-panel 
for i in range(N_LC):
    y_i = np.linspace(-b[i]/2+dy[i]/2,b[i]/2-dy[i]/2,N)
    y.append(y_i)

theta = []                                                              # [-]. Coordinate transformation 
for i in range(N_LC):
    theta_i = np.arccos(-2*y[i]/b[i])
    theta.append(theta_i)

# Propeller positions. Engines are placed symmetrically (laterally) at (b/2) 0.35 and (b/2)*0.7
y_inner = 0.35
y_outer = 0.7

# Defining and computing the axial velocity induced by the propeller. 
def v_axial_propeller(V_0,T,rho,S_p):
    v_a = 0.5*(-V_0 + np.sqrt(V_0**2 + 2*T/(rho*S_p))) #Equation A.29 PhD Veldhuis
    return v_a 

# Defining and computing the radial/swirl velocity induced by the propellers
def v_swirl_propeller(V_0,v_a,omega,r):
    if omega == 0:
        v_swirl = 0
    else:
        v_swirl = (2*V_0*v_a)/(omega*r) #Ferarri 1957
    return v_swirl 

V_a_cruise = []
for i in range(N_LC):
    V_a = v_axial_propeller(V_cruise[i],T_cruise[i],rho_cruise[i],S_p)
    V_a_cruise.append(V_a)
    
    
def n_prop(V_a,V_cruise): 
    return  V_cruise/(V_cruise+V_a)

n_p = [n_prop(V_a_cruise[0],V_cruise[0]),n_prop(V_a_cruise[1],V_cruise[1]),n_prop(V_a_cruise[2],V_cruise[2])]
    

# Creating a matrix of induced propeller velocities according to the propeller placements along the span. 
inducedVelocity = []
for i in range(N_LC): 
    inducedVelocity_i = np.zeros((N,3))
    
    y_lim_inner_Pinner = b[i]/2 * y_inner - R_p
    y_lim_outer_Pinner = b[i]/2 * y_inner + R_p 
    y_lim_inner_Pouter = b[i]/2 * y_outer - R_p 
    y_lim_outer_Pouter = b[i]/2 * y_outer + R_p 

    for j in range(N):
        # Left outboard propeller: Clockwise from behind. 
        if -y_lim_outer_Pouter <= y[i][j] <= - y_lim_inner_Pouter:
            inducedVelocity_i[j,0] = V_a_cruise[i]
            if -y_lim_outer_Pouter <= y[i][j] <= -b[i]/2 * y_outer:
                inducedVelocity_i[j,2] = v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_outer - np.abs(y[i][j])))
            elif -b[i]/2*y_outer <= y[i][j] <= - y_lim_inner_Pouter:
                inducedVelocity_i[j,2] = -v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_outer - np.abs(y[i][j])))
        # Left inboard propeller: Counter clockwise from behind
        elif -y_lim_outer_Pinner <= y[i][j] <= -y_lim_inner_Pinner:
            inducedVelocity_i[j,0] = V_a_cruise[i]
            if -y_lim_outer_Pinner <= y[i][j] <= -b[i]/2 * y_inner:
                inducedVelocity_i[j,2] = -v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_inner - np.abs(y[i][j])))
            elif -b[i]/2 * y_inner <= y[i][j] <= -y_lim_inner_Pinner:
                inducedVelocity_i[j,2] = v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_inner - np.abs(y[i][j])))
                
        # Right inboard propeller: Counter clockwise from behind 
        elif y_lim_inner_Pinner <=  y[i][j] <= y_lim_outer_Pinner:
            inducedVelocity_i[j,0] = V_a_cruise[i]
            if y_lim_inner_Pinner <=  y[i][j] <= b[i]/2 * y_inner:
                inducedVelocity_i[j,2] = v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_inner - np.abs(y[i][j])))
            elif b[i]/2 * y_inner <= y[i][j] <= y_lim_outer_Pinner:
                inducedVelocity_i[j,2] = -v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_inner - np.abs(y[i][j])))
                
        # Right outboard propeller: Clockwise from behind
        elif y_lim_inner_Pouter <=  y[i][j] <= y_lim_outer_Pouter:
            inducedVelocity_i[j,0] = V_a_cruise[i]
            if y_lim_inner_Pouter <=  y[i][j] <= b[i]/2 * y_outer:
                inducedVelocity_i[j,2] = -v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_outer - np.abs(y[i][j])))
            elif b[i]/2 * y_outer <= y[i][j] <= y_lim_outer_Pouter:
                inducedVelocity_i[j,2] = v_swirl_propeller(V_cruise[i],V_a_cruise[i],omega_cruise[i],np.abs(np.abs(b[i]/2) * y_outer - np.abs(y[i][j])))
                
    if PROPELLERS == False:
        inducedVelocity_i = np.zeros((N,3))
    else:
        inducedVelocity_i = inducedVelocity_i
    
    inducedVelocity.append(inducedVelocity_i)
        
# =============================================================================
# Prandtl's Adapted lifting line model with induced propeller velocities. Computing the circulation and induced lift distribution 

# Finding airfoil's clalpha and alpha 0 lift. 
reynolds_data_transposed = np.genfromtxt("reynoldsdistribution.txt")
reynolds_data = np.transpose(reynolds_data_transposed)
Re_min = min(reynolds_data[1])
Re_max = max(reynolds_data[1])

airfoils = ["cal4014l","cal4014l","naca0012"]
clalpha = []
alpha_l0 = []

for i in range(N_LC): 
    airfoil_data_transposed = np.genfromtxt(airfoils[i]+"cla.txt")
    airfoil_data = np.transpose(airfoil_data_transposed)
    airfoil_alpha = airfoil_data[0]
    airfoil_cl = airfoil_data[1]
    airfoil_alpha_fit = airfoil_alpha[2:10]
    airfoil_cl_fit = airfoil_cl[2:10]
    clalpha_i = (airfoil_cl_fit[-1]-airfoil_cl_fit[0])/(airfoil_alpha_fit[-1]-airfoil_alpha_fit[0])*180/np.pi
    clalpha.append(clalpha_i)

    alpha_asfunctionof_cl = interp1d(airfoil_cl,airfoil_alpha)
    alpha_l0_i = np.deg2rad(float(alpha_asfunctionof_cl(0)))
    alpha_l0.append(alpha_l0_i)
    
# Linear twist distribution
twist_distribution = []
for i in range(N_LC):
    twist_distribution_i = []
    twist_distribution_function = interp1d([-b[i]/2,0,b[i]/2],[twist[i],0,twist[i]])
    for j in range(N): 
        twist_distribution_i.append(float(twist_distribution_function(y[i][j])))
    twist_distribution_i = np.array(twist_distribution_i)
    twist_distribution.append(twist_distribution_i)
    
# Calculating geometric angle of attack [rad]
alpha_geo = [] 
for i in range(N_LC): 
    alpha_geo_i = np.deg2rad(alpha_fly[i] + twist_distribution[i])
    alpha_geo.append(alpha_geo_i)
    
# Linear chord distribution 
chord_distribution = []
for i in range(N_LC):
    chord_distribution_i = []
    chord_distribution_function = interp1d([-b[i]/2,0,b[i]/2],[c_tip[i],c_root[i],c_tip[i]])
    for j in range(N): 
        chord_distribution_i.append(chord_distribution_function(y[i][j]))
    chord_distribution_i = np.array(chord_distribution_i)
    chord_distribution.append(chord_distribution_i)

filenames_verification_cl = ["clxflr5lc1.txt","clxflr5lc2.txt","clxflr5lc3.txt"]
filenames_cl_distribution = ["clvsblc1.png","clvsblc2.png","clvsblc3.png"]
filenames_l_distribution = ["lvsblc1.txt","lvsblc2.txt","lvsblc3.txt"]
filenames_verification_ai = ["aixflr5lc1.txt","aixflr5lc2.txt","aixflr5lc3.txt"]
filenames_ai_distribution = ["aivsblc1.png","aivsblc2.png","aivsblc3.png"]

# Creating a system of equations to solve for the Fourier coefficients for each of the load cases. 
for i in range(N_LC):
    print("Evaluating Load Case "+str(i+1))

    # Generating the Fourier coefficient matrix to solve the system. 
    M = np.zeros((N,K))
    for j in range(N):
        for k in range(K): 
            M[j,k] = np.sin((k+1)*theta[i][j])*(np.sin(theta[i][j]) + (k+1)*(chord_distribution[i][j]*clalpha[i]/(4*b[i])))
        
    # Inducing the propeller velocities onto the total velocity vector. 
    V_tot = np.zeros((N,3))
    for j in range(N): 
        V_tot[j,0] = V_cruise[i] + inducedVelocity[i][j][0]
        V_tot[j,2] = inducedVelocity[i][j][2] 
    
    # Creating a normalized vector for each N.
    V_norm = np.zeros((N,1))
    for j in range(N):
        V_norm_j = np.sqrt(V_tot[j][0]**2 +V_tot[j][1]**2 + V_tot[j][2]**2)
        V_norm[j][0] = V_norm_j

    # Generating the solution vector (b_vector)
    B = np.zeros((N,1))
    for j in range(N):
        B[j,0] = (clalpha[i]*chord_distribution[i][j]/(4*b[i]))*(V_cruise[i]/V_norm[j][0])*((alpha_geo[i][j] - alpha_l0[i])+V_tot[j][2]/V_norm[j][0])*np.sin(theta[i][j])
        
    # Solving the system and calculating the Fourier coefficients. Using a least squares solution to solve the system as N =! K. 
    A = lstsq(M,B,rcond=None)[0]

    # Creating a circulation distribution matrix 
    G = []
    for j in range(N):
        total = 0
        for k in range(K): 
            total = total + A[k]*np.sin(theta[i][j]*(k+1))
        total = total * 2 * b[i] * V_cruise[i] 
        G.append(total)

    # Computing the induced angle of attack distribution from the circulation matrix
    alpha_induced = []
    for j in range(N):
        total = 0 
        for k in range(K): 
            total = total + (k+1)*A[k]*(np.sin((k+1)*theta[i][j])/np.sin(theta[i][j]))
        alpha_induced.append(total)
    alpha_induced = np.array(alpha_induced)
    
    # 3D lift coefficient at this angle of attack. 
    C_L = float(A[0][0]*np.pi*AR[i])
        
    # Lift induced drag coefficient at this angle of attack. 
    delta = 0
    for j in range(1,K):
        delta = delta + (j+1)*(A[j][0]/A[0][0])**2
    C_Di = float(np.pi*AR[i]*A[0][0]**2 *(1+delta))

    e1 = C_L**2/(np.pi*C_Di*AR[i])
    e2 = 1/(1+delta)
    # Assert whether AR checks out
    AR_check = C_L**2/(np.pi*C_Di*e2)

    """VERIFICIATION STEP"""
    assert np.isclose(AR[i],float(AR_check),rtol=0.01), "Verification failed: Lifting line method not correctly implemented"
    
    # Computing the Clprandtl. This is a vector containing cl values normalized to the local V so not Vinfinity. 
    clprandtl = np.zeros((N,1))
    for j in range(N): 
        clprandtl[j,0] = 2*G[j][0] / (V_cruise[i] * chord_distribution[i][j])

    # Generating the lift distribution. Kutta-Jouwkowski theorem: Applying it with the freestream velocity instead of propeller induced velocity according to Dr. Sinnige. 
    L = np.zeros((N,1))
    for j in range(N):
        L[j,0] = rho_cruise[i]*V_cruise[i]*G[j][0]
    
  
    # Importing XFLR5's cl vs. y for verification purposes. 
    data_xflr_transposed_cl = np.genfromtxt(filenames_verification_cl[i])
    data_xflr_cl = np.transpose(data_xflr_transposed_cl)
    y_xflr_cl = data_xflr_cl[0]
    cl_xflr = data_xflr_cl[1]
    
    # Importing XFLR5's alpha_induced vs. y for verification purposes. 
    data_xflr_transposed_ai = np.genfromtxt(filenames_verification_ai[i])
    data_xflr_ai = np.transpose(data_xflr_transposed_ai)
    y_xflr_ai = np.array(data_xflr_ai[0])
    ai_xflr = np.array(data_xflr_ai[1])
    ai_xflr = np.deg2rad(ai_xflr)
    
    # Creating propeller plot arrays. 
    theta_p = np.linspace(0,2*np.pi,100)
    x_p_1 = -y_outer*b[i]/2+R_p*np.cos(theta_p)
    x_p_2 = -y_inner*b[i]/2+R_p*np.cos(theta_p)
    x_p_3 = y_outer*b[i]/2+R_p*np.cos(theta_p)
    x_p_4 = y_inner*b[i]/2+R_p*np.cos(theta_p)
    y_p = R_p*np.sin(theta_p)
    
    # Plotting the data. 
    fig = plt.figure(dpi=250)
    plt.plot(y[i],clprandtl, label = "Adapted LL model")
    if XFLR5 == True: 
        plt.scatter(y_xflr_cl,cl_xflr, label = "XFLR5 LL model",color = "indianred")
# =============================================================================
#     
#     plt.plot(x_p_1,y_p,color="black")
#     plt.plot(x_p_2,y_p,color="black")
#     plt.plot(x_p_3,y_p,color="black")
#     plt.plot(x_p_4,y_p,color="black")
#     plt.scatter(-y_outer*b[i]/2,0,color="black")
#     plt.scatter(-y_inner*b[i]/2,0,color="black")
#     plt.scatter(y_outer*b[i]/2,0,color="black")
#     plt.scatter(y_inner*b[i]/2,0,color="black")
#     plt.plot([-b[i]/2,b[i]/2],[0,0],"--",color="black")
# =============================================================================
    plt.xlabel("b [m]")
    plt.ylabel("$C_l$ [-]")
    plt.legend()
    fig.savefig(filenames_cl_distribution[i])
    
    fig2 = plt.figure(dpi =250)
    plt.plot(y[i],-alpha_induced,label = "Adapted LL model")
    if XFLR5 == True: 
        plt.scatter(y_xflr_ai,ai_xflr,label = "XFLR5 LL model",color = "indianred")
    plt.xlabel("b [m]")
    plt.ylabel(r"$\alpha_i$ [$\degree$]")
    plt.legend()
    fig2.savefig(filenames_ai_distribution[i])

    # Formatting the lift distribution for the SMM department. 
    lift_file = open(filenames_l_distribution[i],"w")
    lift_file.close()

    lift_file = open(filenames_l_distribution[i],"w")
    lines = []
    for j in range(N): 
        lines.append(str(L[j][0])+"\t"+str(y[i][j])+"\n")
    lift_file.writelines(lines)

# =============================================================================
# END
print("Propeller-wing interaction finished. See all created files in corresponding directory")













    
        
        






    

        
        
    












