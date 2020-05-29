# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:48:18 2020

@author: Casper

Detailed design phase main airfoil selection. 

ASSUMPTIONS:
    - The entire aircraft generates lift. 
    
"""

# =============================================================================
# Importing relevant modules 
import numpy as np 
from isacalculator import compute_isa 
from scipy.interpolate import interp1d

# =============================================================================
# Functions 
def compute_reynolds(rho,V,L,mu):
    return (rho*V*L)/mu

def compute_cldes(m,rho,V,Lambda):
    return (1.1*m*9.81)/(S*0.5*rho*(V*np.cos(Lambda))**2)

def compute_dydx(x):
    return (1-0)/(max(x)-min(x))

def compute_b(dydx,x):
    return -dydx*min(x)

def compute_caverage(c_r,c_t):
    return 0.5*(c_r+c_t)

def WSM(matrix,weights):
    P = []
    for i in range(len(matrix)):
        total = 0
        for j in range(len(matrix[i])):
            total = total + matrix[i,j]*weights[j]
        P.append(total)
    return P

def winner_WSM(P): 
    index = P.index(max(P))
    return index 

def sensitivity_WSM(matrix,weights):
    P = WSM(matrix,weights)
    sensitivity_matrix = []
    for i in range(len(matrix)-1):
        for j in range(i + 1,len(matrix)):
            A_i = matrix[i]
            A_j = matrix[j]
            sensitivity_matrix_row = []
            for k in range(len(matrix[i])):
                if A_j[k]!=A_i[k]:
                    param = (P[j] - P[i])/(A_j[k]-A_i[k])
                else:
                    param = np.inf
                
                if param >= weights[k]:
                    delta = "N.F."
                    sensitivity_matrix_row.append(delta)
                else:
                    delta = (P[j] - P[i])/(A_j[k]-A_i[k])*(100/weights[k])
                    sensitivity_matrix_row.append(delta)
        
            sensitivity_matrix.append(sensitivity_matrix_row)
                    
    return sensitivity_matrix

def find_critical_criterion(matrix,weights,sensitivity):
    P = WSM(matrix,weights)
    index_winner = winner_WSM(P)
    
    j = 1
    transition_list = [0]
    transition = 0
    total = len(sensitivity)
    n_alternatives = len(matrix)
    
    while total > 0: 
        total = total - (n_alternatives-j)
        transition = transition + (n_alternatives - j)
        transition_list.append(transition)
        j = j + 1 
    
    index_list = []
    for i in range(len(sensitivity)):
        row = i
        for k in range(len(transition_list)-1):
            if transition_list[k] <= row < transition_list[k+1]:
                index_1 = k 
        index_2 = row - transition_list[index_1] + index_1 + 1
        index_list.append([index_1,index_2])
    
    save_indexes = []
    for i in range(len(index_list)):
        if index_list[i][0] == index_winner or index_list[i][1] == index_winner:
            save_indexes.append(i)        
    value = 100
    
    for i in range(len(save_indexes)): 
        for j in range(len(sensitivity[0])): 
            if type(sensitivity[save_indexes[i]][j]) == str:
                 continue
            elif abs(sensitivity[save_indexes[i]][j]) <= abs(value): 
                value = sensitivity[save_indexes[i]][j]
                criterion_index = j 
                alternative_1_index = index_list[save_indexes[i]][0]
                alternative_2_index = index_list[save_indexes[i]][1]
                
    return value, criterion_index, alternative_1_index, alternative_2_index  
# =============================================================================
# Parameters
m = 17.49525
W_S = 122.23
S = (m*9.81)/W_S
V_cruise = 28
V_transition = 14
b = 3 
c_mgc = S/b 
h_transition = 20 
h_cruise = 500 
p_transition,rho_transition,T_transition = compute_isa(h_transition)   
p_cruise,rho_cruise,T_cruise = compute_isa(h_cruise)   
mu_cruise = 17.73e-6                                                    
mu_cruise = 17.73e-6                                                  
mu_transition= 1.7955e-5 
sweep = 25* np.pi / 180

# =============================================================================
# MAIN PROGRAM 

# Computing the Reynolds numbers at cruise and transition 
Re_cruise = compute_reynolds(rho_cruise,V_cruise,c_mgc,mu_cruise)
Re_transition = compute_reynolds(rho_transition,V_transition,c_mgc,mu_cruise)

# Computing the design lift coefficient 
cldes = compute_cldes(m,rho_cruise,V_cruise,sweep)
NACAdigit1 = round(Cldes*20/3,0)

# Clear file from last run
file = open('airfoilanalysisresultsdetailed.txt', 'w')
file.close()

# Opening the file to write to 
file = open("airfoilanalysisresultsdetailed.txt","w")

# Creating a selection of airfoils to be analysed. 
airfoils = [] 
low_speed = [["ag12",6.24],["ag16",7.11],["ag24",8.41],["ag35",8.72],["cal1215j",11.72],["cal2263m",11.72],["cal4014l",10],["e231",12.33],["e374",10.91],["e387",9.07],["rg15",8.92],["s7012",8.75],["s8064",12.33],["s9000",9.01],["sa7035",9.19],["sa7036",9.20],["sd7037",9.20],["sd7080",9.15]]
naca = [["22106",6],["23106",6],["24106",6],["25106",6],["22108",8],["23108",8],["24108",8],["25108",8],["22110",10],["23110",10],["24110",10],["25110",10],["22112",12],["23112",12],["24112",12],["25112",12],["22114",14],["23114",14],["24114",14],["25114",14],["22116",16],["23116",16],["24116",16],["25116",16]]
for i in range(len(naca)):
    airfoils.append(naca[i])
for i in range(len(low_speed)):
    airfoils.append(low_speed[i])
    
    
airfoil_results = []
file_results = []
    
# Analysing the performance data of the airfoils 
for i in range(len(airfoils)): 
    airfoil = airfoils[i]

    airfoil_data = np.genfromtxt(airfoil[0]+".txt")
    airfoil_data_transposed = np.transpose(airfoil_data)
    
    alpha = airfoil_data_transposed[0]
    cl = airfoil_data_transposed[1]
    cm25 = airfoil_data_transposed[3]
    cl_cd = airfoil_data_transposed[8]
    
    cl_asfunctionof_alpha = interp1d(alpha,cl)
    alpha_asfunctionof_cl = interp1d(cl,alpha)
 
    alpha_at_cl0 = float(alpha_asfunctionof_cl(0))
    alpha_at_cldes = float(alpha_asfunctionof_cl(cldes))
    alpha_cruise = alpha_at_cldes
 
    cm25_asfunctionof_alpha = interp1d(alpha,cm25)
    cmac = float(cm25_asfunctionof_alpha(alpha_at_cl0))
 
    cl_cd_asfunctionof_alpha = interp1d(alpha,cl_cd)
    cl_cd_at_cldes = cl_cd_asfunctionof_alpha(alpha_at_cldes)
 
    transition_airfoil_data = np.genfromtxt(airfoil[0]+"t"+".txt")
    transition_airfoil_data_transposed = np.transpose(transition_airfoil_data)
    alpha_transition = transition_airfoil_data_transposed[0]
    cl_transition = transition_airfoil_data_transposed[1]
 
    alpha_asfunctionof_cl_transition = interp1d(cl_transition,alpha_transition)
    cl_asfunctionof_alpha = interp1d(alpha_transition,cl_transition)
    alphas_transition_sampled = np.linspace(alpha_transition[0],alpha_transition[-1],100)
    cl_alpha_samples = cl_asfunctionof_alpha(alphas_transition_sampled)
    cl_max = max(cl_alpha_samples)
    alpha_stall = alpha_asfunctionof_cl_transition(cl_max)
    
    airfoil_results.append([airfoil[0],float(cl_cd_at_cldes),float(cmac),float(alpha_cruise),float(cl_max),float(alpha_stall),airfoil[1]])
    
    result = "\n"+airfoil[0]+": (Cl/Cd)_max @ Cldes = "+str(round(float(cl_cd_at_cldes),1))+", Cmac = "+str(round(float(cmac),4))+", alpha_cruise = "+str(round(float(alpha_cruise),2))+", Cl_max @ transition = "+str(round(float(cl_max),2))+", alpha_stall @transition = "+str(round(float(alpha_stall),1))+" ,(t/c)_max = "+str(round(float(airfoil[1]),1))
    file_results.append(result)
    
file.writelines(file_results)
file.close()


# =============================================================================
# Perform a MCDM. The WSM will be used     

    
    


 































