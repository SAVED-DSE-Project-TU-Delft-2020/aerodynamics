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
from matplotlib import pyplot as plt
import seaborn as sns
sns.set() 

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
NACAdigit1 = round(cldes*20/3,0)

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
# Normalizing Cl/Cd according to a linear scoring function 
clcd = []
for i in range(len(airfoil_results)):
    clcd.append(airfoil_results[i][1])

dydx_clcd = compute_dydx(clcd)
b_clcd = compute_b(dydx_clcd,clcd)
 
clcd_normalised = []

for i in range(len(clcd)):
    clcd_normalised.append(clcd[i]*dydx_clcd+b_clcd)
    
# Normalizing Cmac according to a linear scoring function 
cmac = []
for i in range(len(airfoil_results)):
    cmac.append(airfoil_results[i][2])
 
dydx_cmac = compute_dydx(cmac)
b_cmac = compute_b(dydx_cmac,cmac)

cmac_normalised = []
 
for i in range(len(cmac)):
    cmac_normalised.append(cmac[i]*dydx_cmac+b_cmac)
     
# Normalizing Clmax according to a linear scoring function 
clmax = []
for i in range(len(airfoil_results)):
    clmax.append(airfoil_results[i][4])

dydx_clmax = compute_dydx(clmax)
b_clmax = compute_b(dydx_clmax,clmax)

clmax_normalised = []
 
for i in range(len(clmax)):
    clmax_normalised.append(clmax[i]*dydx_clmax+b_clmax)
    
# Normalizing alphastall according to a linear scoring function 
alphastall = []
for i in range(len(airfoil_results)):
    alphastall.append(airfoil_results[i][5])

dydx_alphastall = compute_dydx(alphastall)
b_alphastall = compute_b(dydx_alphastall,alphastall)

alphastall_normalised = []

for i in range(len(alphastall)):
    alphastall_normalised.append(alphastall[i]*dydx_alphastall+b_alphastall)
    
# Normalizing t/c according to a linear scoring function 
t_c = []
for i in range(len(airfoil_results)): 
    t_c.append(airfoil_results[i][6])

dydx_t_c = compute_dydx(t_c)
b_t_c = compute_b(dydx_t_c,t_c)

t_c_normalised = []
for i in range(len(t_c)):
    t_c_normalised.append(t_c[i]*dydx_t_c+b_t_c)
    
# Creating a normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 42. Number (N) of Criteria (C) = 5.
matrix_normalized_transposed = np.array([clcd_normalised,cmac_normalised,clmax_normalised,alphastall_normalised,t_c_normalised])
matrix_normalized = np.transpose(matrix_normalized_transposed)
 
# Assigning weights to the criteria (THESE HAVE BEEN SELECTED BASED ON ENGINEERING JUDGEMENT!). Weighted Product Model requires the sum of the weights to equal 1.
weights = [0.5,0.05,0.05,0.2,0.2]
criterion = ["Cl/Cd @ Cldes", "Cmac", "Clmax", "Alpha stall", "(t/c)_max"]


# Calculating the scores of each alternative based on the normalized Weighted Sum Model. 
P_WSM = WSM(matrix_normalized,weights)
index = winner_WSM(P_WSM)

print("The best airfoil according to the WSM is the "+airfoils[index][0])
airfoil_WSM = airfoils[index][0]

# WSM AD Sensitivity analysis 
sensitivity_matrix = sensitivity_WSM(matrix_normalized,weights)
critical_criterion_value,critical_criterion_index, alternative_1, alternative_2,  = find_critical_criterion(matrix_normalized,weights,sensitivity_matrix)

# Compiling a .txt file with all selected airfoil data for the other departments
final_file = open("finalresultsV2.txt", "w")
final_file.close()

final_file = open("finalresultsV2.txt","w")
final_file_lines = ["FINAL AIRFOIL SELECTION RESULTS","\n","\n",airfoil_WSM,"\n", "Cl/Cd @ Cldes = "+str(round(float(airfoil_results[index][1]),1)),"\n", "Cmac = "+str(round(float(airfoil_results[index][2]),4)),"\n","Alpha cruise = "+str(round(float(airfoil_results[index][3]),3)),"\n", "Cl_max @ Re_transition = "+str(round(float(airfoil_results[index][4]),1)),"\n","Alpha stall @ Re_transition = "+str(round(float(airfoil_results[index][5]),1)),"\n","(t/c)_max = "+str(airfoils[index][1])]
final_file.writelines(final_file_lines)
final_file_lines_sensitivity = ["\n","The most critical criterion is "+str(criterion[critical_criterion_index])+". This criterion should be changed by "+str(critical_criterion_value)+"% to swap the ranking of "+str(airfoils[alternative_1][0])+" and "+str(airfoils[alternative_2][0])]
final_file.writelines(final_file_lines_sensitivity)
final_file.close()

# =============================================================================
# Plotting the lift curve and the lift drag polar of the winning airfoil at both the cruise and stall Reynolds number 
plot_file_cruise = airfoils[index][0]+".txt"
plot_file_stall = airfoils[index][0]+"t.txt"

cruise_data_transposed = np.genfromtxt(plot_file_cruise)
cruise_data = np.transpose(cruise_data_transposed)
stall_data_transposed = np.genfromtxt(plot_file_stall)
stall_data = np.transpose(stall_data_transposed)

alpha_cruise = cruise_data[0][8:]
cl_cruise = cruise_data[1][8:]
cd_cruise = cruise_data[2][8:]

alpha_stall = stall_data[0][8:]
cl_stall = stall_data[1][8:]
cd_stall = stall_data[2][8:]

liftdragpolar = plt.figure(figsize = (10,5),dpi = 250)
plt.scatter(cd_cruise,cl_cruise, label = "Re = 860000",color = "blue",marker = "^")
plt.plot(cd_cruise,cl_cruise,color = "blue")
plt.scatter(cd_stall,cl_stall, label = "Re = 450000",color = "red")
plt.plot(cd_stall,cl_stall,color = "red")
plt.xlabel("$C_d$ [-]")
plt.ylabel("$C_l$ [-]")
plt.legend()
plt.savefig("liftdragpolar.png")

liftcurve = plt.figure(figsize = (10,5),dpi = 250)
plt.scatter(alpha_cruise,cl_cruise, label = "Re = 860000",color = "blue",marker = "^")
plt.plot(alpha_cruise,cl_cruise,color = "blue")
plt.scatter(alpha_stall,cl_stall, label = "Re = 450000",color = "red")
plt.plot(alpha_stall,cl_stall,color = "red")
plt.xlabel(r"$\alpha$  [-]")
plt.ylabel("$C_l$ [-]")
plt.legend()
plt.savefig("liftcurve.png")

# =============================================================================
# END
print("Airfoil selection finished, see results file in directory")





























