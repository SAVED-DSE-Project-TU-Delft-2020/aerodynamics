""" 
Title: PRELIMINARY AIRFOIL ANALYSIS DESIGN TOOL

This program analyses the performance of 44 NACA 5 digit airfoils. 
It stores the data of this analysis and assesses this data in a WPM and a WSM MCDM analysis. 
A sensitivity analysis of these MCDMs are then performed. 

- TODO: 
    Add verification of the MCDM methods (Based on the source)
    
Author: Casper Kanaar 

""" 
# =============================================================================
# Importing modules 
import numpy as np 
from isacalculator import compute_isa 
from scipy.interpolate import interp1d
import sys
sys.path.append('..')
from parameters import parameters

# =============================================================================
# Functions 
def compute_reynolds(rho,V,L,mu):
    return (rho*V*L)/mu

def compute_cldes(m,rho,V,Lambda):
    return (1.1*m*9.81)/(0.5*rho*(V*np.cos(Lambda))**2)

def compute_dydx(x):
    return (1-0)/(max(x)-min(x))

def compute_b(dydx,x):
    return -dydx*min(x)

def compute_caverage(c_r,c_t):
    return 0.5*(c_r+c_t)

# =============================================================================
#  Parameters 
V_cruise = parameters.Concept_C.Vcruise                                 # [m/s]. Cruise speed. FROZEN from FPP
V_transition = parameters.Concept_C.Vtransition                         # [m/s]. Transition speed. FROZEN from FPP
h_cruise = parameters.Concept_C.hcruise                                 # [m]. Cruise altitude. FROZEN from FPP
h_transition = parameters.Concept_C.htransition                         # [m]. Transition altitude. FROZEN from SIM
p_cruise,rho_cruise,T_cruise = compute_isa(h_cruise)                    # [Pa, kg/m3, K]. Cruise pressure, density, and temperature. FROZEN from FPP
p_transition,rho_transition,T_transition = compute_isa(h_transition)    # [Pa, kg/m3, K]. Cruise pressure, density, and temperature. FROZEN from SIM 
mu_cruise = 17.73e-6                                                    # [Pa.s]. Dynamic viscocity of air at cruise. FROZEN from FPP
mu_transition= 1.7955e-5                                                # [Pa.s]. Dynamic viscocity of air at transition. FROZEN from SIM 
Lambda_AD = parameters.Wing_A.sweepdegLErad                             # [rad]. Main wing sweep angle concept ABD. FROZEN from SMM
Lambda_C = parameters.Wing_C.sweepdegLErad                              # [rad]. Main wing sweep angle concept C. FROZEN from SMM
MTOM = parameters.Concept_A.MTOM                                        # [kg]. Maximum Take-Off Weight. FROZEN from BASELINE
c_r_AD = parameters.Wing_A.rootchord                                    # [m]. Root chord AD. FROZEN from SMM 
c_t_AD = parameters.Wing_A.tipchord                                     # [m]. Tip chord AD. FROZEN from SMM
b = parameters.Wing_C.span                                              # [m]. Wing span. FROZEN from project guide  
WL = parameters.Wing_C.wingloading                                      # [N/m2]. Wing loading. FROZEN from FPP. 
S = parameters.Wing_C.surfaceaera                                       # [m2]. Initial surface area estimate 
c_avg_C = parameters.Wing_C.averagechord                                # [m]. Average chord concept C. FROZEN from SMM

# =============================================================================
# Main

# Compute required parameters for concept A and D.
c_avg_AD = compute_caverage(c_r_AD,c_t_AD)
cldes_AD = compute_cldes(MTOM,rho_cruise,V_cruise,Lambda_AD)
Re_cruise_AD = compute_reynolds(rho_cruise,V_cruise,c_avg_AD,mu_cruise)
Re_transition_AD = compute_reynolds(rho_transition,V_transition,c_avg_AD,mu_transition)

# Compute required parameters for concept C.
cldes_C = compute_cldes(MTOM,rho_cruise,V_cruise,Lambda_C)
Re_cruise_C = compute_reynolds(rho_cruise,V_cruise,c_avg_C,mu_cruise)
Re_transition_C = compute_reynolds(rho_transition,V_transition,c_avg_C,mu_transition)

# =============================================================================
# Analysing all airfoil data 
 
# Clear file from last run
file_AD = open('airfoilanalysisresultsAD.txt', 'w')
file_C = open('airfoilanalysisresultsC.txt','w')
file_AD.close()
file_C.close()

# Opening the file to write to 
file_AD = open("airfoilanalysisresultsAD.txt","w")
file_C = open("airfoilanalysisresultsC.txt","w") 

# Create empty list to store results in 
airfoil_results_AD = []
file_results_AD = []
airfoil_results_C = []
file_results_C = []

# Creating list of airfoils for concept A, D and concept C
airfoils_AD = []
lst1_AD = ["321","331","341","351"]
lst2_AD = ["14","16","18","20"]
airfoils_C = []
lst1_C = ["310","320","330","340","350"]
lst2_C = ["10","12","14","16","18"]


for i in range(len(lst1_AD)):
    for j in range(len(lst2_AD)):
        airfoils_AD.append(lst1_AD[i]+lst2_AD[j])
        
for i in range(len(lst1_C)):
    for j in range(len(lst2_C)):
        airfoils_C.append(lst1_C[i]+lst2_C[j])
 
for i in range(len(airfoils_AD)):
    # Select airfoil
    airfoil = airfoils_AD[i]
   
    airfoil_data = np.genfromtxt(airfoil+".txt")
    airfoil_data_transposed = np.transpose(airfoil_data)
    alpha = airfoil_data_transposed[0]
    cl = airfoil_data_transposed[1]
    cm25 = airfoil_data_transposed[3]
    cl_cd = airfoil_data_transposed[8]

    cl_asfunctionof_alpha = interp1d(alpha,cl)
    alpha_asfunctionof_cl = interp1d(cl,alpha)
 
    alpha_at_cl0 = float(alpha_asfunctionof_cl(0))
    alpha_at_cldes = float(alpha_asfunctionof_cl(cldes_AD))
    alpha_cruise = alpha_at_cldes
 
    cm25_asfunctionof_alpha = interp1d(alpha,cm25)
    cmac = float(cm25_asfunctionof_alpha(alpha_at_cl0))
 
    cl_cd_asfunctionof_alpha = interp1d(alpha,cl_cd)
    cl_cd_at_cldes = cl_cd_asfunctionof_alpha(alpha_at_cldes)
 
    transition_airfoil_data = np.genfromtxt(airfoil+"t"+".txt")
    transition_airfoil_data_transposed = np.transpose(transition_airfoil_data)
    alpha_transition = transition_airfoil_data_transposed[0]
    cl_transition = transition_airfoil_data_transposed[1]
 
    alpha_asfunctionof_cl_transition = interp1d(cl_transition,alpha_transition)
    cl_asfunctionof_alpha = interp1d(alpha_transition,cl_transition)
    alphas_transition_sampled = np.linspace(alpha_transition[0],alpha_transition[-1],100)
    cl_alpha_samples = cl_asfunctionof_alpha(alphas_transition_sampled)
    cl_max = max(cl_alpha_samples)
    alpha_stall = alpha_asfunctionof_cl_transition(cl_max)

    airfoil_results_AD.append([airfoil,float(cl_cd_at_cldes),float(cmac),float(alpha_cruise),float(cl_max),float(alpha_stall)])

    result = "\n"+"NACA"+airfoil+": (Cl/Cd)_max @ Cldes = "+str(round(float(cl_cd_at_cldes),1))+", Cmac = "+str(round(float(cmac),4))+", alpha_cruise = "+str(round(float(alpha_cruise),2))+" , Cl_max @ transition = "+str(round(float(cl_max),2))+" ,alpha_stall @transition = "+str(round(float(alpha_stall),1))
    file_results_AD.append(result)
 
file_AD.writelines(file_results_AD)
file_AD.close()
 
 
for i in range(len(airfoils_C)):
    airfoil = airfoils_C[i]
    
    airfoil_data = np.genfromtxt(airfoil+".txt")
    airfoil_data_transposed = np.transpose(airfoil_data)
    alpha = airfoil_data_transposed[0]
    cl = airfoil_data_transposed[1]
    cm25 = airfoil_data_transposed[3]
    cl_cd = airfoil_data_transposed[8]

    cl_asfunctionof_alpha = interp1d(alpha,cl)
    alpha_asfunctionof_cl = interp1d(cl,alpha)
 
    alpha_at_cl0 = float(alpha_asfunctionof_cl(0))
    alpha_at_cldes = float(alpha_asfunctionof_cl(cldes_C))
    alpha_cruise = alpha_at_cldes
 
    cm25_asfunctionof_alpha = interp1d(alpha,cm25)
    cmac = float(cm25_asfunctionof_alpha(alpha_at_cl0))
 
    cl_cd_asfunctionof_alpha = interp1d(alpha,cl_cd)
    cl_cd_at_cldes = cl_cd_asfunctionof_alpha(alpha_at_cldes)
 
    transition_airfoil_data = np.genfromtxt(airfoil+"t"+".txt")
    transition_airfoil_data_transposed = np.transpose(transition_airfoil_data)
    alpha_transition = transition_airfoil_data_transposed[0]
    cl_transition = transition_airfoil_data_transposed[1]

    alpha_asfunctionof_cl_transition = interp1d(cl_transition,alpha_transition)
    cl_asfunctionof_alpha = interp1d(alpha_transition,cl_transition)
    alphas_transition_sampled = np.linspace(alpha_transition[0],alpha_transition[-1],100)
    cl_alpha_samples = cl_asfunctionof_alpha(alphas_transition_sampled)
    cl_max = max(cl_alpha_samples)
    alpha_stall = alpha_asfunctionof_cl_transition(cl_max)

    airfoil_results_C.append([airfoil,float(cl_cd_at_cldes),float(cmac),float(alpha_cruise),float(cl_max),float(alpha_stall)])
 
    result = "\n"+"NACA"+airfoil+": (Cl/Cd)_max @ Cldes = "+str(round(float(cl_cd_at_cldes),1))+", Cmac = "+str(round(float(cmac),4))+", alpha_cruise = "+str(round(float(alpha_cruise),2))+" , Cl_max @ transition = "+str(round(float(cl_max),2))+" ,alpha_stall @transition = "+str(round(float(alpha_stall),1))
    file_results_C.append(result)
     
file_C.writelines(file_results_C)
file_C.close()

# =============================================================================
# Performing an MCDM for concept A and D (as t/c is a criterion for A and D)
 
# Normalizing Cl/Cd according to a linear scoring function 
clcd_AD = []
for i in range(len(airfoil_results_AD)):
    clcd_AD.append(airfoil_results_AD[i][1])

dydx_clcd_AD = compute_dydx(clcd_AD)
b_clcd_AD = compute_b(dydx_clcd_AD,clcd_AD)
 
clcd_AD_normalised = []

for i in range(len(clcd_AD)):
    clcd_AD_normalised.append(clcd_AD[i]*dydx_clcd_AD+b_clcd_AD)
    
# Normalizing Cmac according to a linear scoring function 
cmac_AD = []
for i in range(len(airfoil_results_AD)):
    cmac_AD.append(airfoil_results_AD[i][2])
 
dydx_cmac_AD = compute_dydx(cmac_AD)
b_cmac_AD = compute_b(dydx_cmac_AD,cmac_AD)

cmac_AD_normalised = []
 
for i in range(len(cmac_AD)):
    cmac_AD_normalised.append(cmac_AD[i]*dydx_cmac_AD+b_cmac_AD)
     
# Normalizing Clmax according to a linear scoring function 
clmax_AD = []
for i in range(len(airfoil_results_AD)):
    clmax_AD.append(airfoil_results_AD[i][4])

dydx_clmax_AD = compute_dydx(clmax_AD)
b_clmax_AD = compute_b(dydx_clmax_AD,clmax_AD)

clmax_AD_normalised = []
 
for i in range(len(clmax_AD)):
    clmax_AD_normalised.append(clmax_AD[i]*dydx_clmax_AD+b_clmax_AD)
    
# Normalizing alphastall according to a linear scoring function 
alphastall_AD = []
for i in range(len(airfoil_results_AD)):
    alphastall_AD.append(airfoil_results_AD[i][5])

dydx_alphastall_AD = compute_dydx(alphastall_AD)
b_alphastall_AD = compute_b(dydx_alphastall_AD,alphastall_AD)

alphastall_AD_normalised = []

for i in range(len(alphastall_AD)):
    alphastall_AD_normalised.append(alphastall_AD[i]*dydx_alphastall_AD+b_alphastall_AD)

tc_AD = [14,16,18,20,14,16,18,20,14,16,18,20,14,16,18,20,] 
# Normalizing t/c according to exponential scoring function (self picked values)
tc_AD_normalised = [0,0.2,0.5,1.0,0,0.2,0.5,1.0,0,0.2,0.5,1.0,0,0.2,0.5,1.0]


# Creating a normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 16. Number (N) of Criteria (C) = 5.
matrix_AD_normalized_transposed = np.array([clcd_AD_normalised,cmac_AD_normalised,clmax_AD_normalised,alphastall_AD_normalised,tc_AD_normalised])
matrix_AD_normalized = np.transpose(matrix_AD_normalized_transposed)
 
# Assigning weights to the criteria (THESE HAVE BEEN SELECTED BASED ON ENGINEERING JUDGEMENT!). Weighted Product Model requires the sum of the weights to equal 1.
weights_AD = [0.4,0.1,0.1,0.05,0.35]

# Calculating the scores of each alternative based on the normalized Weighted Sum Model. 
scores_AD_WSM = []
for i in range(len(matrix_AD_normalized)):
    total = 0
    for j in range(len(matrix_AD_normalized[i])):
        total = total + matrix_AD_normalized[i][j]*weights_AD[j]
    scores_AD_WSM.append(total)

index_max_AD = scores_AD_WSM.index(max(scores_AD_WSM))
print("The best airfoil for Concept A and D according to the WSM is the NACA"+airfoils_AD[index_max_AD])
airfoil_AD_WSM = airfoils_AD[index_max_AD]

# WSM AD Sensitivity analysis 
sensitivity_matrix_AD_WSM = np.zeros((len(matrix_AD_normalized),len(matrix_AD_normalized[i])))
for i in range(len(scores_AD_WSM)-1):
    for j in range(len(weights_AD)):
        if matrix_AD_normalized[15][j] == matrix_AD_normalized[i][j]:
            sensitivity_matrix_AD_WSM[i,j] = 0
        else:
            sensitivity_matrix_AD_WSM[i,j] = (scores_AD_WSM[15] - scores_AD_WSM[i])/(matrix_AD_normalized[15][j]-matrix_AD_normalized[i][j]) * (100/weights_AD[j])
        
# Creating a non normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 20. Number (N) of Criteria (C) = 5
matrix_AD_transposed = np.array([clcd_AD,cmac_AD,clmax_AD,alphastall_AD,tc_AD])
matrix_AD = np.transpose(matrix_AD_transposed)
 
# Computing all ratios according to the Weighted Product Model.
ratio_matrix_AD = np.zeros((len(matrix_AD),len(matrix_AD)))
 
for i in range(len(ratio_matrix_AD)):
    for j in range(len(ratio_matrix_AD)): 
        total = 1
        for k in range(len(matrix_AD[0])): 
            total = total*((matrix_AD[i][k]/matrix_AD[j][k])**weights_AD[k])
        ratio_matrix_AD[i,j] = total

# Checking which airfoil is the best according to the Weighted Product Model
for i in range(len(ratio_matrix_AD)):
    count = 0 
    for j in range(len(ratio_matrix_AD[i])):
        if ratio_matrix_AD[i][j] >= 1: 
            count = count + 1 
    if count == 16:
        index = i
 
print("The best airfoil for Concept A and D according to the WPM is the NACA"+airfoils_AD[index])
airfoil_AD_WPM = airfoils_AD[index]

# WPM AD Sensitivity analysis. 
sensitivity_matrix_AD_WPM = np.zeros((len(matrix_AD),len(matrix_AD[0])))
for i in range(len(matrix_AD)-1):
    for j in range(len(matrix_AD[i])):
        total = 1
        for k in range(len(matrix_AD[i])):
            total = total*(matrix_AD[i][k]/matrix_AD[i+1][k])**weights_AD[k]
        K = np.log10(total)/np.log10(matrix_AD[i][j]/matrix_AD[i+1][j]) * 100 / weights_AD[j]
        sensitivity_matrix_AD_WPM[i,j] = K

# =============================================================================
# Performing an MCDM for concept C (as t/c is not a criterion for C)

# Normalizing Cl/Cd according to a linear scoring function 
clcd_C = []
for i in range(len(airfoil_results_C)):
    clcd_C.append(airfoil_results_C[i][1])

dydx_clcd_C = compute_dydx(clcd_C)
b_clcd_C = compute_b(dydx_clcd_C,clcd_C)

clcd_C_normalised = []

for i in range(len(clcd_C)):
    clcd_C_normalised.append(clcd_C[i]*dydx_clcd_C+b_clcd_C)
     
# Normalizing Clmax according to a linear scoring function 
clmax_C = []
for i in range(len(airfoil_results_C)):
    clmax_C.append(airfoil_results_C[i][4])

dydx_clmax_C = compute_dydx(clmax_C)
b_clmax_C = compute_b(dydx_clmax_C,clmax_C)

clmax_C_normalised = []
 
for i in range(len(clmax_C)):
    clmax_C_normalised.append(clmax_C[i]*dydx_clmax_C+b_clmax_C)
    
# Normalizing alphastall according to a linear scoring function 
alphastall_C = []
for i in range(len(airfoil_results_C)):
    alphastall_C.append(airfoil_results_C[i][5])

dydx_alphastall_C = compute_dydx(alphastall_C)
b_alphastall_C = compute_b(dydx_alphastall_C,alphastall_C)

alphastall_C_normalised = []

for i in range(len(alphastall_C)):
    alphastall_C_normalised.append(alphastall_C[i]*dydx_alphastall_C+b_alphastall_C)

# Creating a normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 20. Number (N) of Criteria (C) = 5
matrix_C_normalized_transposed = np.array([clcd_C_normalised,clmax_C_normalised,alphastall_C_normalised])
matrix_C_normalized = np.transpose(matrix_C_normalized_transposed)

# Assigning weights to the criteria (THESE HAVE BEEN SELECTED BASED ON ENGINEERING JUDGEMENT!). Weighted Product Model requires the sum of the weights to equal 1.
weights_C = [0.6,0.2,0.2]

# Calculating the scores of each alternative based on the normalized Weighted Sum Model. 
scores_C_WSM = []
for i in range(len(matrix_C_normalized)):
    total = 0
    for j in range(len(matrix_C_normalized[i])):
        total = total + matrix_C_normalized[i][j]*weights_C[j]
    scores_C_WSM.append(total)

index_max_C = scores_C_WSM.index(max(scores_C_WSM))
print("The best airfoil for Concept C according to the WSM is the NACA"+airfoils_C[index_max_C])
airfoil_C_WSM = airfoils_C[index_max_C]

# WSM AD Sensitivity analysis 
sensitivity_matrix_C_WSM = np.zeros((len(matrix_C_normalized),len(matrix_C_normalized[i])))
for i in range(len(scores_C_WSM)-1):
    for j in range(len(weights_C)):
        if matrix_C_normalized[24][j] == matrix_C_normalized[i][j]:
            sensitivity_matrix_C_WSM[i,j] = 0
        else:
            sensitivity_matrix_C_WSM[i,j] = (scores_C_WSM[24] - scores_C_WSM[i])/(matrix_C_normalized[24][j]-matrix_C_normalized[i][j]) * (100/weights_C[j])
        
# Creating a non normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 20. Number (N) of Criteria (C) = 5
matrix_C_transposed = np.array([clcd_C,clmax_C,alphastall_C])
matrix_C = np.transpose(matrix_C_transposed)

# Computing all ratios according to the Weighted Product Model.
ratio_matrix_C = np.zeros((len(matrix_C),len(matrix_C)))
 
for i in range(len(ratio_matrix_C)):
    for j in range(len(ratio_matrix_C)): 
        total = 1
        for k in range(len(matrix_C[0])): 
            total = total*((matrix_C[i][k]/matrix_C[j][k])**weights_C[k])
        ratio_matrix_C[i,j] = total

# Checking which airfoil is the best according to the Weighted Product Model
for i in range(len(ratio_matrix_C)):
    count = 0 
    for j in range(len(ratio_matrix_C[i])):
        if ratio_matrix_C[i][j] >= 1: 
            count = count + 1 
    if count == 25:
        index = i
 
print("The best airfoil for Concept C according to the WPM is the NACA"+airfoils_C[index])
airfoil_C_WPM = airfoils_C[index]

# WPM AD Sensitivity analysis. 
sensitivity_matrix_C_WPM = np.zeros((len(matrix_C),len(matrix_C[0])))
for i in range(len(matrix_C)-1):
    for j in range(len(matrix_C[i])):
        total = 1
        for k in range(len(matrix_C[i])):
            total = total*(matrix_C[i][k]/matrix_C[i+1][k])**weights_C[k]
        K = np.log10(total)/np.log10(matrix_C[i][j]/matrix_C[i+1][j]) * 100 / weights_C[j]
        sensitivity_matrix_C_WPM[i,j] = K
# =============================================================================
# Compiling a .txt file with all selected airfoil data for the other departments
final_file = open("finalresults.txt", "w")
final_file.close()

final_file = open("finalresults.txt","w")
final_file_lines_AD = ["FINAL AIRFOIL SELECTION RESULTS","\n","\n","Concept A and D: NACA"+airfoil_AD_WSM,"\n", "Cl/Cd @ Cldes = "+str(round(float(airfoil_results_AD[index_max_AD][1]),1)),"\n", "Cmac = "+str(round(float(airfoil_results_AD[index_max_AD][2]),4)),"\n","Alpha cruise = "+str(round(float(airfoil_results_AD[index_max_AD][2]),3)),"\n", "Cl_max @ Re_transition = "+str(round(float(airfoil_results_AD[index_max_AD][3]),1)),"\n","Alpha stall @ Re_transition = "+str(round(float(airfoil_results_AD[index_max_AD][4]),1))]
final_file.writelines(final_file_lines_AD)
final_file_lines_C = ["\n","\n","Concept C: NACA"+airfoil_C_WSM,"\n", "Cl/Cd @ Cldes = "+str(round(float(airfoil_results_C[index_max_C][1]),1)),"\n", "Cmac = "+str(round(float(airfoil_results_C[index_max_C][2]),4)),"\n","Alpha cruise = "+str(round(float(airfoil_results_C[index_max_C][2]),3)),"\n", "Cl_max @ Re_transition = "+str(round(float(airfoil_results_C[index_max_C][3]),1)),"\n","Alpha stall @ Re_transition = "+str(round(float(airfoil_results_C[index_max_C][4]),1))]
final_file.writelines(final_file_lines_C)
final_file.close()

# ==================================END========================================



        
        
            
 
        

            
         


   
    


    
    

    







