""" 
This program analyses the performance of 40 NACA 5 digit airfoils and performs a MCDM for all concept (A,B,C,D)

Author: Casper Kanaar 

""" 
# =============================================================================
# Importing modules 
import numpy as np 
from isacalculator import compute_isa 
from scipy.interpolate import interp1d

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

# =============================================================================
#  Parameters 
V_cruise = 28                                                           # [m/s]. Cruise speed
V_transition = 18                                                       # [m/s]. Transition speed 
h_cruise = 500                                                          # [m]. Cruise altitude 
h_transition = 20                                                       # [m]. Transition altitude 
p_cruise,rho_cruise,T_cruise = compute_isa(h_cruise)                    # [Pa, kg/m3, K]. Cruise pressure, density, and temperature.
p_transition,rho_transition,T_transition = compute_isa(h_transition)    # [Pa, kg/m3, K]. Cruise pressure, density, and temperature.
mu_cruise = 17.73e-6                                                    # [Pa.s]. Dynamic viscocity of air at cruise
mu_transition= 17.88e-6                                                 # [Pa.s]. Dynamic viscocity of air at transition 
Lambda = 0 / 180 * np.pi                                                # [rad]. Main wing sweep angle 
MTOM = 16.85                                                            # [kg]. Maximum Take-Off Weight
c_avg = 0.3                                                             # [m]. Average chord length

# =============================================================================
# Main

# Compute design lift coefficient 
cldes = compute_cldes(MTOM,rho_cruise,V_cruise,Lambda)

# Compute cruise Reynolds number 
Re_cruise = compute_reynolds(rho_cruise,V_cruise,c_avg,mu_cruise)

# Compute transition Reynolds number
Re_transition = compute_reynolds(rho_transition,V_transition,c_avg,mu_transition)

# =============================================================================
# Analysing all airfoil data 

# Clear file from last run
file_ABD = open('airfoilanalysisresultsABD.txt', 'w')
file_C = open('airfoilanalysisresultsC.txt','w')
file_ABD.close()
file_C.close()

# Opening the file to write to 
file_ABD = open("airfoilanalysisresultsABD.txt","w")
file_C = open("airfoilanalysisresultsC.txt","w") 

# Create empty list to store results in 
airfoil_results_ABD = []
file_results_ABD = []
airfoil_results_C = []
file_results_C = []

# Creating list of airfoils for concept A, B, D and concept C
airfoils_ABD = []
lst1_ABD = ["2211","2311","2411","2511"]
lst2_ABD = ["0","2","4","6","8"]
airfoils_C = []
lst1_C = ["2101","2201","2301","2401","2501"]
lst2_C = ["0","2","4","6"]


for i in range(len(lst1_ABD)):
    for j in range(len(lst2_ABD)):
        airfoils_ABD.append("naca"+lst1_ABD[i]+lst2_ABD[j])
        
for i in range(len(lst1_C)):
    for j in range(len(lst2_C)):
        airfoils_C.append("naca"+lst1_C[i]+lst2_C[j])

for i in range(len(airfoils_ABD)):
    # Select airfoil
    airfoil = airfoils_ABD[i]
    
    airfoil_data = np.genfromtxt(airfoil+".txt")
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

    transition_airfoil_data = np.genfromtxt(airfoil+"transition"+".txt")
    transition_airfoil_data_transposed = np.transpose(transition_airfoil_data)
    alpha_transition = transition_airfoil_data_transposed[0]
    cl_transition = transition_airfoil_data_transposed[1]

    alpha_asfunctionof_cl_transition = interp1d(cl_transition,alpha_transition)
    cl_asfunctionof_alpha = interp1d(alpha_transition,cl_transition)
    alphas_transition_sampled = np.linspace(alpha_transition[0],alpha_transition[-1],100)
    cl_alpha_samples = cl_asfunctionof_alpha(alphas_transition_sampled)
    cl_max = max(cl_alpha_samples)
    alpha_stall = alpha_asfunctionof_cl_transition(cl_max)

    airfoil_results_ABD.append([airfoil,float(cl_cd_at_cldes),float(cmac),float(alpha_cruise),float(cl_max),float(alpha_stall)])

    result = "\n"+airfoil+": (Cl/Cd)_max @ Cldes = "+str(round(float(cl_cd_at_cldes),1))+", Cmac = "+str(round(float(cmac),4))+", alpha_cruise = "+str(round(float(alpha_cruise),2))+" , Cl_max @ transition = "+str(round(float(cl_max),2))+" ,alpha_stall @transition = "+str(round(float(alpha_stall),1))
    file_results_ABD.append(result)

file_ABD.writelines(file_results_ABD)
file_ABD.close()


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
    alpha_at_cldes = float(alpha_asfunctionof_cl(cldes))
    alpha_cruise = alpha_at_cldes

    cm25_asfunctionof_alpha = interp1d(alpha,cm25)
    cmac = float(cm25_asfunctionof_alpha(alpha_at_cl0))

    cl_cd_asfunctionof_alpha = interp1d(alpha,cl_cd)
    cl_cd_at_cldes = cl_cd_asfunctionof_alpha(alpha_at_cldes)

    transition_airfoil_data = np.genfromtxt(airfoil+"transition"+".txt")
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

    result = "\n"+airfoil+": (Cl/Cd)_max @ Cldes = "+str(round(float(cl_cd_at_cldes),1))+", Cmac = "+str(round(float(cmac),4))+", alpha_cruise = "+str(round(float(alpha_cruise),2))+" , Cl_max @ transition = "+str(round(float(cl_max),2))+" ,alpha_stall @transition = "+str(round(float(alpha_stall),1))
    file_results_C.append(result)
    
file_C.writelines(file_results_C)
file_C.close()

# =============================================================================
# Performing an MCDM for concept A and D (as t/c is a criterion for A and D)

# Normalizing Cl/Cd according to a linear scoring function 
clcd_ABD = []
for i in range(len(airfoil_results_ABD)):
    clcd_ABD.append(airfoil_results_ABD[i][1])

dydx_clcd_ABD = compute_dydx(clcd_ABD)
b_clcd_ABD = compute_b(dydx_clcd_ABD,clcd_ABD)

clcd_ABD_normalised = []

for i in range(len(clcd_ABD)):
    clcd_ABD_normalised.append(clcd_ABD[i]*dydx_clcd_ABD+b_clcd_ABD)
    
# Normalizing Cmac according to a linear scoring function 
cmac_ABD = []
for i in range(len(airfoil_results_ABD)):
    cmac_ABD.append(airfoil_results_ABD[i][2])

dydx_cmac_ABD = compute_dydx(cmac_ABD)
b_cmac_ABD = compute_b(dydx_cmac_ABD,cmac_ABD)

cmac_ABD_normalised = []

for i in range(len(cmac_ABD)):
    cmac_ABD_normalised.append(cmac_ABD[i]*dydx_cmac_ABD+b_cmac_ABD)
    
# Normalizing Clmax according to a linear scoring function 
clmax_ABD = []
for i in range(len(airfoil_results_ABD)):
    clmax_ABD.append(airfoil_results_ABD[i][4])

dydx_clmax_ABD = compute_dydx(clmax_ABD)
b_clmax_ABD = compute_b(dydx_clmax_ABD,clmax_ABD)

clmax_ABD_normalised = []

for i in range(len(clmax_ABD)):
    clmax_ABD_normalised.append(clmax_ABD[i]*dydx_clmax_ABD+b_clmax_ABD)
    
# Normalizing alphastall according to a linear scoring function 
alphastall_ABD = []
for i in range(len(airfoil_results_ABD)):
    alphastall_ABD.append(airfoil_results_ABD[i][5])

dydx_alphastall_ABD = compute_dydx(alphastall_ABD)
b_alphastall_ABD = compute_b(dydx_alphastall_ABD,alphastall_ABD)

alphastall_ABD_normalised = []

for i in range(len(alphastall_ABD)):
    alphastall_ABD_normalised.append(alphastall_ABD[i]*dydx_alphastall_ABD+b_alphastall_ABD)

tc_ABD = [10,12,14,16,18,10,12,14,16,18,10,12,14,16,18,10,12,14,16,18] 
# Normalizing t/c according to exponential scoring function (self picked values)
tc_ABD_normalised = [0,0.1,0.25,0.6,1,0,0.1,0.25,0.6,1,0,0.1,0.25,0.6,1,0,0.1,0.25,0.6,1]


# Creating a normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 20. Number (N) of Criteria (C) = 5
matrix_AD_normalized_transposed = np.array([clcd_ABD_normalised,cmac_ABD_normalised,clmax_ABD_normalised,alphastall_ABD_normalised,tc_ABD_normalised])
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
print("The best airfoil for Concept A and D according to the WSM is the "+airfoils_ABD[index_max_AD])
airfoil_AD_WSM = airfoils_ABD[index_max_AD]

# WSM AD Sensitivity analysis 
sensitivity_matrix_AD_WSM = np.zeros((len(matrix_AD_normalized),len(matrix_AD_normalized[i])))
for i in range(len(scores_AD_WSM)-1):
    for j in range(len(weights_AD)):
        if matrix_AD_normalized[19][j] == matrix_AD_normalized[i][j]:
            sensitivity_matrix_AD_WSM[i,j] = 0
        else:
            sensitivity_matrix_AD_WSM[i,j] = (scores_AD_WSM[19] - scores_AD_WSM[i])/(matrix_AD_normalized[19][j]-matrix_AD_normalized[i][j]) * (100/weights_AD[j])
        
# Creating a non normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 20. Number (N) of Criteria (C) = 5
matrix_AD_transposed = np.array([clcd_ABD,cmac_ABD,clmax_ABD,alphastall_ABD,tc_ABD])
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
    if count == 20:
        index = i

print("The best airfoil for Concept A and D according to the WPM is the "+airfoils_ABD[index])
airfoil_AD_WPM = airfoils_ABD[index]

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
# Performing an MCDM for concept B (as t/c is not a criterion for B).
        
# Creating a normalized Measure of Performance Matrx for B. Number (M) of alternatives (A) = 20. Number (N) of Criteria (C) = 4
matrix_B_normalized_transposed = np.array([clcd_ABD_normalised,cmac_ABD_normalised,clmax_ABD_normalised,alphastall_ABD_normalised])
matrix_B_normalized = np.transpose(matrix_B_normalized_transposed)

# Assigning weights to the criteria (THESE HAVE BEEN SELECTED BASED ON ENGINEERING JUDGEMENT!). Weighted Product Model requires the sum of the weights to equal 1.
weights_B = [0.5,0.2,0.2,0.1]

# Calculating the scores of each alternative based on the normalized Weighted Sum Model. 
scores_B_WSM = []
for i in range(len(matrix_B_normalized)):
    total = 0
    for j in range(len(matrix_B_normalized[i])):
        total = total + matrix_B_normalized[i][j]*weights_B[j]
    scores_B_WSM.append(total)

index_max_B = scores_B_WSM.index(max(scores_B_WSM))
print("The best airfoil for Concept B according to the WSM is the "+airfoils_ABD[index_max_B])
airfoil_B_WSM = airfoils_ABD[index_max_B]

# WSM B Sensitivity analysis 
sensitivity_matrix_B_WSM = np.zeros((len(matrix_B_normalized),len(matrix_B_normalized[i])))
for i in range(len(scores_B_WSM)-1):
    for j in range(len(weights_B)):
        if matrix_B_normalized[19][j] == matrix_B_normalized[i][j]:
            sensitivity_matrix_B_WSM[i,j] = 0
        else:
            sensitivity_matrix_B_WSM[i,j] = (scores_B_WSM[19] - scores_B_WSM[i])/(matrix_B_normalized[19][j]-matrix_B_normalized[i][j]) * (100/weights_B[j])

# Creating a non normalized Measure of Performance Matrx for A and D. Number (M) of alternatives (A) = 20. Number (N) of Criteria (C) = 5
matrix_B_transposed = np.array([clcd_ABD,cmac_ABD,clmax_ABD,alphastall_ABD])
matrix_B = np.transpose(matrix_B_transposed)

# Computing all ratios according to the Weighted Product Model.
ratio_matrix_B = np.zeros((len(matrix_B),len(matrix_B)))
 
for i in range(len(ratio_matrix_B)):
    for j in range(len(ratio_matrix_B)): 
        total = 1
        for k in range(len(matrix_B[0])): 
            total = total*((matrix_B[i][k]/matrix_B[j][k])**weights_B[k])
        ratio_matrix_B[i,j] = total

# Checking which airfoil is the best according to the Weighted Product Model
for i in range(len(ratio_matrix_B)):
    count = 0 
    for j in range(len(ratio_matrix_B[i])):
        if ratio_matrix_B[i][j] >= 1: 
            count = count + 1 
    if count == 20:
        index = i

print("The best airfoil for Concept B according to the WPM is the "+airfoils_ABD[index])
airfoil_B_WPM = airfoils_ABD[index]

# WPM AD Sensitivity analysis. 
sensitivity_matrix_B_WPM = np.zeros((len(matrix_B),len(matrix_B[0])))
for i in range(len(matrix_B)-1):
    for j in range(len(matrix_B[i])):
        total = 1
        for k in range(len(matrix_B[i])):
            total = total*(matrix_B[i][k]/matrix_B[i+1][k])**weights_B[k]
        K = np.log10(total)/np.log10(matrix_B[i][j]/matrix_B[i+1][j]) * 100 / weights_B[j]
        sensitivity_matrix_B_WPM[i,j] = K

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
weights_C = [0.5,0.3,0.2]

# Calculating the scores of each alternative based on the normalized Weighted Sum Model. 
scores_C_WSM = []
for i in range(len(matrix_C_normalized)):
    total = 0
    for j in range(len(matrix_C_normalized[i])):
        total = total + matrix_C_normalized[i][j]*weights_C[j]
    scores_C_WSM.append(total)

index_max_C = scores_C_WSM.index(max(scores_C_WSM))
print("The best airfoil for Concept C according to the WSM is the "+airfoils_C[index_max_C])
airfoil_C_WSM = airfoils_C[index_max_C]

# WSM AD Sensitivity analysis 
sensitivity_matrix_C_WSM = np.zeros((len(matrix_C_normalized),len(matrix_C_normalized[i])))
for i in range(len(scores_C_WSM)-1):
    for j in range(len(weights_C)):
        if matrix_C_normalized[19][j] == matrix_C_normalized[i][j]:
            sensitivity_matrix_C_WSM[i,j] = 0
        else:
            sensitivity_matrix_C_WSM[i,j] = (scores_C_WSM[19] - scores_C_WSM[i])/(matrix_C_normalized[19][j]-matrix_C_normalized[i][j]) * (100/weights_C[j])
        
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
    if count == 20:
        index = i

print("The best airfoil for Concept C according to the WPM is the "+airfoils_C[index])
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



        
        
            
 

        

            
         


   
    


    
    

    







