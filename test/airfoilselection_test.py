# -*- coding: utf-8 -*-
"""
This file tests the airfoil selection code. 

@author: Casper Kanaar
"""
# =============================================================================
# Import relevant modules 
import numpy as np 
import sys 
sys.path.append('..')

# =============================================================================
# Defining test parameters 
matrix_test = np.array([[0.3088,0.2897,0.3867,0.1922],[0.2163,0.3458,0.1755,0.6288],[0.4509,0.2473,0.1194,0.0575],[0.0240,0.1172,0.3184,0.1215]])
weights_test = np.array([0.3277,0.3058,0.2876,0.0790])

matrix_test_WPM = np.array([[0.9381,0.3501,0.8811,0.5646],[0.7691,0.4812,0.1679,0.9336],[0.9445,0.1138,0.2219,0.0135],[0.1768,0.0221,0.9462,0.1024]])
weights_test_WPM = np.array([0.4504,0.1231,0.0848,0.3417])

# =============================================================================
# Defining test functions 
# =============================================================================
# Functions 
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
    countlist = []
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
        print(index_1,index_2)
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

# def WPM(matrix,weights):
#     P = []
#     countlist = []
#     for i in range(len(matrix)):
#         for j in range(len(matrix)):
#             if i == j: 
#                 continue
#             else:
#                 lst = [i,j]
#                 if countlist.count([j,i]) == 1:
#                     countlist.append(lst)
#                     continue
#                 else:
#                     A_i = matrix[i]
#                     A_j = matrix[j]
#                     ratio = 1
#                     for k in range(len(A_i)):
#                         ratio = ratio * (A_i[k]/A_j[k])**weights[k]
#                     P.append(ratio)
#    
#   
#     return P, ratio_matrix 
#            
# def winner_WPM(matrix_test_WPM,weights_test_WPM):
#     P, ratio_matrix = WPM(matrix_test_WPM,weights_test_WPM)
#     for i in range(len(ratio_matrix)): 
#         count = 0 
#         for j in range(len(ratio_matrix[i])):
#             if ratio_matrix[i][j] >= 1: 
#                 count = count + 1 
#         if count == len(ratio_matrix):
#             index = i
#     return index 
# =============================================================================
        
# Testing the WSM 
def test_WSM(matrix_test,weights_test):
    P_test = np.array([0.3162,0.2768,0.2621,0.1449])
    P_results = np.array(WSM(matrix_test,weights_test)) 
    for i in range(len(P_test)):
        assert np.isclose(P_test.all(),P_results.all(),rtol = 0.001), "WSM not correctly implemented"
    return print("WSM test completed")

# Testing the winner function of the WSM 
def test_winner_WSM(matrix_test,weights_test): 
    P = WSM(matrix_test,weights_test)
    winner = P[winner_WSM(P)]
    lst1 = np.array([0.3162])
    lst2 = np.array([winner])
    assert np.isclose(lst1.all(),lst2.all(),rtol=0.001), "Winner function WSM not correctly implemented"
    return print("Winner WSM test completed")

# Testing the sensitivty analysis of the WSM
def test_sensitivity_WSM(matrix_test,weights_test): 
    sensitivity_matrix = sensitivity_WSM(matrix_test,weights_test)
    test_sensitivity_matrix = [["N.F.",-229.7,64.8818,-114.1772],[-116.1733,"N.F.",70.35,"N.F."],["N.F.","N.F.","N.F.","N.F."],[-19.1334,48.7901,9.1099,32.5317],["N.F.","N.F.",-320.9,"N.F."],[83.7656,"N.F.",-204.8,2318.10]]
    assert np.shape(sensitivity_matrix) == np.shape(test_sensitivity_matrix), "Sensitivity matrix not the correct shape"
    for i in range(len(sensitivity_matrix)):
        for j in range(len(sensitivity_matrix[i])):
            if type(sensitivity_matrix[i][j]) == str:
                continue 
            else:
                lst1 = np.array([sensitivity_matrix[i][j]])
                lst2 = np.array([test_sensitivity_matrix[i][j]])
                assert np.isclose(lst1.all(),lst2.all(),rtol = 0.001), "WSM Sensitivity matrix not correctly implemented"
    return print("WSM sensitivity test completed")

# =============================================================================
# def test_WPM(matrix_test_WPM,weights_test_WPM): 
#     solution = np.array([1.0192,4.6082,5.3062,4.5216,5.2065,1.1515])
#     P,ratio_matrix = WPM(matrix_test_WPM,matrix_test_WPM)
#     P = np.array(P)
#     assert np.isclose(solution.all(),P.all(),rtol = 0.010), "WPM not correctly implemented"
#     return print("WPM test completed")
# =============================================================================

# =============================================================================
# def test_winner_WPM(matrix_test_WPM,weights_test_WPM):
#     winner = [0]
#     winner_test = [winner_WPM(matrix_test_WPM,weights_test_WPM)]
#     winner = np.array(winner)
#     winner_test = np.array(winner_test)
#     assert np.isclose(winner,winner_test,atol=0), "Winner WPM not correctly implemented"
#     return print("Winner WPM test completed")
# =============================================================================
    
# =============================================================================
# TESTING 
test_WSM(matrix_test,weights_test)
test_winner_WSM(matrix_test,weights_test)
test_sensitivity_WSM(matrix_test,weights_test)

matrix = matrix_test
weights = weights_test
sensitivity = sensitivity_WSM(matrix,weights)
    
 
    

    
    


        

