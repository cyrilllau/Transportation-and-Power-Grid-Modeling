#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 00:54:52 2020

@author: Cyril
"""

import numpy as np
from docplex.mp.model import Model
import pandas as pd
import math
import csv

#########################################################
#################### SET UP DATA ########################
#########################################################

##########Fixed Cost and Charging Load cost#################

A = pd.read_excel('Networks.xlsx','LocationData').values[:,1]
C = pd.read_excel('Networks.xlsx','LocationData').values[:,2]



##########PDN Line para###########

#initialize PDN list
PDN_line_mat = pd.read_excel('Networks.xlsx','PDN_line').values
PDN_line_list = []
for i in range(PDN_line_mat.shape[0]):
    PDN_line_list.append((int(PDN_line_mat[i][0]),int(PDN_line_mat[i][1])))
    
#Resistance, Reactance Capacity Power flow
R_dict = {}
X_dict = {}
Pmax = {}
Qmax = {}
G_dict = {}
for i in range(PDN_line_mat.shape[0]):
    R_dict[PDN_line_list[i]] = PDN_line_mat[i][2]
    X_dict[PDN_line_list[i]] = PDN_line_mat[i][3]
    Pmax[PDN_line_list[i]] = PDN_line_mat[i][4]
    Qmax[PDN_line_list[i]] = PDN_line_mat[i][5]
    G_dict[PDN_line_list[i]] = PDN_line_mat[i][6]


############### PDN Node #################
    
#PDN Node with scenario
P_load_mat = pd.read_excel('Networks.xlsx','PDN_node_Pload').values[:,1:]
Q_load_mat = pd.read_excel('Networks.xlsx','PDN_node_Qload').values[:,1:]


#PDN Node without scenario
W = pd.read_excel('Networks.xlsx','PDN_node_constant').values[:,1]
F = pd.read_excel('Networks.xlsx','PDN_node_constant').values[:,2]




########## TN Node #################
D_mat = pd.read_excel('Networks.xlsx','Demand_TN').values[:,1:][0:]


########## Utility #############
Utility  = pd.read_excel('Networks.xlsx','Utility').values[:,1:]


############# Some Constant Parameter #####################
#Penalty
H = 20
G = 1
#Minimum Utility
Umin = 20
#Voltage bound (need be fixed later)
Vmax = 100000000000000000000
Vmin = 0
#Parameter in the Constrains
Ud= 50
K = 10
M = 500000

##########################
########  Weight #########
##########################

alpha= 0.7
beta = 0.5  #H
gamma = 0.5 #G
phi = 0.5





################# Numbers For Iterarion ###################
I = D_mat.shape[0]    #Number of TN Nodes
#ndata = D_mat.shape[1]  #Number of Scenario
ndata = 5
J = A.shape[0]    #Number of Openning Location
J_plus = J+1    
#N = P_load_mat.shape[0]     #Number of PDN Nodes



###############################################################
#################### Solve the Problem ########################
###############################################################
''' Solve the Problem '''


var_list_location = [1,2,4,5,10,11,13,14,15,19,20,21]
var_list_Jplus = [0,1,2,4,5,10,11,13,14,15,19,20,21]
var_list_PDN = [2,4,6,16,25,20,3,8,23,31,7,9]
var_list_TN = [i for i in range(I)]
var_list_K = [i for i in range(K)]
var_list_PDN_node = [i for i in range(1,35)]
var_list_omega = [i for i in range(ndata)]

model = Model()

U = model.binary_var_matrix(PDN_line_list, var_list_K, name = 'U')
Z = model.binary_var_list(var_list_location, name = 'Z')  #Whether to build a station
X = model.integer_var_list(var_list_location, name = 'X') #How many capacity

S = model.continuous_var_matrix(var_list_omega, var_list_Jplus, name = 'S')
Y = model.integer_var_cube(var_list_omega, var_list_TN,var_list_Jplus, name = 'Y')
P = model.continuous_var_matrix(var_list_omega,PDN_line_list, name = 'P')
Q = model.continuous_var_matrix(var_list_omega,PDN_line_list, name = 'Q')
V = model.continuous_var_matrix(var_list_omega,var_list_PDN_node, name = 'V')
R = model.continuous_var_cube(var_list_omega,PDN_line_list, var_list_K, name = 'R')
T = model.continuous_var_cube(var_list_omega, PDN_line_list, var_list_K, name = 'T')
#PI = model.continuous_var(name = 'PI')

obj_mp =alpha *( sum(A[i]*Z[i] for i in range(J)) +  sum(C[i] * X[i] for i in range(J)) ) \
    +  (sum( G_dict[i] * (sum(U[(i,k)] for k in range(K))) for i in PDN_line_list) )\
        + (1/ndata) * (beta * (H * sum(S[(omega,i)] for i in var_list_Jplus for omega in var_list_omega)) \
    - gamma * G * sum( Utility[i][ var_list_location.index(j) +1 ] * (Y[(omega,i,j)]) for i in var_list_TN for j in var_list_location for omega in var_list_omega))


model.minimize(obj_mp)

X_star = [0] + X
Z_star = [1] + Z

################## ################## ################## 
##################   TN    ####################
################## ################## ################## 

for j in range(J):
    model.add_constraint(X[j] <=  M * Z[j])

# for j in range(J):
#     model.add_constraint(S[j] <=  M * Z[j])
        
    #The chargring capacity cannot be negative

    
for omega in var_list_omega:    
    for j in var_list_Jplus:
            model.add_constraint( (sum(Y[(omega,i,j)] for i in var_list_TN) <= 80000000000* Z_star[var_list_Jplus.index(j)]))
            #model.add_constraint( S[var_list_Jplus.index(j)]  <=  8000000* Z_star[var_list_Jplus.index(j)]  )


for omega in var_list_omega:  
    for j in var_list_Jplus:
            model.add_constraint( (sum(Y[(omega,i,j)] for i in var_list_TN)- S[(omega,j)]) <= X_star[var_list_Jplus.index(j)])
            model.add_constraint( (sum(Y[(omega,i,j)] for i in var_list_TN)- S[(omega,j)]) >= 0 )
            
#(2)
#beta    
for omega in var_list_omega:
    for i in var_list_TN:
        for j in var_list_Jplus:
            for k in var_list_Jplus:
                model.add_constraint( (Utility[i][var_list_Jplus.index(j)] - max(Utility[i][var_list_Jplus.index(k)]-phi*Ud, Umin))*Y[(omega,i,j)] \
                           >= -800000000000*(1-Z_star[var_list_Jplus.index(k)])  )
                    #-800000000000*(1-Z_star[var_list_Jplus.index(k)])
for omega in var_list_omega:  
    for i in var_list_TN:
        model.add_constraint( sum(Y[(omega,i,j)] for j in var_list_Jplus) ==  D_mat[i][omega] ) 
    

################## ################## ################## 
##################   PDN    ####################
################## ################## ################## 


for omega in var_list_omega: 
    for n in var_list_PDN_node:
        if n in var_list_PDN:
            if n != 1:
                index = var_list_PDN.index(n)
                model.add_constraint(sum(P[(omega,line)] for line in PDN_line_list if line[1] == n) \
                                    -sum( P[(omega,line)] for line in PDN_line_list if line[0] == n)\
                                    - W[n-1] * sum( Y[(omega,i, var_list_location[index])] for i in var_list_TN) + W[n-1]*S[(omega,var_list_Jplus[index+1])] \
                                    == P_load_mat[n-1][0]  
                                )
        else:
            #index = var_list_PDN.index(n)
            if n != 1:
                model.add_constraint(sum(P[(omega,line)] for line in PDN_line_list if line[1] == n) \
                                    -sum( P[(omega,line)] for line in PDN_line_list if line[0] == n)\
                                    == P_load_mat[n-1][0]  
                                )
    
    for n in var_list_PDN_node:
        if n in var_list_PDN:
            if n != 1:
                index = var_list_PDN.index(n)
                model.add_constraint(sum(Q[(omega,line)] for line in PDN_line_list if line[1] == n) \
                                    -sum( Q[(omega,line)] for line in PDN_line_list if line[0] == n)\
                                    - 0.5 * sum( Y[(omega,i, var_list_location[index])] for i in var_list_TN) + 0.5 * S[(omega,var_list_Jplus[index+1])] \
                                    == Q_load_mat[n-1][0]  
                            )
        else:
            #index = var_list_PDN.index(n)
            if n != 1:
                model.add_constraint(sum(Q[(omega,line)] for line in PDN_line_list if line[1] == n) \
                                    -sum( Q[(omega,line)] for line in PDN_line_list if line[0] == n)\
                                    == Q_load_mat[n-1][0]  
                                )
    
    model.add_constraint(V[(omega,1)] == 1000000000)
    
    for line in PDN_line_list:
        model.add_constraint(V[(omega,line[0])] +\
                                          sum(R[(omega,line, k)] for k in (var_list_K)) - 
                                          V[(omega,line[1])]  - \
                                          sum(T[(omega,line, k)] for k in (var_list_K)) >= \
                                          2* R_dict[line] * P[(omega,line)] + 2*X_dict[line] * Q[(omega,line)])
        
    for n in var_list_PDN_node:
        model.add_constraint(V[(omega,n)] >= Vmin )
        model.add_constraint(V[(omega,n)] <= Vmax )
    
    for k in range(K):
        for line in PDN_line_list:
            model.add_constraint(-R[(omega,line,k)] >= -Vmax * U[(line, k)])
            model.add_constraint(R[(omega,line,k)] >= 0)
    
    
    for k in range(K):
        for line in PDN_line_list:
            model.add_constraint(-V[(omega,line[0])] + \
                                      R[(omega,line,k)] >= -Vmax*(1-U[(line,k)]))
            model.add_constraint(V[(omega,line[0])] - \
                                      R[(omega,line,k)] >= 0)
                
    
    for k in range(K):
        for line in PDN_line_list:
            model.add_constraint(-T[(omega,line,k)] >= -Vmax * U[(line, k)])
            model.add_constraint(T[(omega,line,k)] >= 0)
    
    
    for k in range(K):
        for line in PDN_line_list:
            model.add_constraint(-V[(omega,line[0])] + \
                                      T[(omega,line,k)] >= -Vmax*(1-U[(line,k)]))
            model.add_constraint(V[(omega,line[0])] - \
                                      T[(omega,line,k)] >= 0)
                
    
    for line in PDN_line_list:
        #For P
        model.add_constraint(P[(omega,line)] >= 0 )
        model.add_constraint(-P[(omega,line)] >= -Pmax[line] * (1+sum(U[(line, k)] for k in range(K)))  )
    
    
    for line in PDN_line_list:
        #For P
        model.add_constraint(Q[(omega,line)] >= 0 )
        model.add_constraint(-Q[(omega,line)] >= -Qmax[line] * (1+sum(U[(line, k)] for k in range(K)))  )
    
    
                    
    



'''''''''''''''''''''
Analysis of result
'''''''''''''''''''''''


#retract the values
sol = model.solve()
X_sol = sol.get_values(X)
X_star = [0] + X_sol 
print(X_sol)
Z_sol = sol.get_values(Z)
Z_star  =  [1] + Z_sol
print(Z_sol)
obj_value = sol.get_objective_value()
S_sol = sol.get_value_dict(S)
U_sol = sol.get_value_dict(U)
P_sol = sol.get_value_dict(P)
Y_sol = sol.get_value_dict(Y)
        
        

#V_sol = sol.get_values(V)
#print(V_sol)
#P_sol = sol.get_value_dict(P)
#print(S_sol[1:])
print("Objective Value: {}".format(obj_value))

Y_sol = sol.get_value_dict(Y)

penalty =  (1/ndata) * (beta * (H * sum(S_sol[(omega,i)] for i in var_list_Jplus for omega in var_list_omega)))
utility = (1/ndata)*gamma * G * sum( Utility[i][ var_list_location.index(j) +1 ] * (Y_sol[(omega,i,j)]) for i in var_list_TN for j in var_list_location for omega in var_list_omega)
# print(sum(S_sol))
# cost = sum(A[i]*Z_sol[i] for i in range(J)) + sum(C[i] * X_sol[i] for i in range(J))  + sum( G_dict[i] * (sum(U_sol[(i,k)] for k in range(K))) for i in PDN_line_list)
# total_utility = G * sum( Utility[i][ var_list_Jplus.index(j)] * Y_sol[(i,j)] for i in var_list_TN for j in var_list_Jplus ) 

Y_sol = sol.get_value_dict(Y)
#print(Y_sol)
# for i,v in Y_sol.items():
#     if v != 0:
#         print(i,v)
demand = {}
location_flow = {}
for i,v in Y_sol.items():
    if v != 0:
        #print(i,v)
        if i[0] == 2:
            if i[2] in location_flow:
                location_flow[i[2]] += v
            else:
                location_flow[i[2]] = v
            
        demand[i] = v
        csv_file = 'Demand.csv'
        with open('dict.csv', 'w', newline="") as csv_file:  
            writer = csv.writer(csv_file)
            for key, value in demand.items():
                writer.writerow([key, value])
                
for i,v in sorted(location_flow.items()):
    print(i,v)
            


cost = alpha *( sum(A[i]*Z_sol[i] for i in range(J)) +  sum(C[i] * X_sol[i] for i in range(J)) ) + \
    (sum( G_dict[i] * (sum(U_sol[(i,k)] for k in range(K))) for i in PDN_line_list) )  


# print("Delta U:", Ud)
# print("Utility:", gamma *total_utility-beta * penalty)
# print("Building Cost: ",alpha*cost)
print("Penalty Cost: ", penalty)
print("Utility:", utility)
print("Building Cost: ",cost)




                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
