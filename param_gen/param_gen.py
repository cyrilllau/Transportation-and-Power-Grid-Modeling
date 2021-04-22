# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 23:04:40 2020

@author: lusha
"""

import pandas as pd, numpy as np, matplotlib.pyplot as plt

if __name__ == '__main__':  
    
    '''generating PDN load'''
    cols = ['node','PL','QL']
    load_df = pd.read_csv('DNdata_org.txt', sep="\t", names = cols, skiprows = 1, header=None)
    print(load_df)
    
    resi_nodes = [2,3,23,24,25,29,30,31,32,33,18,17,16,13,12,22,21,20,19]
    comm_nodes = [4,5,6,7,8,9,10,11,14,15,26,27,28]
    
    cols = ['resi_power']
    resi_df = pd.read_csv('DN_resid_profile.txt', sep="\t", names = cols, header=None)
    print(resi_df)
    resi_base = float(resi_df['resi_power'][0])
    for i in range(len(resi_df)):
        resi_df['resi_power'][i] = float(resi_df['resi_power'][i])/resi_base
    print(resi_df)
    resi_df.plot()
    
    cols =['comm_power']
    comms_df = pd.read_csv('DN_commer_profile.txt', sep="\t", names = cols, header=None)
    print(comms_df)
    comm_base = float(comms_df['comm_power'][0])
    for i in range(len(comms_df)):
        comms_df['comm_power'][i] = float(comms_df['comm_power'][i])/comm_base
    print(comms_df)
    comms_df.plot()
    
    OutputFile = open("output_param.txt", "w+")
    for t in range(24):
        OutputFile.write("t = "+str(t)+"\n")
        for node in resi_nodes:
            mu_p = load_df['PL'][node-1] * resi_df['resi_power'][t]
            sigma = 10
            s = np.random.normal(mu_p, sigma, 20)
            OutputFile.write("Node " + str(node) + " PL scenarios: " + str(s) + "\n")
            mu_q = load_df['QL'][node-1] * resi_df['resi_power'][t]
            sigma = 10
            s = np.random.normal(mu_q, sigma, 20)
            OutputFile.write("Node " + str(node) + " QL scenarios: " + str(s) + "\n")
            
        for node in comm_nodes:
            mu_p = load_df['PL'][node-1] * comms_df['comm_power'][t]
            sigma = 10
            s = np.random.normal(mu_p, sigma, 20)
            OutputFile.write("Node " + str(node) + " PL scenarios: " + str(s) + "\n")
            mu_q = load_df['QL'][node-1] * comms_df['comm_power'][t]
            sigma = 10
            s = np.random.normal(mu_q, sigma, 20)
            OutputFile.write("Node " + str(node) + " QL scenarios: " + str(s) + "\n")
        OutputFile.write("\n")
    
    '''generating TN drive time'''
    resi_nodes = [1,2,3,6,7,8,12,13,16,17,18,19,20,21,24]
    comm_nodes = [4,5,9,10,11,14,15,22,23]
    
    travel_time_base = {}
    origin = 0
    TNfile = open('SiouxFalls_trips.tntp.txt','r')
    lines = TNfile.readlines()
    for line in lines:
        if line.startswith("Origin"):
            origin = int(line.split()[1])
            continue
        if origin:
            des_pairs = line.split(";")
            for des_pair in des_pairs:
                if not des_pair.isspace():
                    des = int(des_pair.split(":")[0])
                    t = float(des_pair.split(":")[1])
                    travel_time_base[(origin, des)] = t
    print(travel_time_base)
    
    cols = ['time','traffic time']
    resi_comm_df = pd.read_csv('resi_comm_traffic_profile.txt', sep="\t", names = cols, header=None)
    resi_comm_df['traffic time'] = resi_comm_df['traffic time'].astype(float)
    for i in range(len(resi_comm_df)):
        resi_comm_df['traffic time'][i] = resi_comm_df['traffic time'][i]/100
    print(resi_comm_df)
    resi_comm_df.plot(x='time',y='traffic time')
    
    
    cols = ['time','traffic time']
    comm_res_df = pd.read_csv('comm_resi_traffic_profile.txt', sep="\t", names = cols, header=None)
    comm_res_df['traffic time'] = comm_res_df['traffic time'].astype(float)
    for i in range(len(comm_res_df)):
        comm_res_df['traffic time'][i] = comm_res_df['traffic time'][i]/100
    print(comm_res_df)
    comm_res_df.plot(x='time',y='traffic time')
    
    
    for t in range(24):
        OutputFile.write("t = "+str(t)+"\n")
        for origin in resi_nodes:
            for des in resi_nodes:
                mu_t =  travel_time_base[(origin, des)] 
                if mu_t == 0:
                    s = [0]*20
                else:
                    sigma = 10
                    s = np.random.normal(mu_t, sigma, 20)
                OutputFile.write("Origin " + str(origin) + " Des: " + str(des) + " traffic time: " + str(s) + "\n")
        
        for origin in resi_nodes:
            for des in comm_nodes:
                mu_t = travel_time_base[(origin, des)]  * resi_comm_df["traffic time"][t]
                if mu_t == 0:
                    s = [0]*20
                else:
                    sigma = 10
                    s = np.random.normal(mu_t, sigma, 20)
                OutputFile.write("Origin " + str(origin) + " Des: " + str(des) + " traffic time: " + str(s) + "\n")

        for origin in comm_nodes:
            for des in comm_nodes:
                mu_t =  travel_time_base[(origin, des)] *2
                if mu_t == 0:
                    s = [0]*20
                else:
                    sigma = 10
                    s = np.random.normal(mu_t, sigma, 20)
                OutputFile.write("Origin " + str(origin) + " Des: " + str(des) + " traffic time: " + str(s) + "\n")
        
        for origin in comm_nodes:
            for des in resi_nodes:
                mu_t = travel_time_base[(origin, des)]  * comm_res_df["traffic time"][t]
                if mu_t == 0:
                    s = [0]*20
                else:
                    sigma = 10
                    s = np.random.normal(mu_t, sigma, 20)
                OutputFile.write("Origin " + str(origin) + " Des: " + str(des) + " traffic time: " + str(s) + "\n")

        
        OutputFile.write("\n")
    OutputFile.close()