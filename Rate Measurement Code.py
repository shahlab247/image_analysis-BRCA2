# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:52:12 2022

@author: ritik
#### Rate Calculator 
"""
import numpy as np
import pandas as pd
#importing sheet
import matplotlib.pyplot as plt

#%%
#Inputting file with autophagasome, autolysosme, cell count and well tags indicating well number and position 
df= pd.read_excel(r")
dfn = df.to_numpy()

timepoints = len(df['ND.T'].value_counts()) #counting the number of unique timepoints 
positions = 3 #INPUT NUMBER OF POSITIONS IMAGED FOR EACH WELL
conds_unique = len(df['ND.M'].value_counts())/positions #calculates number of conditions tested 
baf = timepoints- 3 #indicates the timepoint in the timepoint vector where bafilomycin was added 


#%%
#Averaging out autophagasomes and autolysosomes per cell after removing outliers (distured positions)
AverageAP = dfn[:,8]
AverageAL = dfn[:,9]

avg_values= np.zeros((int((timepoints*conds_unique)),2))
for unique in range(int(conds_unique)):
    for time in range(int(timepoints)):
        idx=int(unique*positions*timepoints + time)
        idx_cal=unique*timepoints+ time
        print("idx:",idx) 
        print("idx_cal:",idx_cal)
        avg_values[idx_cal,0]=np.mean([AverageAP[idx],AverageAP[idx+timepoints],AverageAP[idx+2*timepoints],AverageAP[idx+3*timepoints]]) #averaged (over positions) autophagosomes/cell 
        avg_values[idx_cal,1]=np.mean([AverageAL[idx],AverageAL[idx+timepoints],AverageAL[idx+2*timepoints],AverageAL[idx+3*timepoints]]) #averaged (over positions) autolysosomes/cell 
#%%
#Calculating Rates 
rates = np.zeros((int(conds_unique),5))
for unique in range(int(conds_unique)):
    for time in range(timepoints):
        idx = unique*timepoints+ baf
        rates[unique,0] = (avg_values[idx+2,0]- avg_values[idx,0])/20 #R1
        rates[unique,1] = (avg_values[idx,0]- avg_values[idx-2,0])/20 #dAP/dt_0
        rates[unique,2] = np.subtract(rates[unique,0],rates[unique,1]) #R2 
        rates[unique,3] = (avg_values[idx,1]- avg_values[idx-2,1])/20 #dAL/dt_0
        rates[unique,4] = np.subtract(rates[unique,2],rates[unique,3]) #R3
        
rate_df=pd.DataFrame(rates,columns=['R1','dAP/dt_0','R2','dAL/dt_0','R3'])
rate_df.to_excel(r"G:\.shortcut-targets-by-id\1m1bQR-tnQCfKRd3zgPEENJc3VpOD4fdg\Shared-raw_dataand_experimental_setup\Pancreatic Cancer Cells\Experimental Data\10272022-JQ1-Rap-Rates\PCC20X files\Rates-KN-12hrv3.xlsx")

#%%
#Bar Graph with Rates 

RateGraph = plt.figure()
ax = RateGraph.add_axes([0,0,1,1])
X = np.arange(4)
ax.bar(X + 0.00, rate_df['R1'], color = 'darkorange', width = 0.25)
ax.bar(X + 0.25, rate_df['R2'], color = 'peru', width = 0.25)
ax.bar(X + 0.50, rate_df['R3'], color = 'peachpuff', width = 0.25)
ax.legend(labels=['R1', 'R2', 'R3'], fontsize =15)
ax.set_xticks(X, labels=['JQ1- BRCA2 KO', 'DMSO- BRCA2 KO ', 'JQ1 -Control', 'DMSO- Control' ]) 
ax.set_ylabel('Puncta per cell per minute', fontsize = 15)
ax.set_ylim(-0.2,4)


