import numpy as np
import time, os, sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import fnmatch
import pandas as pd
from scipy import stats
plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600

#%%
#Figure-A

df= pd.read_excel(r"G:\My Drive\Research\nikonTi2\systems autophagy figures\Pancreatic cancer\Data.xlsx",sheet_name='Figure-A')
JQ1_KO= df.iloc[0,1:7].to_numpy()
DMSO_KO=df.iloc[1,1:7].to_numpy()

## plotting figures
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

times= ['R1'.translate(SUB), 'R2'.translate(SUB), 'R3'.translate(SUB)]
conditions= ['JQ1-BRAC2 KO',]
x= np.arange(len(times))/3.5
width= 0.05

c1= [(7,152,171), (251,119,18),(150,211,219),(253, 196, 153)] #enter your color values. 
c2=np.divide(c1,255)

fig, ax = plt.subplots(figsize=(6,5))
rects1 = ax.bar(x - width/2, JQ1_KO[0:3], width, yerr=JQ1_KO[3:], capsize=4, color= c2[2], edgecolor='k',label='BRCA2 KO',zorder=-1)
rects2 = ax.bar(x + 1*width/2, DMSO_KO[0:3], width, yerr=DMSO_KO[3:], capsize=4,color= c2[3],edgecolor='k',label='Control', zorder=-1)

for i in range(len(times)):
    reps=3
    ref=7+i*reps
    ax.scatter(x[i]-width/2 + 0.8*np.random.random(df.iloc[0][ref:ref+reps].size) * width/4, df.iloc[0][ref:ref+reps], color=c2[2],edgecolor='k', zorder=1)
    ax.scatter(x[i]+ 1*width/2 + 0.8*np.random.random(df.iloc[1][ref:ref+reps].size) * width/4, df.iloc[1][ref:ref+reps], color=c2[3],edgecolor='k',zorder=1)

ax.set_ylabel('Rates (Puncta/Cell/Min)')

ax.set_xticks(x)
ax.set_xticklabels(times)
ax.tick_params(axis='y',direction="in", which='both')
ax.autoscale(False)
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
plt.yticks(np.arange(0, 3.6, 0.2)) 
plt.show()

#%%
#Figure-B
df= pd.read_excel(r"G:\My Drive\Research\nikonTi2\systems autophagy figures\Pancreatic cancer\Data.xlsx",sheet_name='Figure-B')

JQ1_KO= df.iloc[0,1:7].to_numpy()
DMSO_KO=df.iloc[1,1:7].to_numpy()
JQ1_BS=df.iloc[2,1:7].to_numpy()
DMSO_BS=df.iloc[3,1:7].to_numpy()


## plotting figures
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

times= ['R1'.translate(SUB), 'R2'.translate(SUB), 'R3'.translate(SUB)]
conditions= ['JQ1-BRAC2 KO','DMSO KO','JQ-1 BS','DMSO-BS']
x= np.arange(len(times))/3.5
width= 0.05

c1= [(7,152,171), (251,119,18),(150,211,219),(253, 196, 153)] #enter your color values. 
c2=np.divide(c1,255)

fig, ax = plt.subplots(figsize=(6,5))
rects1 = ax.bar(x - width/2, JQ1_KO[0:3], width, yerr=JQ1_KO[3:], capsize=4, color= c2[0], edgecolor='k',label='JQ1-BRCA2 KO',zorder=-1)
rects2 = ax.bar(x + 1*width/2, DMSO_KO[0:3], width, yerr=DMSO_KO[3:], capsize=4,color= c2[2],edgecolor='k',label='DMSO-BRCA2 KO', zorder=-1)
rects3 = ax.bar(x + 3*width/2, JQ1_BS[0:3], width, yerr=JQ1_BS[3:], capsize=4,color= c2[1],edgecolor='k',label='JQ1-Control', zorder=-1)
rects4 = ax.bar(x + 5*width/2, DMSO_BS[0:3], width, yerr=DMSO_BS[3:], capsize=4,color= c2[3],edgecolor='k',label='DMSO-Control', zorder=-1)

for i in range(len(times)):
    reps=3
    ref=7+i*reps
    ax.scatter(x[i]-width/2 + 0.8*np.random.random(df.iloc[0][ref:ref+reps].size) * width/4, df.iloc[0][ref:ref+reps], color=c2[0],edgecolor='k', zorder=1)
    ax.scatter(x[i]+ 1*width/2 + 0.8*np.random.random(df.iloc[1][ref:ref+reps].size) * width/4, df.iloc[1][ref:ref+reps], color=c2[2],edgecolor='k',zorder=1)
    ax.scatter(x[i]+ 3*width/2 + 0.8*np.random.random(df.iloc[2][ref:ref+reps].size) * width/4, df.iloc[2][ref:ref+reps], color=c2[1],edgecolor='k',zorder=1)
   
    ax.scatter(x[i]+ 5*width/2 + 0.8*np.random.random(df.iloc[3][ref:ref+reps].size) * width/4, df.iloc[3][ref:ref+reps], color=c2[3],edgecolor='k',zorder=1)
   
ax.set_ylabel('Rates (Puncta/Cell/Min)')
ax.set_xticks(x)
ax.set_xticklabels(times)

ax.tick_params(axis='y',direction="in", which='both')
ax.tick_params(axis='x',direction="in", which='both')
ax.autoscale(False)
ax.legend(loc='upper right', bbox_to_anchor=(1.4, 1))
plt.yticks(np.arange(0, 3.0, 0.2)) 

plt.show()

#%%
#supplemental data
df= pd.read_excel(r"G:\My Drive\Research\nikonTi2\systems autophagy figures\Pancreatic cancer\Data.xlsx",sheet_name='SI-Rapa')


JQ1_BS= df.iloc[0,1:7].to_numpy()
DMSO_BS=df.iloc[1,1:7].to_numpy()

## plotting figures
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
times= ['R1'.translate(SUB), 'R2'.translate(SUB), 'R3'.translate(SUB)]
conditions= ['JQ1-BRAC2 KO','DMSO KO','JQ-1 BS','DMSO-BS']
x= np.arange(len(times))/3.5
width= 0.05

c1= [(7,152,171), (251,119,18),(150,211,219),(253, 196, 153)] #enter your color values. 
c2=np.divide(c1,255)

fig, ax = plt.subplots(figsize=(6,5))

rects1 = ax.bar(x + 1*width/2, JQ1_BS[0:3], width, yerr=JQ1_BS[3:], capsize=4,color= c2[1],edgecolor='k',label='Rapa-Control', zorder=-1)
rects2 = ax.bar(x + 3*width/2, DMSO_BS[0:3], width, yerr=DMSO_BS[3:], capsize=4,color= c2[3],edgecolor='k',label='DMSO-Control', zorder=-1)

for i in range(len(times)):
    ref=7+i*3
    ax.scatter(x[i]+ 1*width/2 + 0.8*np.random.random(df.iloc[0][ref:ref+3].size) * width/4, df.iloc[0][ref:ref+3], color=c2[1],edgecolor='k',zorder=1)
    ax.scatter(x[i]+ 3*width/2 + 0.8*np.random.random(df.iloc[1][ref:ref+3].size) * width/4, df.iloc[1][ref:ref+3], color=c2[3],edgecolor='k',zorder=1)
   
ax.set_ylabel('Rates (Puncta/Cell/Min)')
ax.set_xticks(x)
ax.set_xticklabels(times)
ax.tick_params(axis='y',direction="in", which='both')
ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1))
plt.yticks(np.arange(0, 2.2, 0.2)) 



