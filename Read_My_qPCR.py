# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 12:53:25 2021

@author: Alison
"""

#####################################################


""" Start by changing the following parameters """

workingdir = r"C:\Users\Darwin\Documents\Alison\AL Data\B2P43_Liposomes_345_Fractions\qPCR Liposomes Overnight"
fname = r"C:\Users\Darwin\Documents\Alison\AL Data\B2P43_Liposomes_345_Fractions\qPCR Liposomes Overnight\2021-09-01 Liposomes Overnight -  Quantification Amplification Results_FAM.csv"
AverageDatainTriplicates = True

fname_h = r"C:\Users\Darwin\Documents\Alison\AL Data\B2P43_Liposomes_345_Fractions\qPCR Liposomes Overnight\2021-09-01 Liposomes Overnight -  Headers.csv"

t_per_run = 10.133333333 # minutes

#####################################################


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir(workingdir)

# extract emission data
df = pd.read_csv(fname)
data = df.to_numpy()

# extract header data
dfh = pd.read_csv(fname_h)
headers = dfh.to_numpy()
headers = headers[0]

# X data (cycle no.) sits in 2nd column
cycle = data[:,1] - 1 # minus one because cycle starts at 1 whereas time starts at 0
Time = cycle * t_per_run

# FAM Emission data sits in 3rd and later columns
data = data[:,2:]
H,N = np.shape(data)

# average in triplicates
t = int(N/3)
# data goes here
em_avgs = []
em_std = []
header_avgs = []

for x in range(0,t):
    r = int(x*3)
    em_avgs = np.append(em_avgs, np.mean(data[:,r:r+3],1)).reshape((x+1,H))
    em_std = np.append(em_std, np.std(data[:,r:r+3],1)).reshape((x+1,H))
    header_avgs = np.append(header_avgs, headers[r])

#####################################################

# plot averaged data with error bars

fig1, ax1 = plt.subplots(1, 1)  
  
for x in range(0,t):
    evenly_spaced_interval = np.linspace(0, 1, t)
    colors = [plt.cm.viridis(x) for x in evenly_spaced_interval]
    
    ax1.errorbar(Time, em_avgs[x,:], em_std[x,:], color=colors[x], label = header_avgs[x])
    
    print(x)
    
#plt.ylim((7000, 8000))

ax1.set_ylabel('FAM Emission')
ax1.set_xlabel('Time (minutes)')
ax1.set_title("qPCR CF-LUV's")
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()
#    fig.savefig("Cumulative" + ender + ".png")


#####################################################

# plot each of the 5 sample types in a stacked subplot

fig, axs = plt.subplots(5)
fig.set_size_inches(8, 12)

tt = int(t/3)

for x in range(0,tt):
    
    evenly_spaced_interval = np.linspace(0, 1, tt)
    colors = [plt.cm.viridis(x) for x in evenly_spaced_interval]
    
    for i in range(0,3):
        
        
        dstyle = ['-', '--', ':'] # line styles: solid, dashed, dotted
        leg = ["0.0%", "0.1%", "0.2%"]        
        index = x*3 + i
        print("index =" + str(index))
        
        axs[x].errorbar(Time, em_avgs[index,:], em_std[index,:], color=colors[x], linestyle=dstyle[i], label = leg[i])
        axs[x].set_ylabel('FAM Emission')
        axs[x].set_xlabel('Time (minutes)')
        axs[x].set_xlim(0,250)
        
        temptitle = header_avgs[index]
        subplot_title = temptitle[:-7]
        
        axs[x].set_title(subplot_title)
        axs[x].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.tight_layout()


#####################################################

# normalize each set of frac 4, 5, 6 data to the 0.2% endpoint at cycle 70

# find the correct endpoints
e_index = [8,11,14]
endpoints = []
headers_norm = []
for i in e_index:
    endpoints = np.append(endpoints,em_avgs[i,-1])
    headers_norm = np.append(headers_norm,headers[i])

temp = endpoints
endpoints = []
for i in temp:
    endpoints = np.append(endpoints,[i,i,i])

# section out just the emission data (not the control data sets)
fraction_avgs = em_avgs[-9:]
fraction_std = em_std[-9:]
fraction_headers = header_avgs[-9:]

# divide emission data by the endpoints
em_norm = fraction_avgs / endpoints[:,None]
em_norm_std = fraction_std / endpoints[:,None]

#####################################################

# plot each of the 3 now normalized samples in a stacked subplot
fig, axs = plt.subplots(3)
fig.set_size_inches(8, 7)

L = 3

for x in range(0,L):
    
    for i in range(0,3):
        
        c2 = x+2 # to keep the fraction colors and subplot titles the same, use c2 instead of i
        
        dstyle = ['-', '--', ':'] # line styles: solid, dashed, dotted
        leg = ["0.0%", "0.1%", "0.2%"]        
        index = x*3 + i
        print("index =" + str(index))
        
        axs[x].errorbar(Time, em_norm[index,:], em_norm_std[index,:], color=colors[c2], linestyle=dstyle[i], label = leg[i])
        axs[x].set_ylabel('FAM Emission')
        axs[x].set_xlabel('Time (minutes)')
        axs[x].set_xlim(0,250)
        
        temptitle = header_avgs[e_index[x]]
        subplot_title = temptitle[:-7] + " - Normalized"
        
        axs[x].set_title(subplot_title)
        axs[x].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.tight_layout()


#####################################################

# Plot normalized fractions on one graph

# plot each of the 3 now normalized samples in a stacked subplot
fig, axs = plt.subplots(1)
fig.set_size_inches(8, 4)

L = 3

for x in range(0,L):
    
    for i in range(0,3):
        
        c2 = x+2 # to keep the fraction colors and subplot titles the same, use c2 instead of i
        
        dstyle = ['-', '--', ':'] # line styles: solid, dashed, dotted
        leg = ["0.0%", "0.1%", "0.2%"]        
        index = x*3 + i
        print("index =" + str(index))
        
        axs.errorbar(Time, em_norm[index,:], em_norm_std[index,:], color=colors[c2], linestyle=dstyle[i], label = header_avgs[i])
        axs.set_ylabel('FAM Emission')
        axs.set_xlabel('Time (minutes)')
        axs.set_xlim(0,250)
        
        
axs.set_title('Fractions 4, 5, 6 - Normalized')        
axs.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()



