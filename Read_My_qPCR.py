# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 12:53:25 2021

@author: Alison
"""

#####################################################

"""
To Do:

[X] - Remove unecessary printing of data points in loops
[]  - Add index to data array so that I can handle multiple fluorophores
    [X] - completed for Plot 1
    []  - completed for Plot 2
    []  - completed for Plot 3
    []  - completed for Plot 4
    []  - completed for Plot 5
[]  - Add bar graph for fractions averaged over time
[]  - Add Y/N options for each type of graph at the beginning parameters
[]  - Update "Headers" file input to match a standard 96 well plate

"""


""" Start by changing the following parameters """

wdir    = r"C:\Users\sunsh\Documents\AL Data\B2P51\qPCR"
fname   = r"C:\Users\sunsh\Documents\AL Data\B2P51\qPCR\20210904_CF-LUV_Overnight_AuNP_last6columns -  Quantification Amplification Results_FAM.csv"
fname_h = r"C:\Users\sunsh\Documents\AL Data\B2P51\qPCR\20210904_CF-LUV_Overnight_AuNP_last6columns -  Headers.csv"

AverageDatainTriplicates = True


t_per_run = 10.133333333 # minutes

fluor_name = ['FAM', 'Texas Red', 'Cal Gold 540']

fluor_FAM = True
fluor_TexasRed = True
fluor_CalGold = True

#####################################################

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

# Colors

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

c_FAM = '#2E8B57' #seafoam green
c_TxR = '#A40000' #dark red
c_Cal = '#EDC812' #gold
c_wht = '#FFFFFF' #white

fluor_colors = [c_FAM, c_TxR, c_Cal]

# Get Data

os.chdir(wdir)

# extract emission data
base = fname[:-7]
df = pd.read_csv(fname)
data = df.to_numpy()
H,N = np.shape(data)

fluor = [fluor_FAM, fluor_TexasRed, fluor_CalGold]
N_fluor = sum(1 for x in fluor if x == True)
# remove any unused fluorophores in name list
fluor_index = []
for x in range(0,len(fluor)):
    if fluor[x] == True:
        fluor_index = np.append(fluor_index, x)

# put fluorophore data into empty emission data array
data = np.empty((H,N,N_fluor))
for x in range(0,N_fluor):
    i = int(fluor_index[x])
    temp = pd.read_csv(base + fluor_name[i] + ".csv")
    temp = temp.to_numpy()
    data[:,:,x] = temp

# extract header data
dfh = pd.read_csv(fname_h)
headers = dfh.to_numpy()
headers = headers[0]

# X data (cycle no.) sits in 2nd column
cycle = data[:,1,:] - 1 # minus one because cycle starts at 1 whereas time starts at 0
cycle = cycle[:,0]
Time = cycle * t_per_run

# Emission data sits in 3rd and later columns
data = data[:,2:,:]
H,N,N_fluor = np.shape(data)

if AverageDatainTriplicates == True:
    # average in triplicates
    t = int(N/3)
    # data goes here
    temp_avg = np.empty((H,t,N_fluor))
    temp_std = np.empty((H,t,N_fluor))
    temp_hdr = [] # append strings to an empty array rather than pre-generate an empty array
 
    for i in fluor_index:
        i = int(i)
        for x in range(0,t):
            r = int(x*3)
            temp_avg[:,x,i] = np.mean(data[:,r:r+3,i],1)
            temp_std[:,x,i] = np.std( data[:,r:r+3,i],1)
            temp_hdr = np.append(temp_hdr, headers[r]) # only take the header from the first of the triplicates
    
    # rename variables 
    data = temp_avg
    stdev = temp_std
    headers = temp_hdr
    
#####################################################

# PLOT 1

# plot all data (with error bars if averaged)

for i in range(0,len(fluor)):
    if fluor[i] == True: # for each fluorophore, make a different plot
    
        # start figure for this fluorophore
        fig, ax = plt.subplots(1, 1)
        
        # define color gradient. Every individual experiment gets a color in the gradient
        colors = [colorFader(fluor_colors[i], c_wht, x/t) for x in range(0,t)]
    
        for x in range(0,t):
            if AverageDatainTriplicates == True:
                ax.errorbar(Time, data[:,x,i], stdev[:,x,i], color=colors[x], label=headers[x])
            else:
                ax.plot(Time, data[:,x,i], color=colors[x], label=headers[x])
        
        # set other parameters
        ax.set_ylabel('RFU Emission')
        ax.set_xlabel('Time (minutes)')
        ax.set_title(fluor_name[i])
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        plt.tight_layout()


#####################################################

# PLOT 2

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

# PLOT 3

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

# PLOT 4

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

# PLOT 5

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
        
        axs.errorbar(Time, em_norm[index,:], em_norm_std[index,:], color=colors[c2], linestyle=dstyle[i], label = header_avgs[i])
        axs.set_ylabel('FAM Emission')
        axs.set_xlabel('Time (minutes)')
        axs.set_xlim(0,250)
        
        
axs.set_title('Fractions 4, 5, 6 - Normalized')        
axs.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.tight_layout()



