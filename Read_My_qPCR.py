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
    [X] - completed for Plot 2
    []  - completed for Plot 3
    []  - completed for Plot 4
    []  - completed for Plot 5
[X] - Add bar graph for fractions averaged over time
[]  - Add Y/N options for each type of graph at the beginning parameters
[]  - Update "Headers" file input to match a standard 96 well plate
[]  - Extract Headers and legends directly from headers file (needs a consistent delimiter)

"""


""" Start by changing the following parameters """

wdir    = r"C:\Users\Alison\Documents\AL Data\B2P51\qPCR - 20210910"
fname   = r"C:\Users\Alison\Documents\AL Data\B2P51\qPCR - 20210910\20210910-CF-LUV_AuNP_GT15SWNT -  Quantification Amplification Results_FAM.csv"

fname_h = r"C:\Users\Alison\Documents\AL Data\B2P51\qPCR - 20210910\20210910-CF-LUV_AuNP_GT15SWNT -  Headers.csv"

SampleHeaders = ['20mM HEPES', '2.5uM CF', "CF-LUV's"]


AverageDatainTriplicates = True

N_SampleTypes = 3 # number of different rows (20mM HEPES, 2.5 uM CF, Fraction 1, etc.)

t_per_run = 10.133333333 # minutes

fluor_name = ['FAM', 'Texas Red', 'Cal Gold 540']

fluor_FAM = True
fluor_TexasRed = False
fluor_CalGold = False

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
c_blk = '#000000' #black

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
t = N

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
        colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, t)]
        #colors = [colorFader(fluor_colors[i], c_wht, x/t) for x in range(0,t)]
    
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

# plot each of the 3 sample types in a stacked subplot

# i is figure index (one figure for each fluorophore)
# n is subfigure index (one for each type of sample (HEPES control, CF control, each fraction sample, etc.))
# x is the sub-subfigure index. Within each sample type, how many variations are there (HEPES, 0.1% Triton, 2uL AuNP, 4 uL AuNP, etc.)

dstyle = ['-', '--', ':', '-.'] # line styles: solid, dashed, dotted, dash-dotted
leg = ["Control", "0.5% Triton", "2 uL AuNP", "4 uL AuNP"]       
X = int(t/N_SampleTypes)

for i in range(0,len(fluor)): # for each figure
    if fluor[i] == True: # for each fluorophore, make a different plot

        # start figure for this fluorophore
        fig, axs = plt.subplots(N_SampleTypes, 1)
        fig.set_size_inches(8, 4*N_SampleTypes)

        # define color gradient. Every individual experiment gets a color in the gradient
        #colors = [colorFader(fluor_colors[i], c_blk, x/t) for x in range(0,X)]
        colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, X)]

        for n in range(0,N_SampleTypes): # for each subfigure
        
            print("n = ",str(n))
            
            for x in range(0,X): # for each line in each subfigure
            
                print("x = ",str(x))
            
                t_index = int(n*X + x)     
                
                print("t_index = ",str(t_index))
                
                if AverageDatainTriplicates == True:                    
                    axs[n].errorbar(Time, data[:,t_index,i], stdev[:,t_index,i], color=colors[x], linestyle=dstyle[t_index], label = leg[x])
                else:
                    axs[n].plot(Time, data[:,t_index,i], color=colors[x], linestyle=dstyle[x], label = leg[x])
                
            axs[n].set_ylabel('RFU')
            axs[n].set_xlabel('Time (minutes)')
            
            axs[n].set_title(SampleHeaders[n])
            axs[n].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            
        fig.suptitle(fluor_name[int(fluor_index[i])], fontsize = 36)
        plt.tight_layout()


#####################################################

# PLOT 3

# normalize each set of fraction data to the Triton average at last timepoint

# find the correct endpoints
if AverageDatainTriplicates == True:
    e_index = [(i+2) * t for i in range(0,N_SampleTypes - 2)]
else:
    e_index = [i * X + 2 * X for i in range(0,N_SampleTypes - 2)]

endpoints = np.empty((len(e_index), len(fluor)))
headers_norm = np.empty(len(e_index))

for i in range(0,len(e_index)):
    endpoints[i,:] = np.append(endpoints,data[i,-1,:])
    headers_norm = np.append(headers_norm,headers[i])

# now turn the endpoints list match the size and shape of the data list (X*sample number)
temp = endpoints
endpoints = np.empty((X * N_SampleTypes, len(fluor)))

for i in temp: # for each fluorophore
    
    r = []
    
    for x in range(0,X):
        r = np.append(r,temp[i])

    endpoints = np.append(endpoints,r)

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



###################################################

# PLOT 6

# Bar graph!

# average data along all time points

bar_data = np.average(data, axis=0)

fig, ax = plt.subplots(1)
fig.set_size_inches(8, 4)

X = np.arange(t)

for i in fluor_index:
    i = int(i)
    ax.bar(X + i/4, bar_data[:,i], color = fluor_colors[i], width = 0.25, label = fluor_name[i])
           
plt.xticks(X, headers, rotation=90)
ax.legend()