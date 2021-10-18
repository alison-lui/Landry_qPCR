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

#wdir = r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P92\qPCR"
#fname   = r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P92\qPCR\2021-10-12_CF-LUV_overnight -  Quantification Amplification Results_FAM.csv"
#fname_h = r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P92\qPCR\2021-10-12_CF-LUV_overnight -  Headers.csv"

wdir = r"/Volumes/GoogleDrive/My Drive/Research/Landry Lab Summer Research 2021/AL Data/B2P92/qPCR"
fname = r"/Volumes/GoogleDrive/My Drive/Research/Landry Lab Summer Research 2021/AL Data/B2P92/qPCR/2021-10-12_CF-LUV_overnight -  Quantification Amplification Results_FAM.csv"
fname_h = r"/Volumes/GoogleDrive/My Drive/Research/Landry Lab Summer Research 2021/AL Data/B2P92/qPCR/2021-10-12_CF-LUV_overnight -  Headers.csv"

wdir = r"/Volumes/GoogleDrive/My Drive/Research/Landry Lab Summer Research 2021/AL Data/B2P76/qPCR"
fname = r"/Volumes/GoogleDrive/My Drive/Research/Landry Lab Summer Research 2021/AL Data/B2P76/qPCR/2012-09-23 CF-LUV Overnight -  Quantification Amplification Results_FAM.csv"
fname_h = r"/Volumes/GoogleDrive/My Drive/Research/Landry Lab Summer Research 2021/AL Data/B2P76/qPCR/2012-09-23 CF-LUV Overnight -  Headers.csv"

AverageDatainTriplicates = True

t_per_run = 2 # minutes

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

#####################################################

# Colors

c_FAM = '#2E8B57' #seafoam green
c_TxR = '#A40000' #dark red
c_Cal = '#EDC812' #gold
c_wht = '#FFFFFF' #white
c_blk = '#000000' #black

fluor_colors = [c_FAM, c_TxR, c_Cal]

#####################################################


            
            

# Get Data

os.chdir(wdir)

# extract emission data
base = fname[:-7]
df = pd.read_csv(fname)
data = df.to_numpy()
H,N = np.shape(data)

# list fluorophores
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

    for x in range(0,t):
        r = int(x*3)
        temp_hdr = np.append(temp_hdr, headers[r]) # only take the header from the first of the triplicates

    # rename variables
    data = temp_avg
    stdev = temp_std
    headers = temp_hdr
    
        
# Header Formatting

# split each header by delimeter '-'
headers = [i.split(" - ") for i in headers]
header_index = np.arange(0,len(headers))

# only blank headers don't have additives and are length 1. Add a 'blank' tag to these to make them all the same length
for i in np.arange(0,len(headers)):
    if len(headers[i]) == 1:
        headers[i] = [headers[i][0], 'blank']

# bring together all split headers and add index number into one array like: [main, additive, index]
temp = []
for i in np.arange(0,len(headers)):
    temp = np.append(temp, headers[i])
    temp = np.append(temp, header_index[i])
temp = np.reshape(temp, (len(header_index),3))

# remove any headers that say 'blank'
#for i in reversed(np.arange(0,len(headers))):
#    if temp[i,1] == 'blank':
#        temp = np.delete(temp, i, 0)

# place new data into headers with blanks removed
headers = temp

# identify unique mains and additives in list
headers_main = list(set(headers[:,0]))
headers_adds = list(set(headers[:,1]))

#####################################################

# Function for making plots

def makeplot(data, stdev, main, adds, fi, FIG, AX, colors = 0, labels='adds', leg=True):

    """
    This function takes the 'mains', and 'adds' and creates a subplot of any
    data which matches from the header index. fig and ax are the already made
    fig and axes that we want to put the plot onto
    """

    # loop through list of adds
        # for each of adds values, create a mask on the data for [main, adds[i]]
            # plot this one line

    h_index = np.arange(0,np.shape(data)[1])

    if colors == 0:
        colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, len(adds))]

    for i in np.arange(0,len(adds)):
        ai = adds[i]

        # mask against both main and adds conditions
        mask = (headers[:,0] == main) & (headers[:,1] == ai)

        # stop loop if mask is empty
        if sum(mask) == 0:
            break

        # extract header_indices to plot from headers[mask][2]
        plot_index = int(headers[mask][0][2])

        # define label by default or by design
        if labels == 'adds':
            lb = adds[i]
        elif labels == 'main':
            lb = main[0]
        else:
            lb = labels[i]

        # plot
        y = data[:,plot_index,fi]
        dev = stdev[:,plot_index,fi]
        #ax[t].errorbar(Time, data_norm[t][i,:], stdev_norm[t][i,:], color=colors[i], label = headers_adds[i])
        AX.plot(Time, y, color=colors[i], label = lb)
        AX.fill_between(Time, y+dev, y-dev, color=colors[i], alpha=0.4)

        #AX.errorbar(Time,, , color=colors[i], label = ai)
    if leg == True:
        AX.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
        
    return

#####################################################
#%%

# Set numbering of headers_main and headers_adds

# Start by sorting headers_main and headers_adds lists
headers_main.sort()
headers_adds.sort()


header_check = 'n'
while header_check == 'n':

    print("Main Headers are :")
    print(headers_main)

    header_check = input("Would you like to keep this order? [y/n]: ")
    if header_check == 'y':
        new_headers = headers_main
        break

    print("This is the current order:")

    for i in np.arange(0,len(headers_main)):
        print(str(i) + ": " + headers_main[i])

    # give me the new order
    neworder = input("List the new order in the format of [3,0,2,1] with brackets, commas, and no spaces: ")

    # create a new order
    index = neworder[1:-1].split(",")
    new_headers = []
    for i in index:
        i = int(i)
        new_headers = np.append(new_headers, headers_main[int(i)])
        new_headers = list(new_headers)

    #check the new order
    print("Just to check, this is the new order you wanted?")
    for i in np.arange(0,len(new_headers)):
        print(str(i) + ": " + new_headers[i])
    header_check = input("Is that right? [y/n]: ")

headers_main = new_headers

# Set Adds

print("Ok, thanks!")

adds_check = 'n'
while adds_check == 'n':

    print("Main Adds are :")
    print(headers_adds)

    adds_check = input("Would you like to keep this order? [y/n]: ")
    if adds_check == 'y':
        new_adds = headers_adds
        break

    # this only continues if answer is not 'y'
    print("This is the current order:")

    for i in np.arange(0,len(headers_adds)):
        print(str(i) + ": " + headers_adds[i])

    # give me the new order
    neworder = input("List the new order in the format of [0,3,1,2] with brackets, commas, and no spaces: ")

    # create a new order
    index = neworder[1:-1].split(",")
    new_adds = []
    for i in index:
        i = int(i)
        new_adds = np.append(new_adds, headers_adds[int(i)])
        new_adds = list(new_adds)

    #check the new order
    print("Just to check, this is the new order you wanted?")
    for i in np.arange(0,len(new_adds)):
        print(str(i) + ": " + new_adds[i])
    adds_check = input("Is that right? [y/n]: ")

headers_adds = new_adds


#####################################################

#%%

# PLOT 1 - RAW DATA

# ONE FIGURE ONLY
# ONE SUBFIGURE FOR EACH 'HEADER_MAINS' CATEGORY


for i in fluor_index:
    i = int(i)

    H_len = len(headers_main)

    # start plots with one subplot for each header category
    if H_len > 1:
        fig, axs = plt.subplots(H_len, 1)
    else:
        fig, axs = plt.subplots()
    fig.set_size_inches(8, 4*H_len)

    for m in np.arange(0,len(headers_main)): # for each subplot
        hm = headers_main[m]

        if H_len > 1:
            AX = axs[m]
        else:
            AX = axs

        makeplot(data, stdev, hm, headers_adds, i, fig, AX)
        AX.set_title(hm)
        AX.set_xlabel('Minutes')
        AX.set_ylabel('RFU')

    fig.suptitle(fluor_name[int(fluor_index[i])], fontsize = 36)

    plt.tight_layout()

#####################################################
#%%

# PLOT 2 - RAW DATA

# ONE FIGURE FOR EACH 'HEADER_MAINS' CATEGORY

# plot one figure for each 'header_main' category
# each subplot is a different fluorophore

for m in np.arange(0,len(headers_main)): # for each figure

    HM = headers_main[m]
    F_len = len(fluor_index)
    # start plots with one subplot for each fluorophore
    fig, axs = plt.subplots(F_len)
    fig.set_size_inches(6*F_len, 4)
    if F_len == 1:
        fig.set_size_inches(8, 6)

    for i in fluor_index:
        i = int(i)

        if F_len > 1:
            AX = axs[i]
        else:
            AX = axs

        makeplot(data, stdev, HM, headers_adds, i, fig, AX)
        AX.set_title(fluor_name[i])

        AX.set_xlabel('Minutes')
        AX.set_ylabel('RFU')

    fig.suptitle(HM, fontsize = 24)
    plt.tight_layout()

#####################################################


#%%

# PLOT 3 - RAW DATA

# ONE FIGURE ONLY
# ALL DATA ON THIS FIGURE. ONE COLOR PER 'HEADER_MAINS' CATEGORY

# plot all data on one graph with 4 different colors for base concentrations
fig, AX = plt.subplots()
fig.set_size_inches(6, 4)

for m in np.arange(0,len(headers_main)): # for each concentration

    HM = headers_main[m]

    # loop through list of adds
        # for each of adds values, create a mask on the data for [main, adds[i]]
            # plot this one line

    h_index = np.arange(0,np.shape(data)[1])

    colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, len(headers_main))]

    for i in np.arange(0,len(headers_adds)):
        ai = headers_adds[i]
        # mask against both main and adds conditions
        mask = (headers[:,0] == HM) & (headers[:,1] == ai)

        # stop loop if mask is empty
        if sum(mask) == 0:
            break

        # extract header_indices to plot from headers[mask][2]
        plot_index = int(headers[mask][0][2])

        # plot
        # only label the first data point
        if i == 1:
            AX.plot(Time,data[:,plot_index,0], color=colors[m], label = HM)
        else:
            AX.plot(Time,data[:,plot_index,0], color=colors[m])

    AX.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

AX.set_title('CF Fluorescence of All data')
AX.set_xlabel('Minutes')
AX.set_ylabel('RFU')

plt.tight_layout()

plt.show()


#%%

# CORRECT DATA FOR PHOTOBLEACHING

# FIND THE DATA COLUMN THAT MATCHES HM AND HAS ADDS='HEPES'.
# FROM THIS DATA, FIND THE DELTA THAT IS FROM PHOTOBLEACHING (DATA - DATA[0])
# SUBTRACT THIS DELTA FROM EACH COLUMN THAT ALSO MATCHES HM

controlstring = 'HEPES'

data_ctrl = data.copy()

for HM in headers_main:

    # make mask for header
    mask_HM = headers[:,0] == HM

    # make mask for control data
    mask = (headers[:,0] == HM) & (headers[:,1] == controlstring)

    # subtract control data from all other data with matching header
    # extract header_indices to plot from headers[mask][2]
    ctrl_index = int(headers[mask][0][2])

    # extract the photobleaching control data this corresponds to as 'delta'
    delta = data_ctrl[:,ctrl_index,0] - data_ctrl[0,ctrl_index,0]

    # for every column to fix, subtract the photobleaching data
    for i in np.arange(0,len(headers[mask_HM])): # these are all the columns that need to be corrected

        data_index = int(headers[mask_HM][i][2])

        # regular raw data
        data_ctrl[:,data_index,0] = data_ctrl[:,data_index,0] - delta

#%%

# NORMALIZE DATA

# normalize regular raw data
data_norm  = np.zeros(np.append(np.shape(data),2))
stdev_norm = np.zeros(np.append(np.shape(data),2))

data_norm[:,:,:,0]  = data / data[0,:,:] # normalized to t0
data_norm[:,:,:,1]  = data / data[-1,:,:] # normalized to tf
stdev_norm[:,:,:,0] = stdev / data[0,:,:] # normalized to t0
stdev_norm[:,:,:,1] = stdev / data[-1,:,:] # normalized to tf

# normalized photobleach corrected data
data_norm_ctrl  = np.zeros(np.append(np.shape(data_ctrl),2))
stdev_norm_ctrl = np.zeros(np.append(np.shape(data_ctrl),2))

data_norm_ctrl[:,:,:,0]  = data_ctrl / data_ctrl[0, :,:] # normalized to t0
data_norm_ctrl[:,:,:,1]  = data_ctrl / data_ctrl[-1,:,:] # normalized to tf


#%%
# PLOT 4 - NORMALIZED DATA (NOT PHOTOBLEACH CORRECTED)

# ONE FIGURE FOR EACH 'HEADER_MAINS' CATEGORY

# BOTH NORMALIZATION TO T0 AND TF ARE CALCULATED. ONLY T0 IS PLOTTED

for m in np.arange(0,len(headers_main)): # for each figure

    HM = headers_main[m]
    # start plots with one subplot for each normalization
    fig, axs = plt.subplots()
    fig.set_size_inches(8, 6)

    X = [0] # CHANGE TO [1] FOR NORMALIZING TO TF

    for t in [0]:

        if len(X) == 1:
            AX = axs
        else:
            AX = axs[t]

        makeplot(data_norm[:,:,:,t], stdev_norm[:,:,:,t], HM, headers_adds, 0, fig, AX)
        AX.set_xlabel('Minutes')
        AX.set_ylabel('RFU')
    if t == 0:
        AX.set_title(r'CF Emission normalized to $t_0$')
    elif t==1:
        AX.set_title(r'CF Emission normalized to $t_f$')
    fig.suptitle(HM, fontsize = 24)
    plt.tight_layout()

#%%
# PLOT 5 - NORMALIZED DATA (NOT PHOTOBLEACH CORRECTED)

# ONE FIGURE FOR EACH 'HEADER_ADDS' CATEGORY

for a in np.arange(0,len(headers_adds)): # for each figure

    ai = headers_adds[a]
    # start plots with one subplot for each normalization
    fig, axs = plt.subplots(2)
    fig.set_size_inches(8, 8)

    colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, len(headers_main))]

    for m in np.arange(0,len(headers_main)):
        HM = headers_main[m]

        for t in [0,1]:

            makeplot(data_norm[:,:,:,t], stdev_norm[:,:,:,t], [HM], [ai], 0, fig, axs[t], [colors[m]], labels='main')
            axs[t].set_xlabel('Minutes')
            axs[t].set_ylabel('RFU')
        axs[0].set_title(r'CF Emission normalized to $t_0$')
        axs[1].set_title(r'CF Emission normalized to $t_f$')
    fig.suptitle(ai, fontsize = 24)
    plt.tight_layout()

#%%
# PLOT 6 - RAW DATA

# ONE FIGURE FOR EACH 'HEADER_ADDS' CATEGORY

for a in np.arange(0,len(headers_adds)): # for each figure

    ai = headers_adds[a]
    # start plots with one subplot for each normalization
    fig, axs = plt.subplots()
    fig.set_size_inches(8, 8)

    colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, len(headers_main))]

    for m in np.arange(0,len(headers_main)):
        HM = headers_main[m]

        makeplot(data[:,:,:], stdev[:,:,:], [HM], [ai], 0, fig, axs, [colors[m]], labels='main')
        axs.set_xlabel('Minutes')
        axs.set_ylabel('RFU')
    fig.suptitle(ai, fontsize = 24)
    plt.tight_layout()



#%%

# PLOT 7 - PHOTOBLEACH CORRECTED DATA

# 4 SUBLOTS ON ONE FIGURE

for i in fluor_index:
    i = int(i)

    H_len = len(headers_main)

    # start plots with one subplot for each header category
    if H_len > 1:
        fig, axs = plt.subplots(H_len, 1)
        fig.set_size_inches(8, 4*H_len)
    else:
        fig, axs = plt.subplots()
        fig.set_size_inches(8, 8)

    for m in np.arange(0,len(headers_main)): # for each subplot
        hm = headers_main[m]

        if H_len > 1:
            AX = axs[m]
        else:
            AX = axs

        makeplot(data_ctrl, stdev, hm, headers_adds, i, fig, AX)
        AX.set_title(hm)
        AX.set_xlabel('Minutes')
        AX.set_ylabel('RFU')

    fig.suptitle('Data Corrected for Photobleaching', fontsize = 20)

    plt.tight_layout()

#%%

# PLOT 8 - PHOTOBLEACH CORRECTED DATA

# ONE FIGURE FOR EACH 'HEADER_MAIN' CATEGORY

for m in np.arange(0,len(headers_main)): # for each figure

    HM = headers_main[m]
    F_len = len(fluor_index)
    # start plots with one subplot for each fluorophore
    fig, axs = plt.subplots(F_len)
    fig.set_size_inches(6*F_len, 4)
    if F_len == 1:
        fig.set_size_inches(8, 6)

    for i in fluor_index:
        i = int(i)

        if F_len > 1:
            AX = axs[i]
        else:
            AX = axs

        makeplot(data_ctrl, stdev, HM, headers_adds, i, fig, AX)
        AX.set_title(HM + ' - Corrected for Photobleaching', fontsize = 20)

        AX.set_xlabel('Minutes')
        AX.set_ylabel('RFU')

    plt.tight_layout()

#####################################################

#%%

# PLOT 9 - NORMALIZED AND PHOTOBLEACH CORRECTED DATA

# ONE FIGURE FOR EACH 'HEADERS_MAIN' CATEGORY


for m in np.arange(0,len(headers_main)): # for each figure

    HM = headers_main[m]
    # start plots with one subplot for each normalization
    fig, axs = plt.subplots()
    fig.set_size_inches(8, 6)

    X = [0] # CHANGE TO 1 FOR NORMALIZAING BY TF

    for t in X:

        if len(X) == 1:
            AX = axs
        else:
            AX = axs[t]

        if HM == headers_main[-1]: # the last one, then plot the legend
            makeplot(data_norm_ctrl[:,:,:,t], stdev_norm[:,:,:,t], HM, headers_adds, 0, fig, AX)
        else:
            makeplot(data_norm_ctrl[:,:,:,t], stdev_norm[:,:,:,t], HM, headers_adds, 0, fig, AX, leg='False')
        AX.set_xlabel('Minutes')
        AX.set_ylabel('RFU')
    if t == 0:
        AX.set_title(r'CF Emission normalized to $t_0$' + ' \n and corrected for photobleaching')
    elif t == 1:
        AX.set_title(r'CF Emission normalized to $t_f$' + ' \n and corrected for photobleaching')
    fig.suptitle(HM, fontsize = 24)
    plt.tight_layout()

#%%

# PLOT 10 - NORMALIZED AND PHOTOBLEACH CORRECTED DATA

# ONE SUBFIGURE FOR EACH 'HEADERS_MAIN' CATEGORY

fig, axs = plt.subplots(1,3)
fig.set_size_inches(14, 6)

for m in np.arange(0,len(headers_main)): # for each figure

    HM = headers_main[m]
    # start plots with one subplot for each normalization

    AX = axs[m]


    X = [0] # CHANGE TO 1 FOR NORMALIZAING BY TF

    for t in X:

        if HM == headers_main[-1]: # the last one, then plot the legend
            makeplot(data_norm_ctrl[:,:,:,t], stdev_norm[:,:,:,t], HM, headers_adds, 0, fig, AX)
        else:
            makeplot(data_norm_ctrl[:,:,:,t], stdev_norm[:,:,:,t], HM, headers_adds, 0, fig, AX, leg='False')

        if HM == headers_main[0]: # the first one, give y axis label
            AX.set_ylabel('RFU')
        AX.set_xlabel('Minutes')

        AX.set_title(HM)
        AX.set_ylim([0.4,1.3])
    if t == 0:
        fig.suptitle(r'CF Emission normalized to $t_0$ and corrected for photobleaching')
    elif t == 1:
        fig.suptitle(r'CF Emission normalized to $t_f$ and corrected for photobleaching')


    plt.tight_layout()


#%%

# PLOT 11 - ALL DATA

# ONE SUBFIGURE FOR EACH OF 'RAW', 'NORM', 'NORM+CORRBLEACH'

# ONE COLOR FOR EACH TYPE OF 'HEPES','CF','CF-LUV'

# plot all data on one graph with 4 different colors for base concentrations
fig, axs = plt.subplots(1,3)
fig.set_size_inches(12, 4)

for m in np.arange(0,len(headers_main)): # for each concentration

    HM = headers_main[m]

    # loop through list of adds
        # for each of adds values, create a mask on the data for [main, adds[i]]
            # plot this one line

    h_index = np.arange(0,np.shape(data)[1])

    colors = [plt.cm.viridis(x) for x in np.linspace(0, 1, len(headers_main))]

    for i in np.arange(0,len(headers_adds)):
        ai = headers_adds[i]
        # mask against both main and adds conditions
        mask = (headers[:,0] == HM) & (headers[:,1] == ai)

        # stop loop if mask is empty
        if sum(mask) == 0:
            break

        # extract header_indices to plot from headers[mask][2]
        plot_index = int(headers[mask][0][2])

        # plot
        dd = [data, data_norm[:,:,:,0], data_norm_ctrl[:,:,:,0]]
        
        for d_index in np.arange(0,len(dd)):
            d_index = int(d_index)
            mydata = dd[d_index]
            AX = axs[d_index]            
            
            # only label the first data point
            if i == 1:
                AX.plot(Time,mydata[:,plot_index,0], color=colors[m], label = HM)
            else:
                AX.plot(Time,mydata[:,plot_index,0], color=colors[m])

        AX.set_xlabel('Minutes')
        AX.set_ylabel('RFU')
    AX.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

axs[0].set_title(r'CF Emission')
axs[1].set_title(r'CF Emission normalized to $t_0$')
axs[2].set_title(r'CF Emission normalized to $t_0$' + ' \n and corrected for photobleaching')

plt.tight_layout()
plt.show()
