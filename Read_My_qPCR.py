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

wdir    = r"C:\Users\Alison\Documents\AL Data\B2P61\qPCR"
fname   = r"C:\Users\Alison\Documents\AL Data\B2P61\qPCR\2012-09-14 CF-GUV Overnight -  Quantification Amplification Results_FAM.csv"

fname_h = r"C:\Users\Alison\Documents\AL Data\B2P61\qPCR\2012-09-14 CF-GUV Overnight -  Headers.csv"

AverageDatainTriplicates = True

t_per_run = 2 # minutes

fluor_name = ['FAM', 'Texas Red', 'Cal Gold 540']

fluor_FAM = True
fluor_TexasRed = True
fluor_CalGold = False

#####################################################

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

#####################################################

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
    
    
#####################################################

# Function for making plots

def makeplot(main, adds, fi, FIG, AX):

    """
    This function takes the 'mains', and 'adds' and creates a subplot of any 
    data which matches from the header index. fig and ax are the already made 
    fig and axes that we want to put the plot onto
    """
    
    # loop through list of adds
        # for each of adds values, create a mask on the data for [main, adds[i]]
            # plot this one line
   
    h_index = np.arange(0,np.shape(data)[1])
    
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
    
        # plot
        AX.errorbar(Time,data[:,plot_index,fi], stdev[:,plot_index,fi], color=colors[i], label = ai)
    
    AX.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    return
    
#####################################################

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

#%%

# Set numbering of headers_main and headers_adds



header_check = 'n'
while header_check == 'n':
    
    print("Main Headers are :")
    print(headers_main)
    
    header_check = input("Would you like to keep this order? [y/n]: ")
    if header_check == 'y':
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

# PLOT 1

# plot one figure for each fluorophore
# each subplot is a different 'headers_main' category

for i in fluor_index:
    i = int(i)
    
    H_len = len(headers_main)

    # start plots with one subplot for each header category    
    fig, axs = plt.subplots(H_len,1)
    fig.set_size_inches(8, 4*H_len)

    for m in np.arange(0,len(headers_main)): # for each subplot
        hm = headers_main[m]    
        makeplot(hm, headers_adds, i, fig, axs[m])
        axs[m].set_title(hm)
        
    fig.suptitle(fluor_name[int(fluor_index[i])], fontsize = 36)
    plt.tight_layout()

#####################################################

# PLOT 2

# plot one figure for each 'header_main' category
# each subplot is a different fluorophore

for m in np.arange(0,len(headers_main)): # for each figure
    
    HM = headers_main[m]    
    F_len = len(fluor_index)
    # start plots with one subplot for each fluorophore
    fig, axs = plt.subplots(1, F_len)
    fig.set_size_inches(8*F_len, 4)
    
    for i in fluor_index:
        i = int(i)
        makeplot(HM, headers_adds, i, fig, axs[i])
        axs[i].set_title(fluor_name[i])
        
    fig.suptitle(HM, fontsize = 36)
    plt.tight_layout()

#####################################################


#%%
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