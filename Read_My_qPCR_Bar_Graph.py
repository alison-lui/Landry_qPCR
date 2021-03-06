# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 18:42:58 2021

@author: Alison


This Script takes identically structured fluorescence data from multiple files
(ideally, the same sample assessed on different days), averages over the entire
qPCR timecourse of each well, and plots the following:
    1) a bar graph of the fluorescence in each well over each day. 
        each fluorophore gets its own subplot.
    2) a line plot of FAM / Texas Red fluoresence for each well over the course
        of sample collection days. This should show FAM leaking while Texas Red
        stays about the same.

    if only one fluorophore is selected, graph 2 is not created.
    
This file, ready to commit to main.

"""

# Header Filename
fname_h = r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P95\qPCR\2021-10-28_CF_dilutions -  Headers.csv"

# Data filenames
fnames = [r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P95\qPCR\2021-10-28_CF_dilutions -  Quantification Amplification Results_FAM.csv"]

# Data working directories
wdirs =  [r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P95\qPCR"]


fluor_FAM = True
fluor_TexasRed = False
fluor_CalGold = False

useLatex = False

#######################################

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from sklearn.metrics import r2_score
 
#######################################  
# use LaTeX ?

if useLatex == True:
    # use LaTeX fonts in the plot
    mpl.rcParams['text.usetex'] = True
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica"]})
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

#######################################  

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

fluor_colors = [c_Cal, c_TxR, c_FAM]

########################################

def extractdata(wdir, fname, fluor_index):
    
    " Pull out data for one csv file"
    
    # extract header title
    base = fname[:-7]
    title = os.path.basename(base)
    title = title.split(' -  ')
    title = title[0]
 
   # Get Data
    os.chdir(wdir)
    
    # how many fluorophore files?
    N_fluor = len(fluor_index)
    
    # extract the first emission data so we know what size to make the array
    base = fname[:-7]
    df = pd.read_csv(fname)
    data = df.to_numpy()
    H, N = np.shape(data)
    
    # put fluorophore data into empty emission data array
    data = np.empty((H,N,N_fluor))
    for x in range(0,N_fluor):
        i = int(fluor_index[x])
        temp = pd.read_csv(base + fluor_name[i] + ".csv")
        temp = temp.to_numpy()
        data[:,:,x] = temp
    
    # Emission data sits in 3rd and later columns
    data = data[:,2:,:]
    H,N,N_fluor = np.shape(data)

    return H, N, N_fluor, data, title
    
#####################################################

# organize input data

# extract header data
dfh = pd.read_csv(fname_h)
headers = dfh.to_numpy()
headers = headers[0]

# list fluorophores
fluor_name = ['FAM', 'Texas Red', 'Cal Gold 540']
fluor = [fluor_FAM, fluor_TexasRed, fluor_CalGold]
N_fluor = sum(1 for x in fluor if x == True)
# remove any unused fluorophores in name list
fluor_index = []
for x in range(0,len(fluor)):
    if fluor[x] == True:
        fluor_index = np.append(fluor_index, x)

# pull out the first file name's data
H, N, N_fluor, data, title = extractdata(wdirs[0], fnames[0], fluor_index)

# create empty array for emission and legend data
composite = np.zeros((H,N,N_fluor,len(fnames)))
leg_label = []

# run extractdata and put it into the composite array
for i in np.arange(0,len(fnames)):
    H, N, N_fluor, composite[:,:,:,i], templabel = extractdata(wdirs[i], fnames[i], fluor_index)
    leg_label = np.append(leg_label, templabel)
    
#####################################################
#%% FIGURE 1

# Create bar graph with one subplot for each fluorophore

# average data along all time points
bar_composite = np.average(composite, axis=0)

fig, axs = plt.subplots(N_fluor)
fig.set_size_inches(8, 4*N_fluor)


X = np.arange(N)

for i in fluor_index:
    i = int(i)
    
    colors = [colorFader(fluor_colors[i], c_wht, mix=x) for x in np.linspace(0, 1, len(fnames)+2)]
    
    if N_fluor == 1:
        AX = axs
    else:
        AX = axs[i]
    
    
    for f in np.arange(0,len(fnames)):
        AX.bar(X + f/(len(fnames)+1), bar_composite[:,i,f], color = colors[f], width = 0.25, label = leg_label[f])
    
    AX.set_xticks(X + 1/(len(fnames)+1)) 
    AX.set_xticklabels(headers, rotation=90)
    
    AX.set_title(fluor_name[i])
    AX.legend() 
           
plt.tight_layout()

###############################################
#%% FIGURE 2

# If only one data set given, plot two fluorophores side by side on bar graph

if len(fnames) == 1:
    
    fig, axs = plt.subplots()
    for f in fluor_index:
        f = int(f)
        axs.bar(X + f*0.25, bar_composite[:,f,0], color = fluor_colors[f], width = 0.25, label = fluor_name[f])

    axs.set_xticks(X + 1/(len(fnames)+1)) 
    axs.set_xticklabels(headers, rotation=90)
    
    axs.set_title(leg_label[0])
    axs.legend() 
           
    plt.tight_layout()


###############################################

#%% FIGURE 3

# Calculate ratio of CF to DF

DF = [1,0.8,0.6,0.4,0.2,0.1,0,0,0,0,0,0,0]

# subtract the average of the zero values 
zeroaverage = np.mean(bar_composite[-7:])
barsubtracted = bar_composite - zeroaverage
barsubtracted = barsubtracted.flatten()

# fit all
xdatafull = DF[0:7]
ydatafull = barsubtracted[0:7]
fit = np.polyfit(xdatafull, ydatafull, deg=1, full=False)
yfitfull = [x * fit[0] + fit[1] for x in xdatafull]

# make a line using the first 5 points
xdata = DF[2:7]
ydata = barsubtracted[2:7]
fit = np.polyfit(xdata, ydata, deg=1, full=False)
yfit = [x * fit[0] + fit[1] for x in xdata]
rr = r2_score(ydata, yfit)


fig,ax = plt.subplots()

ax.scatter(DF,bar_composite, color='b')
#ax.plot(xdatafull,yfitfull, color='b')
ax.scatter(DF,barsubtracted, color='k')
ax.plot(xdata,yfit, color='k')

ax.set_xlabel('Dilution Factor')
ax.set_ylabel('RFU')
stringlabel = [r"$R^2$ = " + str(np.round(rr,decimals=3)) 
                + '\n y = ' + str(np.round(fit[0],decimals=0)) 
                + '*x + ' + str(np.round(fit[1],decimals=3))]
ax.annotate(stringlabel[0], xy=(np.mean(DF)*1.5,np.mean(ydata)))

fig.tight_layout()

plt.show()
    

#%% FIGURE 4


if N_fluor == 1:
    sys.exit()

# Compute FAM RFU / Texas Red RFU and scatter plot for each fraction 

# first mask out anyones that have:
# FAM above 60000 (oversaturated)
# Texas Red below 2500 (undersaturated)
mask = (bar_composite[:,0,0] < 60000) & (bar_composite[:,1,0] > 2500)

ratios = np.divide(bar_composite[:,0,:], bar_composite[:,1,:])


# plot with each sample fraction is a different color

colors = [plt.cm.rainbow(x) for x in np.linspace(0, 1, len(headers[mask]))]

fig, ax = plt.subplots()

XN = np.arange(0,len(headers[mask]))

ax.set_xticks(XN+0.5) 
ax.set_xticklabels(headers[mask], rotation=90)

for n in XN:
    x = np.linspace(n, n+1, len(fnames))
    y = ratios[mask][n,:]
    ax.plot(x, y, color=colors[n])
    plt.gca().get_xticklabels()[n].set_color(colors[n])
    


ax.set_title('FAM / Texas Red Fluorescence over the days')
ax.set_ylabel('FAM / Texas Red Fluorescence')

plt.tight_layout()

        
    
    
    
    
    
    
    
    
    
    

    
    
