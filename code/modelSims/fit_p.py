#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to call the boostrap_p function to fit the data to the LV model using only
one parameter: the value of p in the geometric distribution, which represents the 
fraction of the starting community with rRNA = 1.

Written by Clare Abreu for Abreu & Dal Bello et al, Science Advances (2023)
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
plt.style.use("ggplot")
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 18
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
import os
import sys

modelFns_path = './modelFns.py'
sys.path.append(os.path.dirname(os.path.expanduser(modelFns_path)))
import modelFns as model

data_path = '../../data/'
sim_df = pd.read_csv(f'{data_path}simData/all_simulated_WMCNs_deltafit.csv')
sim_df.set_index('Temperature',inplace=True)
metadata = pd.read_csv(f'{data_path}Mega_Table_metadata.csv')
metadata.set_index('Sample',inplace=True)
p_range=np.unique(sim_df['p'])

# %%
def round_to_nearest(num, num_list):
    # Function to round mean estimated result to nearest simulated parameter
    min_diff = float('inf')
    nearest_num = None
    for n in num_list:
        diff = abs(num - n)
        if diff < min_diff:
            min_diff = diff
            nearest_num = n
    return nearest_num

#%%
# Fit the AOT dataset:
AOT_metadata = metadata.filter(like='AOT', axis=0)
AOT_metadata = AOT_metadata[AOT_metadata['Filter Fraction']=='a']   # Include only free-living filter fraction
AOT_metadata = AOT_metadata[~AOT_metadata['Temperature'].isna()]
AOT_metadata = AOT_metadata[~AOT_metadata['WMCN'].isna()]
AOT_metadata = AOT_metadata[['WMCN','Temperature','Latitude']]
sample_size=len(AOT_metadata)
num_trials = 100      # 1000 trials were used for Fig. 2; can reduce to decrease computational time
[AOT_mu,AOT_SD,AOT_mean_min_error] = model.bootstrap_p(p_range,sample_size,num_trials,AOT_metadata,sim_df)

# Plot the result:
plt.figure(figsize=(10,7))
plt.scatter(AOT_metadata['Latitude'],AOT_metadata['WMCN'],c = AOT_metadata['Temperature'])
ax = plt.gca()
ax.set_facecolor("white")
plt.setp(ax.spines.values(), color='k')

AOT_plot_p = round_to_nearest(AOT_mu, p_range)
sim_df_AOT = sim_df[sim_df['p']==AOT_plot_p]
sim_WMCN = []
lats=[]
temps=[]
for k in range(len(AOT_metadata)):
    this_temp = np.round(AOT_metadata.iloc[k].Temperature,2)
    sim_WMCN.append(sim_df_AOT.loc[this_temp, 'WMCN'])
    lats.append(AOT_metadata.iloc[k].Latitude)
    temps.append(this_temp)
    
d = {'Latitude': lats, 'WMCN': sim_WMCN}
AOT_ave = pd.DataFrame(data=d)
grouped_df = AOT_ave.groupby("Latitude")
mean_df = grouped_df.mean()
mean_df = mean_df.reset_index()

plt.plot(mean_df["Latitude"],mean_df['WMCN'],'k')
plt.xlabel('Latitude')
plt.ylabel('MCN')
plt.title('AOT')
plt.show()

#%%
# Fit the LMO dataset:
LMO_metadata = metadata.filter(like='LMO', axis=0)
LMO_metadata = LMO_metadata[LMO_metadata['Filter Fraction']=='3-0.2']
LMO_metadata = LMO_metadata[~LMO_metadata['Temperature'].isna()]
LMO_metadata = LMO_metadata[~LMO_metadata['WMCN'].isna()]
LMO_metadata = LMO_metadata[['WMCN','Temperature','Date']]
sample_size=len(LMO_metadata)
num_trials = 100      # 1000 trials were used for Fig. 2; can reduce to decrease computational time
[LMO_mu,LMO_SD,LMO_mean_min_error] = model.bootstrap_p(p_range,sample_size,num_trials,LMO_metadata,sim_df)

# Plot the result:
LMO_metadata = LMO_metadata.sort_values(by='Date')
LMO_metadata.set_index('Date',inplace=True)
LMO_metadata.index = pd.to_datetime(LMO_metadata.index)

plt.figure(figsize=(12,7))
plt.scatter(LMO_metadata.index,LMO_metadata['WMCN'],c = LMO_metadata['Temperature'])
ax = plt.gca()
ax.set_facecolor("white")
plt.setp(ax.spines.values(), color='k')

LMO_plot_p = round_to_nearest(LMO_mu, p_range)
sim_df_LMO = sim_df[sim_df['p']==LMO_mu]
sim_WMCN = []
dates=[]
for k in range(len(LMO_metadata)):
    this_temp = np.round(LMO_metadata.iloc[k].Temperature,2)
    sim_WMCN.append(sim_df_LMO.loc[this_temp, 'WMCN'])
    dates.append(LMO_metadata.index[k])
    
d = {'Date': dates, 'WMCN': sim_WMCN}
LMO_ave = pd.DataFrame(data=d)

plt.plot(LMO_ave["Date"],LMO_ave['WMCN'],'k')
plt.xlabel('Date')
plt.ylabel('MCN')
plt.title('LMO')
plt.show()


#%%
# Fit the TARA dataset:
TARA_metadata = metadata.filter(like='TARA', axis=0)
TARA_metadata = TARA_metadata[TARA_metadata['Latitude']<0.001]
TARA_metadata = TARA_metadata[TARA_metadata['Latitude']>-40.001]
TARA_metadata = TARA_metadata[~TARA_metadata['Temperature'].isna()]
sample_size=len(TARA_metadata)
num_trials = 100      # 1000 trials were used for Fig. 2; can reduce to decrease computational time
[TARA_mu,TARA_SD,TARA_mean_min_error] = model.bootstrap_p(p_range,sample_size,num_trials,TARA_metadata,sim_df)

# Plot the result:
TARA_metadata = TARA_metadata[~TARA_metadata['WMCN'].isna()]
TARA_metadata = TARA_metadata[['WMCN','Temperature','Depth']]
plt.figure(figsize=(10,7))
plt.scatter(TARA_metadata['Depth'],TARA_metadata['WMCN'],c = TARA_metadata['Temperature'],facecolor='w')
ax = plt.gca()
ax.set_facecolor("white")
ax.set_facecolor("white")
plt.setp(ax.spines.values(), color='k')

TARA_plot_p = round_to_nearest(TARA_mu, p_range)
sim_df_TARA = sim_df[sim_df['p']==TARA_plot_p]
sim_WMCN = []
depths=[]
for k in range(len(TARA_metadata)):
    this_temp = np.round(TARA_metadata.iloc[k].Temperature,2)
    sim_WMCN.append(sim_df_TARA.loc[this_temp, 'WMCN'])
    depths.append(TARA_metadata.iloc[k].Depth)
    
d = {'Depth': depths, 'WMCN': sim_WMCN}
TARA_ave = pd.DataFrame(data=d)
grouped_df = TARA_ave.groupby("Depth")
mean_df = grouped_df.mean()
mean_df = mean_df.reset_index()
mean_df.to_csv('/Users/clare/Dropbox/WMCN macroecological patterns/WMCN analyses/Fig2 data from Clare/Fits.5.22/TARA_fit.csv',index=False)


plt.plot(mean_df["Depth"],mean_df['WMCN'],'k')
plt.xlabel('Depth')
plt.ylabel('MCN')
plt.title('TARA')
plt.show()