#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to simulate the linear consumer resource model for any number of 
species and resources over a range of temperatures and nutrient supply
concentrations, and plot the weighted mean max growth rate at each 
temperature and nutrient concentration. This is the kind of plot shown in 
Figure S4. 
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

#%% Simulate linear consumer resource model

num_species = 15 # Number of species
num_resource = 10 # Number of resources
num_hours = 1000 # Number of simulated hours per trial. Need enough for equilibration.
start_density = .1 # Starting density of all species
delta = 0.1 # Global death/dilution rate
num_trials = 250 # Number of trials
T = [276,278,280,282,284,286,288,290,292,294,296,298]#,295]
T_C = [x - 273 for x in T]
r_conc = [1,1.5,2,2.5,3,3.5,4,4.5] # Supply concentration of all resources

# Simulate the model:
[WMGR_rm_ave,xticks,yticks] = model.get_WMGR_ave(num_species,num_resource,num_hours,start_density,delta,num_trials,T,r_conc)

# Plot the result:
fig, ax = plt.subplots(1,1)
img = ax.imshow(WMGR_rm_ave,extent=[-1,1,-1,1],origin='lower')

x_label_list = T_C
y_label_list = r_conc
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(x_label_list)
ax.set_yticklabels(y_label_list)
ax.set_xlabel('Temperature (C)')
ax.set_ylabel('Nutrient Supply')
fig.colorbar(img)
plt.show()