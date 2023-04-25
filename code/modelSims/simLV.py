#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to simulate the Lotka-Volterra model for 100 species over a range
of temeperatures and plot weighted mean copy number at each temperature.
This is the kind of plot shown in Figure 2B. Below, the following can be 
edited: the mean and SD of their interaction strength, the global death rate, 
and the fraction of the starting community with rRNA = 1. Copy numbers 1-10 
are included and drawn from a geometric distribution in the imported function 
Arr_100species_ivp.

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

#%% Simulate LV model
num_trials = 500
p=0.8 # This represents the fraction of the starting community with rRNA copy number = 1
delta = 0.05 # Global death rate
mu = 0.5 # Mean competition coeffcient (<1 required for coexistence)
sigma = 0.25 # SD of competition coefficients
[T,WMCN,CN_rel_abunds] = model.Arr_100species_ivp(delta,mu,sigma,num_trials,p)
T_C = T
T_C[:] = [T_Cn - 273 for T_Cn in T_C] # Convert Kelvin to Celsius

# Plot every trajectory:
for o in range(num_trials):
    plt.plot(T_C,WMCN[o,:],c='b',alpha=0.1)
    plt.scatter(T_C,WMCN[o,:],c='b',alpha=0.1)
    
# Plot the average:
WMCN_ave = np.mean(WMCN,0)
plt.plot(T_C,WMCN_ave,c='b',alpha=1,linewidth=5.0)

plt.xlabel("Temperature")
plt.ylabel("MCN")
plt.title("100 species, 10 copy numbers, %d trials" %(num_trials))
ax = plt.gca()
ax.set_facecolor("white")
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.show()