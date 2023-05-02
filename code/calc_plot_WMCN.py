#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script calculates the weighted mean copy number (MCN, or a measure of the 
distribution of growth rates, of the free-living datasets, with the exception
of LMO, which is too large to include in Github. You can download it from Dryad
here: https://datadryad.org/stash/share/LLaE8Wv1Wmo6f0M0nlL7S7U2kIfdj1OH9LHfnmKv8mA

The WMCN vs. temperature trend is then plotted for these datasets. Note that
these plots do show a negative trend as predicted by our model, but they do 
not include the multivariate analysis included in the paper, which controls
for the effects of nutrients and other environmental variables. See 
the R notebook analiz_environmental_data.Rmd for this.

Written by Clare Abreu for Abreu & Dal Bello et al, Science Advances (2023)
"""

import matplotlib.pyplot as plt
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
# Import the function and data:
calc_fns_path = './calc_WMCN_WMGR_functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(calc_fns_path)))
import calc_WMCN_WMGR_functions as calc
data_path = '../data/genericData/'
all_data = ['AOT_free','TARA','Med','SPOT'] # Add 'LMO_free' here if you've downloaded the abundance data from Dryad

# Calculate WMCN and plot the results
for this_data in all_data:
    calc_df = calc.calc_wmcn_only(this_data,data_path)
    plt.scatter(calc_df['Temperature'],calc_df['WMCN'])
    plt.ylabel('WMCN')
    plt.xlabel('Temperature')
    plt.title(this_data)
    plt.show()