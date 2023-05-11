# ocean_temperature_microbes

This repository contains the data and code used in "Warmer temperatures favor slower-growing bacteria in natural marine communities," by Abreu & Dal Bello et al (Science Advances, 2023: https://www.science.org/doi/10.1126/sciadv.ade8352).

The repository contains two folders: code and data.

The data folder includes raw taxa data from marine datasets (with the exception of the large LMO dataset, which is available at Dryad: https://datadryad.org/stash/share/LLaE8Wv1Wmo6f0M0nlL7S7U2kIfdj1OH9LHfnmKv8mA ), environmental metadata from these datasets, ribosomal RNA operon (rRNA) copy number estimates for taxa downloaded from rrnDB (https://rrndb.umms.med.umich.edu/), computed weighted mean copy numbers of samples in these datasets, and simulated data generated from Lotka-Volterra model simulations, which we used to fit the model to the datasets.

The code folder contains functions that can map taxa data to rRNA copy numbers and calculate weighted mean copy numbers, as well as a subfolder with code that simulates the Lotka-Volterra model and consumer resource model and plots the results (including fits of the LV model to the datasets).
