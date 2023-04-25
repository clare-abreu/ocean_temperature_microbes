This folder contains the code used in Abreu & Dal Bello et al (Science Advances, 2023).

The subfolder # modelSims contains functions and code for simulating the Lotka-Volterra model and a consumer resource model.

The files include here are:

# assign_copynums_functions.py:
These functions assign copy numbers to taxa at the highest-resolution taxonomic level available from the rrnDB, using the mean from that table. Most datasets do not contain species classifications, but the first function will handle species if they are included. The second function does not, but is otherwise identical.


# calc_WMCN_WMGR_functions.py:
These functions calculate WMCN (and WMGR) for each sample in a dataset. The datasets used for these functions are in the Generic Data folder of data/.


# analiz_environmental_data.Rmd:
These notebook contains the workflow to use Generalized Additive Models and compute the effect of available environmental predictors on WMCN and WMGR. This notebook reproduce the results reported in Fig. 3 and Table S1. The imported dataset is available in the Generic Data folder as X.


