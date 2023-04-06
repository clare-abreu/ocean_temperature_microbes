This folder contains the code used in Abreu & Dal Bello et al (Science Advances, 2023).

The included files are:


assign_copynums_functions.py:
These functions assign copy numbers to taxa at the highest-resolution taxonomic level available from the rrnDB, using the mean from that table. Most datasets do not contain species classifications, but the first function will handle species if they are included. The second function does not, but is otherwise identical.


calc_WMCN_WMGR_functions.py:
These functions calculate WMCN (and WMGR) for each sample in a dataset. The datasets used for these functions are in the Generic Data folder of data/.
