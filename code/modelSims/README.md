This folder contains code for simulating the Lotka-Volterra (LV) model, as well as a linear consumer resource (CR) model (the mathematical details of these models are described in the Methods section and in the Supplementary Text).

The included files are:


# modelFns.py:
This file contains all of the functions used to simulate the LV and CR models, as well as a function to fit the LV simulations with the best value of p, which is the geometric distribution parameter representing the fraction of the starting community with rRNA copy number = 1, to each of the datasets. 


# fit_p.py:
This file fits each dataset to the model and plots the results.


# simLV.py:
This file simulates the Lotka-Volterra model over a range of temperatures for 100 species. The weighted mean copy number, which is a measure of the distribution of fast and slow growers in the equilibrated community, is plotted at each temperature, showing that warmer temperatures favor slower growers.


# simCR.py:
This file simulates the consumer resource model over a range of temperatures and nutrient supply concentrations for any number of species and resources. The weighted mean maximum growth rate, which is a measure of the distribution of fast and slow growers in the equilibrated community, is plotted at each environmental condition, showing that warmer temperatures favor slower growers.
