#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for Lotka-Volterra (LV) and consumer resources simulations,
as well as other functions calling the simulation functions, and a 
bootstrapping function to fit real data to the LV model

Written by Clare Abreu for Abreu & Dal Bello et al, Science Advances (2023)
"""


def dX_dt_GLV(t, X, r, alpha, delta, N):
# LV Generalized Species Number
# This function can simulate any number of species and outputs a time series of species
# abundances. It should be called like this: 
# X_final = solve_ivp(dX_dt_GLV, t_span, X_initial, args=(r,alpha,delta,N)), where:
    # X_initial: length-N 1D array of initial abundances, 
    # t_span: 1D array of starting/ending time steps (e.g. [0, num_hours]),
    # r: length-N 1D array of maximum growth rates,
    # alpha: NxN 2D array of competition coefficients (with diagonal set to 1),
    # delta: the global death rate, and
    # N: the number of species

    import numpy as np
    
    comp_sum = np.zeros(N)
    for k in range(N):
        for j in range(N):
            comp_sum[k] += alpha[k,j]*X[j]
            
    return_list = []
    
    for k in range(N):
        return_list.append(X[k]*(r[k]*(1-comp_sum[k])-delta))
        
    return return_list # return_list is the 1D array of species abundances




def Arr_100species_ivp(delta,mu,sigma,num_trials,p):
# This function calls the above function, dX_dt_GLV, to integrate the Lotka-Volterra
# model for 100 species over a range of temperatures and outputs the weighted mean
# copy number over that temperature range, along with the relative abundances of each
# copy number.
# It requires the following:
    # delta: the global death rate,
    # mu: the mean competition coefficient (alpha)
    # sigma: the standard deviation of the competition coefficients (alpha)
    # num_trials: the number of trials
    # p: geometric dist. parameter; the fraction of rRNA copy number = 1 in the starting community
    
    import numpy as np
    import random
    from scipy.integrate import solve_ivp
    

    T = [278,281,284,287,290,293,296,299] # Temperature range in Kelvin
    N = 100                               # Number of species
    num_hours = 800     # Number of hours per simulation-- decrease to decrease computational time
    t_span = [0,num_hours]

    max_cn = 10         # Maximum rRNA copy number
    cn_range = np.linspace(1,max_cn,max_cn)

    # Create geometric distribution of randomly drawn rRNA copy numbers:
    geom_dist = []
    for n in cn_range:
        geom_dist.append(p*(1-p)**n)
    geom_weights = geom_dist/np.sum(geom_dist)

    # Create empty arrays to hold results of weighted mean copy number and 
    # relative abundances across all trials and temperatures:
    WMCN = np.zeros([num_trials,len(T)])
    CN_rel_abunds = np.zeros([num_trials,max_cn,len(T)])

    for o in range(num_trials):

        # Randomly draw growth rates from uniform distribution:
        #R = np.random.randint(max_cn, size=N)+1

        # Randomly draw growth rates from geometric distribution:
        R = np.array(random.choices(cn_range, weights=geom_weights, k=N))       

        Gcoeff = 6500 # Arrhenius model parameter when r = aR*exp(-G/T); this is equivalent to r = aR*exp(-E/kT) but with different units
        G = np.ones(N)*Gcoeff

        Xo = 0.01 + np.zeros(N)

        r=[]
        a = 1e9
        for t in range(len(T)):
            r.append(a*R*np.exp(-G/T[t]))

        # Randomly draw competition coefficients from normal distribution:
        alpha = np.random.normal(mu,sigma,[N,N])
        for k in range(N):
            alpha[k,k]=1

        copynums = R

        for t in range(len(T)):
            
            # Integrate:
            Xs_ivp = solve_ivp(dX_dt_GLV, t_span, Xo, args=(r[t],alpha,delta,N))
            Xs = np.transpose(Xs_ivp.y)
        
            # Solve for relative abundances of each copy number:
            Biomass = np.sum(Xs[len(Xs)-1,:])
            for k in range(N):
                WMCN[o,t] += copynums[k]*(Xs[len(Xs)-1,k]/Biomass)
            
            rel_abunds = Xs[len(Xs)-1,:]/Biomass
            for c in range(max_cn):
                CN_rel_abunds[o][c,t] = np.sum(rel_abunds[R==(c+1)])
        
    return(T,WMCN,CN_rel_abunds)





def bootstrap_p(p_range,sample_size,num_trials,df,sim_df):
# This function fits the data to the LV model simulation results by bootstrapping
# the datasets and fitting them to the best value of p, the geometric distribution
# parameter.
# It requires the following:
    # p_range: the range of values of p used
    # sample_size: the bootstrapping sample size- should be equal to the dataset size
    # num_trials: number of trials
    # df: metadata dataframe containing 3 columns: Sample, WMCN, Temperature
    # sim_df: simulation result dataframe containing 3 columns: Temperature, p, WMCN
    
    import numpy as np
    import random
    
    WMCNs = df['WMCN'].tolist()
    temps = df['Temperature'].tolist()
    p_est = np.zeros([num_trials])
    min_error = np.zeros([num_trials])

    for o in range(num_trials):
        track_rms_error = np.zeros([len(p_range)])
        for k in range(sample_size):
            random_index = random.randrange(len(WMCNs))    # Pick a random WMCN
            random_temp = temps[random_index]              # Find corresponding temp
            random_temp_round = np.round(random_temp,2)    # Round temp to 2 decimal points
            random_WMCN = WMCNs[random_index]              # Match to WMCN
            this_df = sim_df.loc[random_temp_round]        # Call simulated WMCNs for that temp

            for j in range(len(p_range)):
                # Loop through p, calculate error between actual WMCN and simulated WMCN
                this_WMCN = this_df[this_df['p']==np.round(p_range[j],4)]['WMCN'].loc[random_temp_round]
                track_rms_error[j] += (random_WMCN - this_WMCN)**2

        rms_error = np.sqrt(track_rms_error/sample_size).tolist()
        p_est[o] = p_range[rms_error.index(min(rms_error))]
        min_error[o] = rms_error[rms_error.index(min(rms_error))]
        
    mean_min_error = np.mean(min_error)

    return(np.mean(p_est),np.std(p_est),mean_min_error) # Outputs the mean estimate of p, its SD, and a RMS error




def dX_dt_GenSR(t, X, r, delta, Co, num_species, num_resource):
# Linear consumer resource model, any number of species, resources
# This function can simulate any number of species and resources and outputs
# abundances of species and resources. It should be called like this: 
# sol = solve_ivp(dX_dt_GenSR, t_span, X_initial, args=(r,delta,Co,num_species,num_resource)), where:
    # X_initial: length-(num_species+num_resource) 1D array of initial species and resource abundances,
    # t_span: 1D array of starting/ending time steps (e.g. [0, num_hours]),
    # r: length-num_species 1D array of maximum growth rates,
    # delta: the global death rate,
    # Co: length-num_resource 1D array of initial resource abundances,
    # num_species: number of species
    # num_resource: number of resources

    import numpy as np

    C = np.zeros(num_resource)
    for k in range(num_resource):
        C[k] = X[k+num_species]
    
    rpc = np.zeros([num_species,num_resource])
    for k in range(num_species):
        for j in range(num_resource):
            rpc[k,j] = r[k,j]*C[j]
            
    rp_all_C = np.zeros([num_species])
    for k in range(num_species):
        rp_all_C[k] = np.sum(rpc[k,:])
        
    r_c = np.zeros([num_species,num_resource])
    for k in range(num_species):
        for j in range(num_resource):
            r_c[k,j] = r[k,j]*X[k]
            
    r_all_N = np.zeros([num_resource])
    for j in range(num_resource):
        r_all_N[j] = np.sum(r_c[:,j])
        
    return_list = []
    for k in range(num_species):
        return_list.append(X[k]*(rp_all_C[k] - delta))
        
    for j in range(num_resource):
        return_list.append(delta*(Co[j]-C[j]) - C[j]*r_all_N[j])
    
    return return_list # return_list is the 1D array of species abundances followed by resource abundances




def get_WMGR_ave(num_species,num_resource,num_hours,start_density,delta,num_trials,T,r_conc):
# This function calls the above function, dX_dt_GenSR, to integrate the linear consumer resource
# model for any number of species and resources over a range of temperatures and resource supply 
# concentrations and outputs the weighted mean maximum growth rate over that temperature range.
# It requires the following:
    # num_species: number of species
    # num_resources: number of resources
    # num_hours: number of hours per simulation
    # start_density: starting abundance of all species
    # delta: the global death rate,
    # num_trials: the number of trials
    # T: 1D array of temperature range in Kelvin
    # r_conc: 1D array of range of resource supply concentrations
    
    import numpy as np
    from scipy.integrate import solve_ivp
    
    t_span = [0,num_hours]
    No = start_density + np.zeros(num_species+num_resource)
    T_C = [x - 273 for x in T] # Convert Kelvin to Celsius
    # Create empty array to hold weighted mean maximum growth rates for all simulations:
    WMGR_rm = np.zeros([num_trials, len(r_conc), len(T)])
    
    for o in range(num_trials):

        Gcoeff = 3864   # Arrhenius parameter
        G = np.ones([num_species,num_resource])*Gcoeff
        # Draw the growth rates for each species on each resource from a uniform distribution:
        rcoeff = np.random.uniform(.5e5,3e5,[num_species,num_resource])

        # Set the growth rate at which the WMGR will be "measured":
        rm = rcoeff*np.exp(-G/T[len(T)-1])/num_resource

        for c in range(len(r_conc)):
            # Set the resource supply concentration:
            Co = r_conc[c]  + np.zeros(num_resource)
            
            # Set the initial resource abundances equal to the resource supply concentration:
            for j in range(num_resource):
                No[j+num_species] = Co[0]

            for h in range(len(T)):
                r = rcoeff*np.exp(-G/T[h])/num_resource
                
                # Integrate:
                sol = solve_ivp(dX_dt_GenSR, t_span, No, args=(r,delta,Co,num_species,num_resource))
                # Calculate WMGRs:
                num_steps = len(sol.t)
                Biomass = np.sum(sol.y[0:num_species,num_steps-1])
                WMGR_rm[o,c,h] = np.sum(np.sum(rm,1)*(sol.y[0:num_species,num_steps-1]/Biomass))

    WMGR_rm_ave = np.mean(WMGR_rm,0)
    
    # Create tick marks for easy heatmap plotting:
    xticks=np.arange(-1-1/len(T_C),1+1/len(T_C),2/(len(T_C)))
    xticks=xticks[1:]
    yticks=np.arange(-1-1/len(r_conc),1+1/len(r_conc),2/(len(r_conc)))
    yticks=yticks[1:]
    
    return [WMGR_rm_ave,xticks,yticks]

