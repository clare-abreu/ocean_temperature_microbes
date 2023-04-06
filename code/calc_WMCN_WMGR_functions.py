#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 22:35:56 2023

@author: clare
"""
def calc_wmcn_only(this_data,data_path):
    # This function calculates WMCN for each sample in a dataset. The datasets
    # used for this function are in the Generic Data folder of data/.
    
    import numpy as np
    import pandas as pd
    
    df_tax = pd.read_csv(f'{data_path}{this_data}_generic_data.csv', index_col=0)
    df_meta = pd.read_csv(f'{data_path}{this_data}_generic_metadata.csv', index_col=0)

    All_samples = df_meta.index.tolist()
    All_OTUs = df_tax.index.tolist()

    df_meta['WMCN'] = np.nan
    df_meta['WMCN, no SAR11'] = np.nan
    df_meta['WMCN, no CN1'] = np.nan
    df_meta['WMCN, Heterotrophs'] = np.nan

    #Exclude SAR11_clade from Order:
    df_tax_nosar11 = df_tax[df_tax['Order']!='SAR11_clade'].copy()
    df_tax_nosar11 = df_tax_nosar11[df_tax_nosar11['Order']!='SAR11 clade'].copy()
    df_tax_nosar11 = df_tax_nosar11[df_tax_nosar11['Order']!='Pelagibacterales'].copy()
    df_tax_noCN1 = df_tax[df_tax['Copy Number']>1].copy()
    df_tax_hetero = df_tax[df_tax['Phototroph']==0].copy()

    # Compute WMCNs
    df_tax_c = df_tax[~df_tax['Copy Number'].isna()]
    df_tax_nosar11_c = df_tax_nosar11[~df_tax_nosar11['Copy Number'].isna()]
    df_tax_noCN1_c = df_tax_noCN1[~df_tax_noCN1['Copy Number'].isna()]
    df_tax_hetero_c = df_tax_hetero[~df_tax_hetero['Copy Number'].isna()]
    
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN = np.sum(df_tax_c[this_sample])/np.sum(df_tax_c[this_sample]/df_tax_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN'] = this_WMCN

    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN_nosar11 = np.sum(df_tax_nosar11_c[this_sample])/np.sum(df_tax_nosar11_c[this_sample]/df_tax_nosar11_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN, no SAR11'] = this_WMCN_nosar11
            
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN_noCN1 = np.sum(df_tax_noCN1_c[this_sample])/np.sum(df_tax_noCN1_c[this_sample]/df_tax_noCN1_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN, no CN1'] = this_WMCN_noCN1
            
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN_hetero = np.sum(df_tax_hetero_c[this_sample])/np.sum(df_tax_hetero_c[this_sample]/df_tax_hetero_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN, Heterotrophs'] = this_WMCN_hetero
            
    return(df_meta)





def calc_wmgr_wmcn(this_data,data_path):
    
    # This function calculates both WMCN and WMGR for each sample in a dataset. The datasets
    # used for this function are in the Generic Data folder of data/.
    
    import numpy as np
    import pandas as pd
    
    if this_data == 'SPT':
        df_tax = pd.read_csv(f'{data_path}SPT_generic_data.csv', index_col=0, dtype={'Copy Number Classification Level': str})
        df_tax['Copy Number Classification Level'] = df_tax['Copy Number Classification Level'].fillna("")
        df_meta=pd.read_csv(f'{data_path}SPT_generic_metadata.csv', dtype={'Unnamed: 0': str})
        df_meta.rename({'Unnamed: 0': 'Sample'}, axis=1, inplace=True)
        df_meta.set_index('Sample', inplace=True)
    else:
        df_tax = pd.read_csv(f'{data_path}{this_data}_generic_data.csv', index_col=0)
        df_meta = pd.read_csv(f'{data_path}{this_data}_generic_metadata.csv', index_col=0)

    All_samples = df_meta.index.tolist()
    All_OTUs = df_tax.index.tolist()

    df_meta['WMGR'] = np.nan
    df_meta['WMGR, no SAR11'] = np.nan
    df_meta['WMGR, copio'] = np.nan
    df_meta['WMGR, Heterotrophs'] = np.nan
    df_meta['WMCN'] = np.nan
    df_meta['WMCN, no SAR11'] = np.nan
    df_meta['WMCN, no CN1'] = np.nan
    df_meta['WMCN, Heterotrophs'] = np.nan

    #Exclude SAR11_clade from Order:
    df_tax_nosar11 = df_tax[df_tax['Order']!='SAR11_clade'].copy()
    df_tax_nosar11 = df_tax_nosar11[df_tax_nosar11['Order']!='SAR11 clade'].copy()
    df_tax_nosar11 = df_tax_nosar11[df_tax_nosar11['Order']!='Pelagibacterales'].copy()
    df_tax_noCN1 = df_tax[df_tax['Copy Number']>1].copy()
    df_tax_copio = df_tax[df_tax['Growth Rate']>0.1386].copy()
    df_tax_hetero = df_tax[df_tax['Phototroph']==0].copy()

    # Compute WMGRs
    df_tax_r = df_tax[~df_tax['Growth Rate'].isna()]
    df_tax_nosar11_r = df_tax_nosar11[~df_tax_nosar11['Growth Rate'].isna()]
    df_tax_copio_r = df_tax_copio[~df_tax_copio['Growth Rate'].isna()]
    df_tax_hetero_r = df_tax_hetero[~df_tax_hetero['Growth Rate'].isna()]
    
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax_r.columns.tolist():
            this_WMGR = np.sum(df_tax_r[this_sample]*df_tax_r['Growth Rate'])/np.sum(df_tax_r[this_sample])
            df_meta.loc[this_sample,'WMGR'] = this_WMGR

    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax_r.columns.tolist():
            this_WMGR_nosar11 = np.sum(df_tax_nosar11_r[this_sample]*df_tax_nosar11_r['Growth Rate'])/np.sum(df_tax_nosar11_r[this_sample])
            df_meta.loc[this_sample,'WMGR, no SAR11'] = this_WMGR_nosar11

    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax_r.columns.tolist():
            this_WMGR_copio = np.sum(df_tax_copio_r[this_sample]*df_tax_copio_r['Growth Rate'])/np.sum(df_tax_copio_r[this_sample])
            df_meta.loc[this_sample,'WMGR, copio'] = this_WMGR_copio
            
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax_r.columns.tolist():
            this_WMGR_hetero = np.sum(df_tax_hetero_r[this_sample]*df_tax_hetero_r['Growth Rate'])/np.sum(df_tax_hetero_r[this_sample])
            df_meta.loc[this_sample,'WMGR, Heterotrophs'] = this_WMGR_hetero

    # Compute WMCNs
    df_tax_c = df_tax[~df_tax['Copy Number'].isna()]
    df_tax_nosar11_c = df_tax_nosar11[~df_tax_nosar11['Copy Number'].isna()]
    df_tax_noCN1_c = df_tax_noCN1[~df_tax_noCN1['Copy Number'].isna()]
    df_tax_hetero_c = df_tax_hetero[~df_tax_hetero['Copy Number'].isna()]
    
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN = np.sum(df_tax_c[this_sample])/np.sum(df_tax_c[this_sample]/df_tax_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN'] = this_WMCN

    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN_nosar11 = np.sum(df_tax_nosar11_c[this_sample])/np.sum(df_tax_nosar11_c[this_sample]/df_tax_nosar11_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN, no SAR11'] = this_WMCN_nosar11
            
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN_noCN1 = np.sum(df_tax_noCN1_c[this_sample])/np.sum(df_tax_noCN1_c[this_sample]/df_tax_noCN1_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN, no CN1'] = this_WMCN_noCN1
            
    for k in range(len(All_samples)):
        this_sample = All_samples[k]
        if this_sample in df_tax.columns.tolist():
            this_WMCN_hetero = np.sum(df_tax_hetero_c[this_sample])/np.sum(df_tax_hetero_c[this_sample]/df_tax_hetero_c['Copy Number'])
            df_meta.loc[this_sample,'WMCN, Heterotrophs'] = this_WMCN_hetero
            
    return(df_meta)