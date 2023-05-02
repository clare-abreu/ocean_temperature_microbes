#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
These functions assign copy numbers to taxa at the highest-resolution taxonomic 
level available from the rrnDB, using the mean from 
rrnDB-5.6_pantaxa_stats_NCBI_WITH_SAR11_CLADES.csv. Most datasets 
do not contain species classifications, but the first function will handle 
species if they are included. The second function does not, but is 
otherwise identical.

Written by Clare Abreu for Abreu & Dal Bello et al, Science Advances (2023)
"""


def assign_copynums_with_species(data):
    # This function assigns copy numbers to taxa at the highest-resolution taxanomic
    # level available from the rrnDB, using the mean from that table. Most datasets
    # do not contain species classifications, but this function will handle species
    # if they are included. The following function does not, but is otherwise identical.
    
    import numpy as np
    import pandas as pd

    # Remove extra spaces:
    data['Phylum'] = data['Phylum'].str.strip()
    data['Class'] = data['Class'].str.strip()
    data['Order'] = data['Order'].str.strip()
    data['Family'] = data['Family'].str.strip()
    data['Genus'] = data['Genus'].str.strip()
    data['Species'] = data['Species'].str.strip()
    
    NCBI_copy_nums = pd.read_csv('../data/rrnDB-5.6_pantaxa_stats_NCBI_WITH_SAR11_CLADES.csv')
    NCBI_copy_nums.set_index('taxid', inplace=True)
    species_nums = NCBI_copy_nums['rank'] == 'species'
    genus_nums = NCBI_copy_nums['rank'] == 'genus'
    family_nums = NCBI_copy_nums['rank'] == 'family'
    order_nums = NCBI_copy_nums['rank'] == 'order'
    class_nums = NCBI_copy_nums['rank'] == 'class'
    phylum_nums = NCBI_copy_nums['rank'] == 'phylum'
    NCBI_species = NCBI_copy_nums[genus_nums]
    NCBI_genus = NCBI_copy_nums[genus_nums]
    NCBI_family = NCBI_copy_nums[family_nums]
    NCBI_order = NCBI_copy_nums[order_nums]
    NCBI_class = NCBI_copy_nums[class_nums]
    NCBI_phylum = NCBI_copy_nums[phylum_nums]
    
    data['Copy Number'] = np.nan
    #data['#OTU_ID'] = np.nan
    data['Copy Number Classification Level'] = np.nan
    species_list = data['Species']
    genus_list = data['Genus']
    family_list = data['Family']
    order_list = data['Order']
    class_list = data['Class']
    phylum_list = data['Phylum']
    to_sequence = [] # List of OTU nums that don't have copy num match
    unknown_levels = [] # List of levels that don't have copy num match
    sequences = [] # List of OTU 16S sequences that don't have copy num match

    for k in range(len(data)):
        OTU_ID = data['OTU ID'].iloc[k]
        
        if not pd.isnull(species_list.iloc[k]):
            that_row = NCBI_species[NCBI_species['name']==species_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Species'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_genus[NCBI_genus['name']==genus_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Genus'
                    else:
                        print('Warning: Two entries found')
                else:
                    that_row = NCBI_family[NCBI_family['name']==family_list.iloc[k]]
                    if len(that_row)>0:
                        if len(that_row) == 1:
                            data['Copy Number'].iloc[k] = float(that_row['mean'])
                            data['Copy Number Classification Level'].iloc[k] = 'Family'
                        else:
                            print('Warning: Two entries found')
                    else:
                        that_row = NCBI_order[NCBI_order['name']==order_list.iloc[k]]
                        if len(that_row)>0:
                            if len(that_row) == 1:
                                data['Copy Number'].iloc[k] = float(that_row['mean'])
                                data['Copy Number Classification Level'].iloc[k] = 'Order'
                            else:
                                print('Warning: Two entries found')
                        else:
                            that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
                            if len(that_row)>0:
                                if len(that_row) == 1:
                                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                                    data['Copy Number Classification Level'].iloc[k] = 'Class'
                                else:
                                    print('Warning: Two entries found')
                                    
                            else:
                                that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                                if len(that_row)>0:
                                    if len(that_row) == 1:
                                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                                        data['Copy Number Classification Level'].iloc[k] = 'Class'
                                    else:
                                        print('Warning: Two entries found')
                                    
                                else:
                                    #print('Error: No classification found')
                                    data['Copy Number'].iloc[k] = 'nan'
                                    uknown_OTU = OTU_ID
                                    to_sequence.append(uknown_OTU)
                                    unknown_levels.append(1)
                                    #sequences.append(data['Sequence'][k])
        
        
    
        elif not pd.isnull(genus_list.iloc[k]):
            that_row = NCBI_genus[NCBI_genus['name']==genus_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Genus'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_family[NCBI_family['name']==family_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Family'
                    else:
                        print('Warning: Two entries found')
                else:
                    that_row = NCBI_order[NCBI_order['name']==order_list.iloc[k]]
                    if len(that_row)>0:
                        if len(that_row) == 1:
                            data['Copy Number'].iloc[k] = float(that_row['mean'])
                            data['Copy Number Classification Level'].iloc[k] = 'Order'
                        else:
                            print('Warning: Two entries found')
                    else:
                        that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
                        if len(that_row)>0:
                            if len(that_row) == 1:
                                data['Copy Number'].iloc[k] = float(that_row['mean'])
                                data['Copy Number Classification Level'].iloc[k] = 'Class'
                            else:
                                print('Warning: Two entries found')
                        else:
                            that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                            if len(that_row)>0:
                                if len(that_row) == 1:
                                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                                    data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                                else:
                                    print('Warning: Two entries found')
                            else:
                                #print('Error: No classification found')
                                data['Copy Number'].iloc[k] = 'nan'
                                uknown_OTU = OTU_ID
                                to_sequence.append(uknown_OTU)
                                unknown_levels.append(1)
                                #sequences.append(data['Sequence'][k])
                
        
        elif not pd.isnull(family_list.iloc[k]):
            that_row = NCBI_family[NCBI_family['name']==family_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Family'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_order[NCBI_order['name']==order_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Order'
                    else:
                        print('Warning: Two entries found')
                else:
                    that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
                    if len(that_row)>0:
                        if len(that_row) == 1:
                            data['Copy Number'].iloc[k] = float(that_row['mean'])
                            data['Copy Number Classification Level'].iloc[k] = 'Class'
                        else:
                            print('Warning: Two entries found')
                    else:
                        that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                        if len(that_row)>0:
                            if len(that_row) == 1:
                                data['Copy Number'].iloc[k] = float(that_row['mean'])
                                data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                            else:
                                print('Warning: Two entries found')
                        else:
                            #print('Error: No classification found')
                            data['Copy Number'].iloc[k] = 'nan'
                            uknown_OTU = OTU_ID
                            to_sequence.append(uknown_OTU)
                            unknown_levels.append(2)
                            #sequences.append(data['Sequence'][k])
        
        elif not pd.isnull(order_list.iloc[k]):
            that_row = NCBI_order[NCBI_order['name']==order_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Order'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Class'
                    else:
                        print('Warning: Two entries found')
                else:
                    that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                    if len(that_row)>0:
                        if len(that_row) == 1:
                            data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                            data['Copy Number'].iloc[k] = float(that_row['mean'])
                        else:
                            print('Warning: Two entries found')
                    else:
                        #print('Error: No classification found')
                        data['Copy Number'].iloc[k] = 'nan'
                        uknown_OTU = OTU_ID
                        to_sequence.append(uknown_OTU)
                        unknown_levels.append(3)
                        #sequences.append(data['Sequence'][k])
        
        
        elif not pd.isnull(class_list.iloc[k]):
            that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Class'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                    else:
                        print('Warning: Two entries found')
                else:
                    #print('Error: No classification found')
                    data['Copy Number'].iloc[k] = 'nan'
                    uknown_OTU = OTU_ID
                    to_sequence.append(uknown_OTU)
                    unknown_levels.append(4)
                    #sequences.append(data['Sequence'][k])

        elif not pd.isnull(phylum_list.iloc[k]):
            that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                else:
                    print('Warning: Two entries found')
            else:
                #print('Error: No classification found')
                data['Copy Number'].iloc[k] = 'nan'
                uknown_OTU = OTU_ID
                to_sequence.append(uknown_OTU)
                unknown_levels.append(5)
                #sequences.append(data['Sequence'][k])

        else:
            #print('Error: No classification found')
            data['Copy Number'].iloc[k] = 'nan'
            uknown_OTU = OTU_ID
            to_sequence.append(uknown_OTU)
            unknown_levels.append(6)
            #sequences.append(data['Sequence'][k])


    return(data,to_sequence,unknown_levels)












def assign_copynums(data):
    # This function is identical to the previous one, assign_copynums_with_species,
    # with the one difference that it does not handle species classifications.
    
    import numpy as np
    import pandas as pd

    # Remove extra spaces:
    data['Phylum'] = data['Phylum'].str.strip()
    data['Class'] = data['Class'].str.strip()
    data['Order'] = data['Order'].str.strip()
    data['Family'] = data['Family'].str.strip()
    data['Genus'] = data['Genus'].str.strip()
    #data['Species'] = data['Species'].str.strip()
    
    NCBI_copy_nums = pd.read_csv('../data/rrnDB-5.6_pantaxa_stats_NCBI_WITH_SAR11_CLADES.csv')
    NCBI_copy_nums.set_index('taxid', inplace=True)
    species_nums = NCBI_copy_nums['rank'] == 'species'
    genus_nums = NCBI_copy_nums['rank'] == 'genus'
    family_nums = NCBI_copy_nums['rank'] == 'family'
    order_nums = NCBI_copy_nums['rank'] == 'order'
    class_nums = NCBI_copy_nums['rank'] == 'class'
    phylum_nums = NCBI_copy_nums['rank'] == 'phylum'
    NCBI_genus = NCBI_copy_nums[genus_nums]
    NCBI_family = NCBI_copy_nums[family_nums]
    NCBI_order = NCBI_copy_nums[order_nums]
    NCBI_class = NCBI_copy_nums[class_nums]
    NCBI_phylum = NCBI_copy_nums[phylum_nums]
    
    data['Copy Number'] = np.nan
    #data['#OTU_ID'] = np.nan
    data['Copy Number Classification Level'] = np.nan
    genus_list = data['Genus']
    family_list = data['Family']
    order_list = data['Order']
    class_list = data['Class']
    phylum_list = data['Phylum']
    to_sequence = [] # List of OTU nums that don't have copy num match
    unknown_levels = [] # List of levels that don't have copy num match
    sequences = [] # List of OTU 16S sequences that don't have copy num match

    for k in range(len(data)):
        OTU_ID = data['OTU ID'].iloc[k]
    
        if not pd.isnull(genus_list.iloc[k]):
            that_row = NCBI_genus[NCBI_genus['name']==genus_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Genus'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_family[NCBI_family['name']==family_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Family'
                    else:
                        print('Warning: Two entries found')
                else:
                    that_row = NCBI_order[NCBI_order['name']==order_list.iloc[k]]
                    if len(that_row)>0:
                        if len(that_row) == 1:
                            data['Copy Number'].iloc[k] = float(that_row['mean'])
                            data['Copy Number Classification Level'].iloc[k] = 'Order'
                        else:
                            print('Warning: Two entries found')
                    else:
                        that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
                        if len(that_row)>0:
                            if len(that_row) == 1:
                                data['Copy Number'].iloc[k] = float(that_row['mean'])
                                data['Copy Number Classification Level'].iloc[k] = 'Class'
                            else:
                                print('Warning: Two entries found')
                        else:
                            that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                            if len(that_row)>0:
                                if len(that_row) == 1:
                                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                                    data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                                else:
                                    print('Warning: Two entries found')
                            else:
                                #print('Error: No classification found')
                                data['Copy Number'].iloc[k] = 'nan'
                                uknown_OTU = OTU_ID
                                to_sequence.append(uknown_OTU)
                                unknown_levels.append(1)
                                #sequences.append(data['Sequence'][k])
                
        
        elif not pd.isnull(family_list.iloc[k]):
            that_row = NCBI_family[NCBI_family['name']==family_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Family'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_order[NCBI_order['name']==order_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Order'
                    else:
                        print('Warning: Two entries found')
                else:
                    that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
                    if len(that_row)>0:
                        if len(that_row) == 1:
                            data['Copy Number'].iloc[k] = float(that_row['mean'])
                            data['Copy Number Classification Level'].iloc[k] = 'Class'
                        else:
                            print('Warning: Two entries found')
                    else:
                        that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                        if len(that_row)>0:
                            if len(that_row) == 1:
                                data['Copy Number'].iloc[k] = float(that_row['mean'])
                                data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                            else:
                                print('Warning: Two entries found')
                        else:
                            #print('Error: No classification found')
                            data['Copy Number'].iloc[k] = 'nan'
                            uknown_OTU = OTU_ID
                            to_sequence.append(uknown_OTU)
                            unknown_levels.append(2)
                            #sequences.append(data['Sequence'][k])
        
        elif not pd.isnull(order_list.iloc[k]):
            that_row = NCBI_order[NCBI_order['name']==order_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Order'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Class'
                    else:
                        print('Warning: Two entries found')
                else:
                    that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                    if len(that_row)>0:
                        if len(that_row) == 1:
                            data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                            data['Copy Number'].iloc[k] = float(that_row['mean'])
                        else:
                            print('Warning: Two entries found')
                    else:
                        #print('Error: No classification found')
                        data['Copy Number'].iloc[k] = 'nan'
                        uknown_OTU = OTU_ID
                        to_sequence.append(uknown_OTU)
                        unknown_levels.append(3)
                        #sequences.append(data['Sequence'][k])
        
        
        elif not pd.isnull(class_list.iloc[k]):
            that_row = NCBI_class[NCBI_class['name']==class_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Class'
                else:
                    print('Warning: Two entries found')
            else:
                that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
                if len(that_row)>0:
                    if len(that_row) == 1:
                        data['Copy Number'].iloc[k] = float(that_row['mean'])
                        data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                    else:
                        print('Warning: Two entries found')
                else:
                    #print('Error: No classification found')
                    data['Copy Number'].iloc[k] = 'nan'
                    uknown_OTU = OTU_ID
                    to_sequence.append(uknown_OTU)
                    unknown_levels.append(4)
                    #sequences.append(data['Sequence'][k])

        elif not pd.isnull(phylum_list.iloc[k]):
            that_row = NCBI_phylum[NCBI_phylum['name']==phylum_list.iloc[k]]
            if len(that_row)>0:
                if len(that_row) == 1:
                    data['Copy Number'].iloc[k] = float(that_row['mean'])
                    data['Copy Number Classification Level'].iloc[k] = 'Phylum'
                else:
                    print('Warning: Two entries found')
            else:
                #print('Error: No classification found')
                data['Copy Number'].iloc[k] = 'nan'
                uknown_OTU = OTU_ID
                to_sequence.append(uknown_OTU)
                unknown_levels.append(5)
                #sequences.append(data['Sequence'][k])

        else:
            #print('Error: No classification found')
            data['Copy Number'].iloc[k] = 'nan'
            uknown_OTU = OTU_ID
            to_sequence.append(uknown_OTU)
            unknown_levels.append(6)
            #sequences.append(data['Sequence'][k])


    return(data,to_sequence,unknown_levels)
