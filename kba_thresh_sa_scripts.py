#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import libraries
import os
import datetime
import requests
import io
from glob import glob

import earthpy as et
import numpy as np
import pandas as pd

# create function to write targets.csv files, for each threshold test value

# NOT SURE IF THIS WILL ULTIMATELY BE NEEDED IN ITS CURRENT FORM, AS THE 
#TARGETS INPUT FILE MIGHT ONLY BE USED BY THE QGIS ADDIN CLUZ, RATHER THAN 
# MARXAN/MARXANCONPY ITSELF

def create_targets_files(eco, thresholds_test, eco_info, path):
 # YOUR CODE HERE
     """creates the targets.csv files needed for Marxan analysis 
     (?? only when using CLUZ add-in in QGIS ??).

     Parameters
     ----------
     eco : str
     name of ecosystem that will be analyzed by Marxan

     thresholds_test : list
     list of threshold values to be tested for each ecosystem

     eco_info : dataframe
     source of info for each ecosystem, with columns 'OID' (Unique ID number),
     'Name' (ecosystem name), Type (number representing RLE Status), Size of 
     Ecosystem (units of area measurement) and the Current IUCN Threshold 
     value, based upon ecosystem's RLE status

     path : filepath
     filepath to ecosystem subdirectory where targets files will be saved

     Returns
     -------
     returned_data : csv
     csv files are saved to ecosystem directories, one file for each threshold
     value to be tested
     """
     for val in thresholds_test:
            target_info = {'Id': [eco_info.loc[eco]['OID']], 
                           'Name': [eco], 
                           'Type': [eco_info.loc[eco]['Type']], 
                           'sq_km': [eco_info.loc[eco]['US_km2']],
                           'iucn_th': [eco_info.loc[eco]['Current_IUCN_TH']]}
            target_df = pd.DataFrame(data=target_info).set_index('Id')
            target_df['Target'] = (target_df['sq_km'] * target_df['iucn_th'])
            target_df['Target'] = (val * target_df['Target'])
            target_df.drop(["sq_km", "iucn_th"], axis = 1, inplace = True)
            outpath = os.path.join(path, 'targets_' +str(val) + '.csv')
            target_df.to_csv(outpath)
     return 
   

