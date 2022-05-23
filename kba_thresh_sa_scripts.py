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

def create_cluz_targets_files(eco, thresholds_test, eco_info, path):
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
            csv = target_df.to_csv(outpath)
     return csv
   
    
def get_marxan_input_files(eco, files_to_get):
     """
     Currently this formula will find the input files Lana created using the
     ArcMarxan Toolbox plugin in ArcGIS, which have been stored to the assets
     directory of our GitHub repository.  We hope this may be a placeholder
     function, to be replaced with functions that might create these files 
     directly using the opensource code available from the opensource QMarxan 
     Toolbox plugin for QGIS.

     Parameters
     ----------
     eco : str
     the abbreviated one word short name used for ecosystem, identifies a
     subdirectory of the marxan_input directory
     
     files_to_get : list
     list of filenames to retrieve from the marxan_input/eco directory of 
     the repo

     -------
     returned_data : the specified dat files, saved to eco/input local 
     directory
     """
     inputfile_ls = files_to_get
     
     for file in inputfile_ls:
        urltext = "https://raw.githubusercontent.com/csandberg303/kba-threshold-sensitivity-analysis/main/assets/data/marxan_input/"
        url = urltext + eco + "/" + file
        # downloading the info from file stored on github
        fileinfo = requests.get(url).content
        # Reading the downloaded content and turning it to a pandas dataframe
        fileinfo_df = pd.read_csv(io.StringIO(fileinfo.decode('utf-8')))
        filename = file
        output = fileinfo_df.to_csv(file)
     return output
   


# In[ ]:




