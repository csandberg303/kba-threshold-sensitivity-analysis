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


# In[2]:


# FORMULA TO GET LANA'S INPUT FILES FROM THE REPO TO LOCAL DIRECTORY
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
        fileinfo_df = pd.read_csv(io.StringIO(fileinfo.decode('utf-8')),
                                 index_col=False).squeeze("columns")
        filename = file
        output = fileinfo_df.to_csv(file, index=False)
     return output




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


# In[3]:


# THESE FUNCTIONS ARE TAKEN/ADAPTED FROM QMARXAN TOOLBOX ALGORITHM CODE

# taken from lines ~98-104 of QMarxan algorithm

# formatAsME - format as Marxan Exponent format like
        #              Input File Editor
        #
def formatAsME(inVal):
    outStr = "%.14E" % float(inVal)
    parts = outStr.split('E')
    sign = parts[1][:1]
    exponent = "%04d" % float(parts[1][1:])
    outStr = parts[0] + 'E' +  sign + exponent
    return(outStr)



# TO CREATE INPUT.DAT (LINES 128-183 OF ALGORITHM FILE)
def create_input_dat(dest):
    """
     To create the input.dat file that stores processing parameters

     Parameters
     ----------
     dest : str
     directory input.dat file will be saved to

     other parameters will be added to replace the default initial values
     that are included in the QMarxan code

     -------
     returned_data : the input.dat file

    """
    output = os.path.join(dest,'qm_input.dat')
    f = open(output, 'w')
    f.write("Input file for Annealing program.\n")
    f.write('\n')
    f.write('This file generated for KBA Threshold\n')
    f.write('Analysis project using code from\n')
    f.write('QMarxan Toolbox 2.0\n')
    f.write('created by Apropos Information Systems Inc.\n')
    f.write('\n')
    f.write("General Parameters\n")
    f.write("PROP %s\n" % formatAsME(0.5))
    f.write("RANDSEED -1\n")
    f.write("NUMREPS 100\n")
    f.write('\n')
    f.write("Annealing Parameters\n")
    f.write("NUMITNS 1000000\n")
    f.write("STARTTEMP %s\n" % formatAsME(-1.0))
    f.write("COOLFAC %s\n" % formatAsME(-1.0))
    f.write("NUMTEMP 10000\n")
    f.write("")
    f.write("Cost Threshold\n")
    f.write("COSTTHRESH %s\n" % formatAsME(0.0))
    f.write("THRESHPEN1 %s\n" % formatAsME(0.0))
    f.write("THRESHPEN2 %s\n" % formatAsME(0.0))
    f.write("\n")
    f.write("Input Files\n")
    f.write("INPUTDIR input\n")
    f.write("SPECNAME spec.dat\n")
    f.write("PUNAME pu.dat\n")
    f.write("PUVSPRNAME puvsp.dat\n")
    f.write("BOUNDNAME bound.dat\n")
    f.write("MATRIXSPORDERNAME puvsp_sporder.dat\n")
    f.write("\n")
    f.write("Save Files\n")
    f.write("SCENNAME output\n")
    f.write("SAVERUN 3\n")
    f.write("SAVEBEST 3\n")
    f.write("SAVESUMMARY 3\n")
    f.write("SAVESCEN 3\n")
    f.write("SAVETARGMET 3\n")
    f.write("SAVESUMSOLN 3\n")
    f.write("SAVELOG 3\n")
    f.write("SAVESNAPSTEPS 0\n")
    f.write("SAVESNAPCHANGES 0\n")
    f.write("SAVESNAPFREQUENCY 0\n")
    f.write("OUTPUTDIR output\n")
    f.write("\n")
    f.write("Program control\n")
    f.write("RUNMODE 1\n")
    f.write("MISSLEVEL %s\n" % formatAsME(0.95))
    f.write("ITIMPTYPE 1\n")
    f.write("HEURTYPE -1\n")
    f.write("CLUMPTYPE 0\n")
    f.write("VERBOSITY 2\n")
    f.write("SAVESOLUTIONSMATRIX 3\n")
    f.write("\n")
    f.close()
    return output


# In[ ]:
