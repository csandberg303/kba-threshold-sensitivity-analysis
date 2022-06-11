#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import libraries
import os
import datetime
import io
import requests
import shutil
from glob import glob

# from qgis.core import *

import earthpy as et
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rasterio.crs import CRS
from rasterio.plot import plotting_extent
import rioxarray as rxr


# In[2]:


# create function to get Lana's input files (created with ArcGIS) from the
# repo to local directory
def get_marxan_input_files(eco, files_to_get, scen_id):
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
     the abbreviated one word short name used for ecosystem being analyzed;
     identifies a subdirectory of the timestamped marxan run directory

     files_to_get : list
     list of filenames to retrieve from the marxan_input/eco directory of
     the repo

     -------
     returned_data : the specified dat files, saved to eco/input local
     directory
     """
     inputfile_ls = files_to_get

     for file in inputfile_ls:
        urltext = ("https://raw.githubusercontent.com/csandberg303/"
                   "kba-threshold-sensitivity-analysis/main/assets/data/"
                   "marxan_input/")
        url = urltext + eco + "/" + file
        # downloading the info from file stored on github
        fileinfo = requests.get(url).content
        # Reading the downloaded content and turning it to a pandas dataframe
        fileinfo_df = pd.read_csv(io.StringIO(fileinfo.decode('utf-8')),
                                 index_col=False).squeeze("columns")
        filename = file
        output = fileinfo_df.to_csv(file, index=False)
        print(scen_id + ": " + file + " successfully copied from url")
     return output


# In[3]:


# write formula to get the shapefile and rasters that have been saved to the
#''kba_thres_sa/shp_hex' and 'kba_thres_sa/r_tif' local directories

# Currently I've manually copied Lana's ArcGIS files to these locations, using
# the naming convention 'eco.shp' for hexfiles and 'eco.tif for the rasters.

# IN THE FUTURE, the .shp & .tif files may be created and placed in the
# 'shp_hex' and 'r_tif' directories using code rather than ArcGIS, but this
# 'get_source_files' formula will still function to copy the needed files into
# the 'eco' directory when the 'eco/input' directories are created.
def get_source_files(path, eco, scen_id):
    """
    path : str
    local directory where the shapefiles or rasters are stored

    eco : str
    the abbreviated one word short name used for ecosystem being analyzed;
    identifies a subdirectory of the timestamped marxan run directory

    scen_id : str
    scenario id, info to be included as prefix on generated output

    """
    source_file_ls = glob(os.path.join(path, eco + '*'))
    if source_file_ls == []:
        print("no files found in " + path + "with expected name " + eco + "?")
    else:
        for file in source_file_ls:
            shutil.copy(file, os.getcwd())
            print(scen_id + ": "+ os.path.basename(file) +
                  " copied successfully")

    return print(scen_id + ": finished copying source files from " + path)


# In[4]:


# set crs of shp and tif to ESPG 5070 and save as new files

def set_source_files_crs (path, eco, espg, scen_id):
    """
    To set crs of shp and tif to ESPG 5070, add columns to shp and save as new
    files

    Parameters
    ----------
    path : str
    path to local 'kba_thresh_sa' directory where 'hex_shp' and 'r_tif'
    directories are stored

    eco : str
    the abbreviated one word short name used for ecosystem being analyzed;
    identifies a subdirectory of the timestamped marxan run directory

    espg : str
    espg number (we're using ESPG:5070)

    scen_id : str
    scenario id, info to be included as prefix on generated output

    -------
    returned_data : updated shp and tif


    """
    # open the shp and tif files saved at given path
    shp_data_path = os.path.join(path, "hex_shp", eco + '.shp')
    shp_layer = gpd.read_file(shp_data_path)

    tif_data_path = os.path.join(path, "r_tif", eco + '.tif')
    tif_layer = rxr.open_rasterio(tif_data_path, masked=True).squeeze()

    # reproject CRS of shp
    shp_layer_crs = shp_layer.to_crs(epsg=espg)

    # create shp output variable
    shp_espg_file = shp_layer_crs.to_file(eco + "_espg_" + espg + ".shp",
                                          index=False)

    # reproject CRS of tif
    # create a rasterio crs object
    crs_espg = CRS.from_string('EPSG:' + espg)
    # reproject tif using the crs object
    tif_layer_crs = tif_layer.rio.reproject(crs_espg)

    # create path that new tif file will be saved to
    tif_layer_crs_path = os.path.join(os.getcwd(),
                                      eco + "_espg_" + espg + ".tif")
    # create tif output variable
    tif_espg_file = tif_layer_crs.rio.to_raster(tif_layer_crs_path)

    # save the reprojected .shp and .tif files
    output = (shp_espg_file, tif_espg_file)

    print(scen_id + ": finished set_source_files_crs")
    return output


# In[5]:


# function to create pu.dat file

def create_pu_dat(eco, path, scen_id):
    """
    To create the pu.dat file that stores information about planning units in
    hex grid

    Parameters
    ----------
    eco : str
    the abbreviated one word short name used for ecosystem being analyzed;
    identifies a subdirectory of the timestamped marxan run directory

    path : str
    local directory where 'hex_shp' directory is stored

    scen_id : str
    scenario id, info to be included as prefix on generated output

    -------
    returned_data : the pu.dat input file

    """

    # open hex.shp file with set crs
    shp_crs_path = glob(os.path.join(path, eco + '_espg_*.shp'))

    # create df based on hexfile.shp
    shp_crs_layer = gpd.read_file(shp_crs_path[0])

    # create new column in .shp for 'id'
    shp_crs_layer.insert(0, 'id', range(1, 1 + len(shp_crs_layer)))

    # set values in column 'Cost' to 1, and column 'Status' to = 0
    shp_crs_layer["Cost"] = 1
    shp_crs_layer["Status"] = 0

    # create pu.dat file
    pu_dat = shp_crs_layer[["id", "Cost", "Status"]].set_index("id")
    output = pu_dat.to_csv('pu.dat')
    print(scen_id + ": pu.dat file successfully created")
    return output


# In[6]:


# # funtion to create spec.dat file

# def create_spec_dat(info_df, eco, prop=0.3, spf=1, minclump=False):
#     """
#     To create the spec.dat file, which stores information about ecosytem to be
#     analyzed in marxan run

#     Parameters
#     ----------
#     info_df : df
#     dataframe of ecosystem info, including 'Short_Name', 'US_km2' and
#     'Current_IUCN_TH' columns

#     eco : str
#     the abbreviated one word short name used for ecosystem being analyzed;
#     identifies a subdirectory of the timestamped marxan run directory

#     prop : float
#     The proportion of the total amount of the feature which must be included
#     in the solution; must be between 0 and 1 (tutorial suggests 0.3)

#     spf : int
#     species penalty factor

#     minclump : bool
#     determines if additional field 'target2' should be included in results,
#     to show the 'Minimum clump size for the representation of conservation
#     features in the reserve system'. Default set to 'False'.  If 'True',
#     'target2' column is added, with the minimum clump size calculated as the
#     ecosystem's extent in meters ('US_km2' * 1,000,000) by it's current IUCN
#     Threshold value ('Current_IUCN_TH' = .05 if CR or EN, or 0.10 if VU)
#     -------

#     returned_data : the spec.dat input file

#     """
#     # set columns of spec.dat, if minclump parameter is False
#     if minclump == False:
#         data = [{'id': 1, 'prop': prop, 'spf': 1, 'name': eco}]
#     # include add'l 'target2' column in file if minclump parameter is True
#     else:
#         target2 = info_df.at[eco,'Current_IUCN_TH'] * (
#             info_df.at[eco,'US_km2'] * 1000000)
#         data = [{'id': 1, 'prop': prop, 'target2': target2, 'spf': 1,
#                  'name': eco}]
#     # set index, and save file as 'spec.dat'
#     spec_dat = pd.DataFrame(data).set_index('id')
#     output = spec_dat.to_csv('spec.dat')
#     print(eco + ": spec.dat file successfully created")
#     return output


# In[7]:


# 2nd try - funtion to create spec.dat file

def create_spec_dat_v2(info_df, prop, target2, spf, eco, scen_id):
    """
    To create the spec.dat file, which stores information about ecosytem to be
    analyzed in marxan run

    Parameters
    ----------
    info_df : df
    dataframe of ecosystem info, including 'Short_Name', 'US_km2' and
    'Current_IUCN_TH' columns

    prop : float
    The proportion of total ecosystem area that must be included in solution

    target2: float
    minimum clumpsize of area, in order to be included in solution (*KBA*)

    spf : int
    species penalty factor

    eco : str
    the abbreviated one word short name used for ecosystem being analyzed;
    identifies a subdirectory of the timestamped marxan run directory

    scen_id : str
    scenario id, info to be included as prefix on generated output
    -------

    returned_data : the spec.dat input file

    """
#     # set columns of spec.dat, if minclump parameter is False
#     if minclump == False:
#         data = [{'id': 1, 'prop': prop, 'spf': spf, 'name': eco}]
#     # include add'l 'target2' column in file if minclump parameter is True
#     else:

    # this is where the KBA threshold can be tested -
    # KBA @ 100% = info_df.at[eco,'Current_IUCN_TH'] * (info_df.at[eco,'US_km2'] * 1000000
    # iterate by - test_threshold = [1.0, 0.75, 0.50, 0.25]
    # ITERATION WILL HAPPEN OUTSIDE FORMULA - SOMEWHERE IN MAIN LOOP
    # target2 = (target2 * 1000000) / 0.30
    data = [{'id': 1,
             'prop': prop,
             'spf': spf,
             'target2': target2,
             'name': eco
            }]
    # set index, and save file as 'spec.dat'
    spec_dat = pd.DataFrame(data).set_index('id')
    output = spec_dat.to_csv('spec.dat')
    print(scen_id + ": spec.dat file successfully created")
    return output


# In[8]:


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
    the abbreviated one word short name used for ecosystem being analyzed;
    identifies a subdirectory of the timestamped marxan run directory

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
           output = target_df.to_csv(outpath)
    print(eco + ": targets files created for each test threshold value")
    return output


# In[9]:


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
def create_input_dat(dest, blm, scen_id):
    """
    To create the input.dat file that stores processing parameters

    Parameters
    ----------
    dest : str
    directory input.dat file will be saved to

    prop : float
    must be a number between 0 and 1; represents the proportion of PU to be
    included in the initial reserve (default value is 0.5)

    scen_id : str
    scenario id, info to be included as prefix on generated output files

    other parameters will be added to replace the default initial values
    that are included in the QMarxan code

    -------
    returned_data : the input.dat file

    """
    output = os.path.join(dest,'input.dat')
    f = open(output, 'w')
    f.write("Input file for Annealing program.\n")
#     f.write('\n')
    f.write('This file generated for KBA Threshold\n')
    f.write('Analysis project using code from\n')
    f.write('QMarxan Toolbox 2.0\n')
    f.write('created by Apropos Information Systems Inc.\n')
#     f.write('\n')
    f.write("General Parameters\n")
    f.write("BLM " + str(blm) + "\n") # Boundary Length Modifier
    f.write("PROP %s\n" % formatAsME(0.5)) # Proportion of PU selected 1st run
    f.write("RANDSEED -1\n") # Random seed number
    f.write("NUMREPS 100\n") # Num of repeat runs (or solutions)
#     f.write('\n')
    f.write("Annealing Parameters\n")
    f.write("NUMITNS 1000000\n") # Num of iterations for annealing
    f.write("STARTTEMP %s\n" % formatAsME(-1.0)) # start temp for annealing
    f.write("COOLFAC %s\n" % formatAsME(-1.0)) # cooling factor for annealing
    f.write("NUMTEMP 10000\n") # num of temp decreases for annealing
    f.write("")
    f.write("Cost Threshold\n")
    f.write("COSTTHRESH %s\n" % formatAsME(0.0)) # cost threshold
    f.write("THRESHPEN1 %s\n" % formatAsME(0.0)) # size of cost thresh penalty
    f.write("THRESHPEN2 %s\n" % formatAsME(0.0)) # shp of cost thresh penalty
#     f.write("\n")
    f.write("Input Files\n")
    f.write("INPUTDIR input\n") # name of dir containing input files
    f.write("SPECNAME spec.dat\n") # Conservation Feature File
    f.write("PUNAME pu.dat\n") # Planning Unit File
    f.write("PUVSPRNAME puvsp.dat\n") # PU vs Conservation Feature File
    f.write("BOUNDNAME bound.dat\n") # Boundary Length File
#     f.write("BLOCKDEFNAME blockdef.dat\n") # Block Definition File
#     f.write("MATRIXSPORDERNAME puvsp_sporder.dat\n") # PUVSPR ordered by SP
    f.write("SCENNAME " + scen_id + "\n") # Scenario name for saved output
    f.write("SAVERUN 3\n") # Save each run (1-.dat, 2-.txt, 3-.csv)
    f.write("SAVEBEST 3\n") # Save the best run (1-.dat, 2-.txt, 3-.csv)
    f.write("SAVESUMMARY 3\n") # Save summary info (1-.dat, 2-.txt, 3-.csv)
    f.write("SAVESCEN 3\n") # Save scenario info (1-.dat, 2-.txt, 3-.csv)
    f.write("SAVETARGMET 3\n") # Save targets met information
    f.write("SAVESUMSOLN 3\n") # Save summed solution info (1 dat,2 txt,3 csv)
    f.write("SAVEPENALTY 3\n") # Save computed feature penalties
    f.write("SAVELOG 3\n") # Save log files (1-.dat, 2-.txt, 3-.csv)
    f.write("SAVESNAPSTEPS 0\n") # Save snapshots of each n steps
    f.write("SAVESNAPCHANGES 0\n") # Save snapshots after every n change
    f.write("SAVESNAPFREQUENCY 0\n") # Frequency of snapshots if used
    f.write("SAVESOLUTIONS MATRIX 3\n") # Save all runs in a single matrix
    f.write("OUTPUTDIR output\n") # name of dir containing output files
#     f.write("\n")
    f.write("Program control\n")
    f.write("RUNMODE 1\n") # Run option
    f.write("MISSLEVEL %s\n" % formatAsME(1.0)) # Species missing proportion
    f.write("ITIMPTYPE 1\n") # Iterative improvement
    f.write("HEURTYPE -1\n") # Heuristic
    f.write("CLUMPTYPE 0\n") # Clumping rule
    f.write("VERBOSITY 3\n") # Screen output
#     f.write("\n")
    f.close()
    print(os.path.basename(dest) + ": input.dat created successfully")
    return output


# In[10]:


# TO CREATE summary of the marxan run, incl info from the output file scenario
# details (sen.dat), input variables within the workflow (ALSO ADD AREA &
# OUTPUT STATS)
def create_mxrun_summary(dest, espg, prop, blm, target2, spf, scen_id, eco, df):
    """
    To create an output summary, showing local variables and info from sen.dat
    output file

    Parameters
    ----------
    dest : str
    path to the 'eco' subdirectory

    espg : str
    espg number (we're using ESPG:5070)

    prop : float - USE BLM HERE INSTEAD
    The proportion of the total amount of the feature which must be included
    in the solution; must be between 0 and 1 (tutorial suggests 0.3)

    target2 : float
    the min acceptable clump size - SHOULD EQUAL KBA (5 or 10 % x test thresh)

    spf : int
    species penalty factor

    scen_id : str
    scenario id, info to be included as prefix on generated output files

    df : df
    provided df with info about ecosystem's RLE status and area extent

    other parameters may be added to replace the default initial values
    that are included in the QMarxan code

    -------
    returned_data : the input.dat file

    """
    # display info from sen.dat output file
    sen_path = glob(os.path.normpath(os.path.join(dest, "output", "*_sen.*")))
    sen_df = pd.read_table(sen_path[0], header=None)
    sen_l1 = sen_df[0].iloc[0]
    sen_l2 = sen_df[0].iloc[1]
    sen_l3 = sen_df[0].iloc[2]
    sen_l4 = sen_df[0].iloc[3]
    sen_l5 = sen_df[0].iloc[4]
    sen_l6 = sen_df[0].iloc[5]
    sen_l7 = sen_df[0].iloc[6]
    sen_l8 = sen_df[0].iloc[7]
    sen_l9 = sen_df[0].iloc[8]
    sen_l10 = sen_df[0].iloc[9]
    sen_l11 = sen_df[0].iloc[10]
    sen_l12 = sen_df[0].iloc[11]
    sen_l13 = sen_df[0].iloc[12]
    sen_l14 = sen_df[0].iloc[13]
    sen_l14 = sen_df[0].iloc[14]
    sen_l15 = sen_df[0].iloc[15]

    output = os.path.join(dest, 'output', scen_id + '_mxrun_summary.dat')
    f = open(output, 'w')
    f.write("Scenario Details\n")
    f.write(sen_l1 + "\n")
    f.write(sen_l2 + "\n")
    f.write(sen_l3 + "\n")
    f.write(sen_l4 + "\n")
    f.write(sen_l5 + "\n")
    f.write(sen_l6 + "\n")
    f.write(sen_l7 + "\n")
    f.write(sen_l8 + "\n")
    f.write(sen_l9 + "\n")
    f.write(sen_l10 + "\n")
    f.write(sen_l11 + "\n")
    f.write(sen_l12 + "\n")
    f.write(sen_l13 + "\n")
    f.write(sen_l14 + "\n")
    f.write(sen_l15 + "\n")
    f.write("\n")
     # display info from stored variables in workflow
    blmstr = str(blm)
    propstr = str(prop)
    tgt2str = str(target2)
    spfstr = str(spf)
    f.write('Variables Set Locally -\n')
    f.write('ESPG value for raster and shapefile: ' + espg + "\n")
    f.write('input file variables:\n')
    f.write('BLM: ' +  blmstr + "\n")
    f.write('prop: ' +  propstr + "\n")
    f.write('target2: ' +  tgt2str + " (" + str(target2/1000000) + " km2)\n")
    f.write('Species Penalty Factor: ' + spfstr +'\n')
    f.write('scen_id: ' + scen_id + '\n')
    f.write("\n")

    # display data from eco df
    f.write('Spatial Extent of Ecosystem & KBA Thresholds\n')
    f.write('eco: ' + eco + '\n')
    f.write('US_km2: ' + str(df.at[eco,'US_km2']) + '\n')
    f.write('RLE_FINAL: ' + df.at[eco,'RLE_FINAL'] + '\n')
    f.write('Current_IUCN_TH: ' + str(df.at[eco,'Current_IUCN_TH']) + '\n')
    f.write('KBA @ 1.00 IUCN TH: ' + str(df.at[eco,'US_km2']*df.at[eco,'Current_IUCN_TH']) + " km2\n")
    f.write('KBA @ 0.75 IUCN TH: ' + str(0.75*(df.at[eco,'US_km2']*df.at[eco,'Current_IUCN_TH'])) + " km2\n")
    f.write('KBA @ 0.50 IUCN TH: ' + str(0.50*(df.at[eco,'US_km2']*df.at[eco,'Current_IUCN_TH'])) + " km2\n")
    f.write('KBA @ 0.25 IUCN TH: ' + str(0.25*(df.at[eco,'US_km2']*df.at[eco,'Current_IUCN_TH'])) + " km2\n")
    f.close()
    print(os.path.basename(dest) + (': mxrunsummary created successfully\n'
                                    'End run: ') + os.path.basename(dest) +
                                    '\n')
    return output


# In[11]:


# create function to create plot of best solution from Marxan output, and
# also save the shapefile merged with best solution as new file

def get_bestshp_and_bestplot(eco, path, espg, scen_id):
    """
    plots an image to show what hex cells were selected in the best run
    (also saves the shapefile merged with '_best.csv' needed to produce plot)

    Parameters
    ----------
    eco : str
    the abbreviated one word short name used for ecosystem being analyzed;
    identifies a subdirectory of the timestamped marxan run directory

    path : filepath
    filepath to ecosystem subdirectory

    espg : str
    espg number (we're using ESPG:5070)

    scen_id : str
    scenario id, info to be included as prefix on generated output files

    Returns
    -------
    output : list
    will generate two files named in output list
    (shapefile merged with best solution, and plot of best solution)
    """
    # Open raster data, set plotting extent
    raster_path = os.path.normpath(os.path.join(path,
                               "source_data",
                               eco + "_espg_" + espg + ".tif"))
    raster_layer = rxr.open_rasterio(raster_path, masked=True).squeeze()
    raster_extent = plotting_extent(raster_layer,
                                    raster_layer.rio.transform())

    # open shapefile created in the 'set_source_files_crs' function
    shp_path = os.path.normpath(os.path.join(
        path, "source_data", eco + "_espg_" + espg + ".shp"))
    shp_layer = gpd.read_file(shp_path)

    # open '_best' file created by Marxan and saved to 'output' directory
    globfile = glob(os.path.normpath(os.path.join(path, 'output', '*_best*')))

    if globfile == []:
        output = print (
            "ERROR: best run file not found - check output/log. \nWill need "
            "to resolve error and rerun Marxan if not completed successfully")

    else:
        best_run_path = globfile[0]
        best_run = pd.read_csv(best_run_path)

        # merge best_run df to shp layer
        shp_layer.insert(0, 'PUID', range(1, 1 + len(shp_layer)))
        shp_layer = shp_layer.merge(best_run, on='PUID')

        # open 'puvsp.dat' and merge with shp layer to get 'amount' from puvsp
        puvsp_path = os.path.normpath(os.path.join(path, 'input', 'puvsp.dat'))
        puvsp = pd.read_csv(puvsp_path)
        puvsp = puvsp.rename(columns={'pu': 'PUID'})
        shp_layer = shp_layer.merge(puvsp, on='PUID')

        fig_title_metr = shp_layer.query("SOLUTION == 1")['amount'].sum()/1000000
        ftm_string = str(fig_title_metr)

        # save merged shp as new file
        shp_w_best_and_amt = shp_layer.to_file(eco + "_w_best.shp", index=False)
        print (eco + '.shp merged with ' + scen_id +
               "_best.csv and puvsp.dat, saved as " + eco +
               "_w_best_and_amt.shp file")
        print ('preparing plots...')

        # create visualization showing hexcell selection from best run solution
        fig, ax = plt.subplots(figsize=(10, 10))
        shp_layer.plot(column='SOLUTION', cmap='tab20', ax=ax, alpha=0.65)
        ax.imshow(raster_layer, cmap='jet', extent=raster_extent,
                  interpolation='nearest')
        ax.set(title= scen_id + ': best run solution' +
               '\nTotal Selected Ecoystem = ' + ftm_string + ' sq km')
        ax.set_axis_off()

        best_plot = plt.savefig('best_plot.png', facecolor='w', edgecolor='k',
                                dpi=1200)
        plt.close(fig)
        print (scen_id + ": best plot saved as .png")

        output = (shp_w_best_and_amt, best_plot)
    return output


# In[12]:


# create function to create plot of summed solution from Marxan output, and
# also save the shapefile merged with best solution as new file

def get_ssolnshp_and_ssolnplot(eco, path, espg, scen_id):
    """
    plots an image to show hex cell selection frequency
    (also saves the shapefile merged with '_ssoln.csv' needed to produce plot)

    Parameters
    ----------
    eco : str
    the abbreviated one word short name used for ecosystem being analyzed;
    identifies a subdirectory of the timestamped marxan run directory

    path : filepath
    filepath to ecosystem subdirectory

    espg : str
    espg number (we're using ESPG:5070)

    scen_id : str
    scenario id, info to be included as prefix on generated output files

    Returns
    -------
    output : list
    will generate two files named in output list
    (shapefile merged with summed solution, and plot of summed solution)
    """
    # Open raster data, set plotting extent
    raster_path = os.path.normpath(os.path.join(path,
                               "source_data",
                               eco + "_espg_" + espg + ".tif"))
    raster_layer = rxr.open_rasterio(raster_path, masked=True).squeeze()
    raster_extent = plotting_extent(raster_layer,
                                    raster_layer.rio.transform())

    # open shapefile created in the 'set_source_files_crs' function
    shp_path = os.path.normpath(os.path.join(
        path, "source_data", eco + "_espg_" + espg + ".shp"))
    shp_layer = gpd.read_file(shp_path)

    # open '_ssoln' file created by Marxan and saved to 'output' directory
    globfile = glob(os.path.normpath(os.path.join(path, 'output',
                                                  '*_ssoln*')))
    if globfile == []:
        output = print (
            "ERROR: summed solutions file not found - check output/log. \n"
            "Will need to resolve error and rerun Marxan if not completed "
            "successfully")
    else:
        ssoln_path = globfile[0]
        ssoln = pd.read_csv(ssoln_path)
        ssoln = ssoln.rename(columns={'planning_unit': 'PUID'})

        # merge ssoln df to shp layer
        shp_layer.insert(0, 'PUID', range(1, 1 + len(shp_layer)))
        shp_layer = shp_layer.merge(ssoln, on='PUID')
        shp_w_ssoln = shp_layer.to_file(eco + "_w_ssoln.shp", index=False)
        print (eco + '.shp merged with ' + scen_id + "_ssoln.csv, saved as " + eco
               + "_w_ssoln.shp file")
        print ('preparing plots...')

        # create visualization showing hexcell selection from summed solution
        fig, ax = plt.subplots(figsize=(10, 10))
        shp_layer.plot(column='number', cmap='viridis', ax=ax, alpha=0.65)
        ax.imshow(raster_layer, cmap='jet', extent=raster_extent,
                  interpolation='nearest')
        ax.set(title= scen_id + ': summed solution' +
               '\n(hex cell selection frequency)')
        ax.set_axis_off()

        ssoln_plot = plt.savefig('ssoln_plot.png', facecolor='w', edgecolor='k',
                                dpi=1200)
        plt.close(fig)
        print (scen_id + ": ssoln plot saved as .png")

        output = (shp_w_ssoln, ssoln_plot)
    return output


# FUNCTION TO CREATE PUVSP.DAT & SPEC.DAT
# (PUVSP_SPORDER.DAT not needed, as currently we are looking at only one ecosystem at a time)
#
# * find code to replicate the Zonal Histogram tool in QGIS (Processing Toolbox > Raster Analysis > Zonal Histogram)
# https://docs.qgis.org/3.22/en/docs/user_manual/processing_algs/qgis/rasteranalysis.html#zonal-histogram
# import processing
# processing.run("algorithm_id", {parameter_dictionary})
#
# I think this comes from pyQGIS, but I don't think that is included with the earth-analytics-python-env??
#
# Here's info on running pyQGIS in Jupyter
# https://lerryws.xyz/posts/PyQGIS-in-Jupyter-Notebook
#
# here's my initial guess, using parameters copied from QGIS 'Zonal Histograms' log -
# from qgis.core import processing
# processing.run("native:zonalhistogram", { 'COLUMN_PREFIX' : '', 'INPUT_RASTER' : 'F:/NatureServe/LanasData/raster/foothill_r.tif', 'INPUT_VECTOR' : 'F:/NatureServe/LanasData/marxan_prep/foothill/pulayercws.shp', 'OUTPUT' : 'F:/NatureServe/pulayerfeatures.shp', 'RASTER_BAND' : 1 })
#
# * add column, multiplying pixel count by raster pixel area variable (our data's is 900, 30m x 30m = 900 sq m/pixel), this gives total extent of ecosystem within each individual planning unit hex of the hexfile.shp (ex. if pixelcount = 5, area = 4500, or 5 x 900)
#
# * use this as the source info for qmarxan 'export feature files' function
# input parameters copied from QGIS -
# Input parameters:
# { 'FEAT_FIELDS' : ['7147'], 'OUT_DIR' : 'F:\\NatureServe\\524 QM test\\input', 'PU_FIELD' : 'PUID', 'PU_LAYER' : 'F:/NatureServe/pulayerfeatures.shp' }
#
