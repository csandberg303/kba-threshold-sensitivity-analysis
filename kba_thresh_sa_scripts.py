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
import numpy as np
import pandas as pd


# In[2]:


# create function to get Lana's input files (created with ArcGIS) from the
# repo to local directory
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
        print(eco + ": " + file + " successfully copied from url")
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
def get_source_files(path, eco):
    """
    path : str
    local directory where the shapefiles or rasters are stored

    eco : str
    the abbreviated one word short name used for ecosystem
    """
    source_file_ls = glob(os.path.join(path, eco + '*'))
    if source_file_ls == []:
        print("no files found in " + path + "with expected name " + eco + "?")
    else:
        for file in source_file_ls:
            shutil.copy(file, os.getcwd())
            print(eco + ": "+ os.path.basename(file) + " copied successfully")


# In[4]:


# FUNCTION TO CREATE PU.DAT FILE
# create df based on hexfile.shp
# 3 columns - puid, cost & status
# No. of rows to equal number of hex cells in the hexfile.shp
# save as a csv .dat

def create_pu_dat(eco, path):
    """
     To create the pu.dat file that stores information about planning units in
     hex grid

     Parameters
     ----------
     eco : str
     name of ecosystem that will be analyzed by Marxan

     path : str
     local directory where 'hex_shp' directory is stored

     -------
     returned_data : the pu.dat input file

    """

    # open hex.shp file
    shp_data_path = os.path.join(path, "hex_shp", eco + '.shp')
    shp_layer = gpd.read_file(shp_data_path)

    # Reproject CRS to ESPG 5070
    shp_layer_5070 = shp_layer.to_crs(epsg='5070')

    # create columns for 'pu_id', cost' (value = 1) and 'status' (value = 0)
    shp_layer_5070.insert(0, 'pu_id', range(1, 1 + len(shp_layer_5070)))
    shp_layer_5070["cost"] = 1
    shp_layer_5070["status"] = 0

    # create pu.dat file
    pu_dat = shp_layer_5070[["pu_id", "cost", "status"]].set_index("pu_id")
    output = pu_dat.to_csv('pu2.dat')
    print(eco + ": pu.dat file successfully created")
    return output


# In[5]:


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
     print(eco + ": targets files created for each test threshold value")
     return csv


# In[6]:


# WRITE FUNCTION TO CREATE PUVSP.DAT & SPEC.DAT
# (PUVSP_SPORDER.DAT not needed, as currently we are working with just one one
# ecosystem at a time, rather than looking at mulitiple conservation features
# in a single run)

# find code to replicate the Zonal Histogram tool in QGIS
# (Processing Toolbox > Raster Analysis > Zonal Histogram)
# https://docs.qgis.org/3.22/en/docs/user_manual/processing_algs/qgis/
# rasteranalysis.html#zonal-histogram
# import processing processing.run("algorithm_id", {parameter_dictionary})

# I think the 'import processing' code suggestion comes from pyQGIS, but I
# don't think that is included with the earth-analytics-python-env??

# Here's info on running pyQGIS in Jupyter
# https://lerryws.xyz/posts/PyQGIS-in-Jupyter-Notebook

# here's my initial guess at some of the pyQGIS code, using parameters copied
# from the QGIS 'Zonal Histograms' log window -
# from qgis.core import processing processing.run(
#     "native:zonalhistogram", {
#         'COLUMN_PREFIX' : '',
#         'INPUT_RASTER' : 'F:/NatureServe/LanasData/raster/foothill_r.tif',
#         'INPUT_VECTOR' : 'F:/NatureServe/LanasData/marxan_prep/foothill/pulayercws.shp',
#         'OUTPUT' : 'F:/NatureServe/pulayerfeatures.shp',
#         'RASTER_BAND' : 1
#     })

# the result will show the number of raster pixels within each hexgrid cell.

# add new column, multiplying this pixel count by the raster pixel area
# variable, to give the total extent of ecosystem within each individual
# planning unit hex cell.
# the pixel area can be determined by looking at the raster saved to the
# 'source_data' directory in an earlier formula. Maybe a new function should
# be written to get this value, and that function would be called within this
# current 'create_puvsp_dat' function?
# (our data's pixel area is 900, as 30m x 30m = 900 sq m/pixel. If the
# pixelcount = 5, area = 4500, or 5 x 900)

# use this table (dataframe?) as the source info for qmarxan 'export feature
# files' function (WHICH WOULD ALSO CREATE THE SPEC.DAT FILE)

# input parameters copied from QGIS 'export_features_files' log window -
# Input parameters: {
#     'FEAT_FIELDS' : ['7147'],
#     'OUT_DIR' : 'F:\NatureServe\524 QM test\input',
#     'PU_FIELD' : 'PUID',
#     'PU_LAYER' : 'F:/NatureServe/pulayerfeatures.shp'
# }


# In[7]:


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
def create_input_dat(dest, prop, scen_name):
    """
     To create the input.dat file that stores processing parameters

     Parameters
     ----------
     dest : str
     directory input.dat file will be saved to

     prop : float
     must be a number between 0 and 1; represents the proportion of PU to be
     included in the initial reserve (default value is 0.5)

     scen_name : str
     scenario name, to be included as ID on generated output files

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
    f.write("BLM 1\n") # Boundary Length Modifier
    f.write("PROP %s\n" % formatAsME(prop)) # Proportion of PU (or sub TARGET)
    f.write("RANDSEED -1\n") # Random seed number
    f.write("NUMREPS 100\n") # Num of repeat runs (or solutions)
    f.write('\n')
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
    f.write("\n")
    f.write("Input Files\n")
    f.write("INPUTDIR input\n") # name of dir containing input files
    f.write("SPECNAME spec.dat\n") # Conservation Feature File
    f.write("PUNAME pu.dat\n") # Planning Unit File
    f.write("PUVSPRNAME puvsp.dat\n") # PU vs Conservation Feature File
    f.write("BOUNDNAME bound.dat\n") # Boundary Length File
    f.write("BLOCKDEFNAME blockdef.dat\n") # Block Definition File
    f.write("MATRIXSPORDERNAME puvsp_sporder.dat\n") # PUVSPR ordered by SP
    f.write("SCENNAME " + scen_name + "\n") # Scenario name for saved output
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
    f.write("\n")
    f.write("Program control\n")
    f.write("RUNMODE 1\n") # Run option
    f.write("MISSLEVEL %s\n" % formatAsME(1.0)) # Species missing proportion
    f.write("ITIMPTYPE 1\n") # Iterative improvement
    f.write("HEURTYPE -1\n") # Heuristic
    f.write("CLUMPTYPE 0\n") # Clumping rule
    f.write("VERBOSITY 3\n") # Screen output
    f.write("\n")
    f.close()
    return output


# # Write formula to create 'bound.dat' file
# (# code below taken from 'qmarxan_toolbox_algorithm.py')
#
# # # Constants used to refer to parameters and outputs. They will be
# # # used when calling the algorithm from another algorithm, or when
# # # calling from the QGIS console.
#
# # PU_LAYER = 'PU_LAYER'
# # PU_FIELD = 'PU_FIELD'
# # BND_METHOD = 'BND_METHOD'
# # BND_TREAT = 'BND_TREAT'
# # BND_VALUE = 'BND_VALUE'
# # CALC_FIELD = 'CALC_FIELD'
# # CALC_METHOD = 'CALC_METHOD'
# # TOL = 'TOL'
# # OUT_DIR = 'OUT_DIR'
#
# # def create_bound_dat(???self, config???):
# #         """
# #         Here we define the inputs and output of the algorithm, along
# #         with some other properties.
# #         """
# #         # pu layer
# #         self.addParameter(
# #                 self.PU_LAYER,
# #                 self.tr('Planning unit layer (source for bound.dat file)'),
# #                 [QgsProcessing.TypeVectorPolygon]
# #             )
# #         )
# #         # pu id
# #         self.addParameter(
# #             QgsProcessingParameterField(
# #                 self.PU_FIELD,
# #                 self.tr('Planning unit id field'),
# #                 parentLayerParameterName=self.PU_LAYER,
# #                 type=QgsProcessingParameterField.Numeric
# #             )
# #         )
# #         #
# #         # advanced settings
# #         #
# #         #  bnd method
# #         bndMethod = QgsProcessingParameterEnum(
# #             self.BND_METHOD,
# #             self.tr('Boundary method (how lengths between planning units will be set)'),
# #             options = ["Single","Measured","Weighted","Field"],
# #             defaultValue = 0
# #         )
# #         bndMethod.setFlags(bndMethod.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
# #         self.addParameter(bndMethod)
# #         # bnd treatment
# #         bndTreatment = QgsProcessingParameterEnum(
# #             self.BND_TREAT,
# #             self.tr('Boundary treatment (how values for PUs on perimeter of study area will be set)'),
# #             options = ["Full Value","Half Value","Exclude"],
# #             defaultValue = 0
# #         )
# #         bndTreatment.setFlags(bndTreatment.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
# #         self.addParameter(bndTreatment)
# #         # single value
# #         bndValue = QgsProcessingParameterNumber(
# #             self.BND_VALUE,
# #             self.tr('Boundary value (value for all boundaries regardless of measured length)'),
# #             type=QgsProcessingParameterNumber.Integer,
# #             minValue=0,
# #             defaultValue=1,
# #             optional=True
# #         )
# #         bndValue.setFlags(bndValue.flags() | QgsProcessingParameterDefinition.FlagAdvanced )
# #         self.addParameter(bndValue)
# #         # calculation field
# #         calcField = QgsProcessingParameterField(
# #             self.CALC_FIELD,
# #             self.tr('Calculation field (field to weight or assign boundary lengths)'),
# #             parentLayerParameterName=self.PU_LAYER,
# #             type=QgsProcessingParameterField.Numeric,
# #             optional = True
# #         )
# #         calcField.setFlags(calcField.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
# #         self.addParameter(calcField)
# #         # calculation method
# #         calcMethod = QgsProcessingParameterEnum(
# #             self.CALC_METHOD,
# #             self.tr('Calculation method (how to assign boundary length if values between adjacent planning units differ)'),
# #             options = ["Mean","Maximum","Minimum"],
# #             defaultValue = 0,
# #             optional = True
# #         )
# #         calcMethod.setFlags(calcMethod.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
# #         self.addParameter(calcMethod)
# #         # rounding precision
# #         tolerance = QgsProcessingParameterEnum(
# #             self.TOL,
# #             self.tr('Export precision tolerance (in map units)'),
# #             options = ["100","10","1","0.1","0.01","0.001","0.0001","0.00001"],
# #             defaultValue = 3
# #         )
# #         tolerance.setFlags(tolerance.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
# #         self.addParameter(tolerance)
#
# #         # select output folder
# #         defDir = os.path.join(os.path.expanduser('~'),'marxanproj1','input')
# #         self.addParameter(
# #             QgsProcessingParameterFolderDestination(
# #                 self.OUT_DIR,
# #                 self.tr('Marxan input folder (place to write bound.dat file)'),
# #                 defDir,
# #                 optional=False
# #             )
# #         )
# #

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
