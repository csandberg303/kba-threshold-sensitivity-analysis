# ea-2022-kba-sensitivity-analysis
Sensitivity Analysis of KBA Ecosystem Criteria. A Study Project in Collaboration with NatureServe.

# Members 
* [Christy Sandberg](https://github.com/csandberg303)
* [Eric Nutt](https://github.com/ericrnutt5)
* [Lana Kurakina](https://github.com/S-Kur)
* [Elsa Culler - mentor](https://github.com/eculler)

# License
We are using the standard Apache 2.0 License.

# Purpose of the project
Systematic Conservation Planning is a field of research that uses GIS and Python for spatial analysis, with the purpose of strategically identifying areas that meet a given conservation goal. This occurs by bringing collected data together into an iterative workflow to define areas that successfully meet the required criteria. Depending on the conservation goal, data can be analyzed from a diverse range of sources including species populations, vegetative land cover, nesting sites, seasonal feeding patterns or even the location of historic cultural heritage sites. Results of the analysis are shared with key stakeholders and decision makers, enabling them to efficiently direct limited resources towards solutions with the lowest costs and greatest chance for long term ecological success.

The International Union for Conservation of Nature (IUCN) has established a framework of criteria and thresholds for identifying Key Biodiversity Areas (KBAs) to conserve threatened species and threatened ecosystems. Our project will use Python to analyze vegetative land cover in the Continental United States (CONUS), and take a closer look at the spatial thresholds set by the IUCN for Threatened Ecosystems. We will perform a sensitivity analysis to see what KBAs can be found for a series of ecosystems using 100%, 75%, 50% and 25% of the current IUCN thresholds. 

Another factor we will look at is the natural spatial pattern of each ecosystem (linear, patchy or matrix-forming), and how that pattern may affect how KBAs are identified. Perhaps a lower threshold would more effectively trigger KBA identification in linear ecosystems, such as those found along rivers and coasts. If so, spatial pattern might be a characteristic to be included in future IUCN Guidelines for KBA identification.

In addition to providing information specific to our question of how the IUCN thresholds perform in defining KBAs for our series of CONUS ecoystems, We hope that providing well documented code using the marxancopy package will have the extended benefit of allowing others to explore our results and find new applications for the code using their own source data.  This code will efficiently create the needed directory structure and input files for Marxan analysis, which is more commonly accessed through the GUI interface of Marxan Connect or ArcGIS/QGIS.

# Installation instructions
### Installation of earth analytics python environment
We are using earth-analytics-enviornment provided and maintained by Earth Lab, University of Colorado. Please follow this [link](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/) for installation instructions.

### Installation of marxanconpy python package:
To install the marxanconpy python package, complete the steps below:
* activate your earth-analytics-conda environment
* install igraph and wxPython libraries
* AFTER that install marxanconpy -- NOTE -- There will be a small error message with a warning in the end, but it was installed successfully.

Here is the command line script necessary for the installation above:
* conda activate earth-analytics-python
* pip install igraph
* pip install wxPython
* pip install marxanconpy


# Required Tools and Packages
* ArcGIS and ArcMarxan Toolbox, Version 2.0.2; Available at https://aproposinfosystems.com/
* and/or QGIS and QMarxan Toolbox, Version 2.0.1
* matplotlib
* numpy
* geopandas
* xarray
* rioxarray

# Data
* Raster dataset featuring US ecosystems at 30 m spatial resolution, provided by NatureServe.
* Shapefile with a network of 7 sq. mile hexagons covering North America, also provided by NatureServe.
* Processed shapefiles for each indiviadual ecosystem, contain network of hexes showing spatial distribution of an ecosystem, unique hex IDs, and quantity of raster cells (of the initial dataset) which got a hexagon overlaps.

# Data Storage
* Sample data to run the workflow is provided in the assets > data > raster & shape hex folders. Rasters for three ecosystems and 7mi2 hex shapefile overlays of the three ecosystem rasters are made available there for the process of testing the workflow. These are all of the necessary files to run the marxanconpy workflow.


# Workflow


![Workflow](https://github.com/csandberg303/kba-threshold-sensitivity-analysis/blob/main/assets/figures/KBA_Analysis_Workflow.jpg)
