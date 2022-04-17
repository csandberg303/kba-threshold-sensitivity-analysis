# ea-2022-kba-sensitivity
Sensitivity Analysis of KBA Ecosystem Criteria. A Study Project in Collaboration with NatureServe.

# Members 
[Christy Sandberg](https://github.com/csandberg303)
[Eric Nutt](https://github.com/ericrnutt5)
[Lana Kurakina](https://github.com/S-Kur)
[Elsa Culler - mentor](https://github.com/eculler)

# Purpose of the project
We will perform a sensitivity analysis of the current thresholds set by the IUCN, looking to see what KBAs can be found for each ecosystem type at 100%, 75%, 50% and 25% of the current thresholds.

Another factor we will look at is the natural spatial pattern of each ecosystem type (linear, patchy or matrix-forming), and how that pattern may affect how KBAs are identified.

Perhaps a lower threshold would more effectively trigger KBA identification in linear ecosystems, such as those found along rivers and coasts. If so, spatial pattern might be a characteristic to be included in future IUCN Guidelines.

# Installation instructions
We are using earth-analytics-enviornment provided and maintained by Earth Lab, University of Colorado. Please follow this [link](https://www.earthdatascience.org/workshops/setup-earth-analytics-python/) for installation instructions.

# Data
* Raster dataset featuring US ecosystems at 30 m spatial resolution, provided by NatureServe.
* Shapefile with a network of 7 sq. mile hexagons covering North America, also provided by NatureServe.
* Processed shapefiles for each indiviadual ecosystem, contain network of hexes showing spatial distribution of an ecosystem, unique hex IDs, and quantity of raster cells (of the initial dataset) which got a hexagon overlaps.

# Data Storage
* Currently data is stored in our shared Google Drive. User will have to download data and save it in earth-analytics/data/natureserve


# Workflow


![Workflow](https://github.com/S-Kur/ea-2022-kba-sensitivity/blob/main/assests/figures/KBA_Analysis_Workflow.jpg)
