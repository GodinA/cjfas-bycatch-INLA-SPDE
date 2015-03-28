# Background
This repo contains the necessary data and R codes to reproduce the analysis of the publication titled: 'Applying Bayesian spatiotemporal models to fisheries bycatch in the Canadian Arctic' published in the *Canadian Journal of Fisheries and Aquatic Sciences*, 2015, 72(2): 186-197, 10.1139/cjfas-2014-0159 available at: http://www.nrcresearchpress.com/doi/abs/10.1139/cjfas-2014-0159#.VRalgEZTTZg

# Contents of this directory
This directory contains two folders `Data` and `Rcodes`. 
* `Data`
  * `gridpred.txt` -  Grid (4 x 4 km) used in prediction
  *  `fackdata.txt` - Greenland shark bycatch data (counts) 
  *  `etopo1_bedrock.xyz` - NOAA ETOPO1 Bathymetry of Baffin Bay from https://www.ngdc.noaa.gov/mgg/global/global.html
  *  `EEZ` - Folder containing the shapefile for the Canadian Northwest Atlantic Exclusive Economic Zone (EEZ)
* `Rcodes`
  * `AllmodelsDICs.R` - R codes to run all spatiotemporal model tested (See Table 2 and 3 in paper)
  * `Suppl_Material.R` - R codes to fit the best model, predict bycatch hotspots, and create most figures included in the paper.

# Notes
As per data sharing agreement with Fisheries and Ocean Canada, `fackdata.txt` is not the original dataset used in the paper, but is similar - results will be slightly different.

# Other resources
For further information on R-INLA, please refer to their website: http://www.r-inla.org/
Information on the hurdle and mixture likelihood used in our modeling are available: http://www.math.ntnu.no/inla/r-inla.org/doc/likelihood/zeroinflated.pdf

