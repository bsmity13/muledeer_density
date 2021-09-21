## Code for "An Eulerian Perspective on Green-Wave Surfing"

_Authors_:  

  * Tatum Del Bosco

  * Brian J. Smith <a itemprop="sameAs" content="https://orcid.org/0000-0002-0531-0492" href="https://orcid.org/0000-0002-0531-0492" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

  * Tal Avgar <a itemprop="sameAs" content="https://orcid.org/
0000-0002-8764-6976" href="https://orcid.org/
0000-0002-8764-6976" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

Submitted to *Ecology*.

## Data

### Camera data  

The raw camera trap data that go into this analysis are stored in the comma-delimited file `data/muledeer_density.csv`. <div style = 'color:red'>TBD: </div>An archived copy of this CSV is stored XX. Full metadata are available with the archived version, but we have included a brief description of each field here:  

  * `site` -- integer indicating the unique site number.
  * `x` -- x-coordinate of the site in UTMs (EPSG: 32612).
  * `y` -- y-coordinate of the site in UTMs (EPSG: 32612).
  * `date` -- date photos were recorded.
  * `doy` -- day-of-year as an integer (e.g., 1 = Jan 1; 365 = Dec 31).
  * `t` -- time index for the statistical model, starting at 1 on 15 March 2019 and incrementing by 1 each day.
  * `snow` -- boolean (0/1) indicating whether snow was present at that site on that date.
  * `camera` -- camera type (`cuddeback` or `reconyx`).
  * `start_date` -- first date camera was active at that site.
  * `ndvi` -- normalized difference vegetation index at that site on that date.
  * `irg` -- instantaneous rate of greenup at that site on that date.
  * `grass_forb` -- percent cover of grasses and forbs at that site.
  * `tree` -- percent cover of trees at that site.
  * `elev` -- elevation of the site in feet.
  * `live_pres` -- boolean (0/1) indicating whether livestock were present at that site on that date.
  * `end_date` -- last date camera was active at that site.
  * `photos` -- total number of deer recorded in photos (within the conduits, see Fig. 2 in main text) by that camera on that date. Values of `NA` indicate dates before `start_date` or after `end_date`.
  * `peak_irg` -- date when IRG peaked at that site.
  * `days_to_peak` -- difference in days between current date (`date`) and the date of peak IRG (`peak_irg`).
  * `elev_doy` -- composite variable equal to `elev` times `doy`
  * `elev_doy2` -- composite variable equal to `elev` times `doy^2`
  * `ndvi_grass_forb` -- composite variable equal to `ndvi` times `grass_forb`
  * `irg_grass_forb` -- composite variable equal to `irg` times `grass_forb`
  
### Geospatial data

Most of the geospatial data needed for this analysis have already been attached to the site data in `data/muledeer_density.csv`, as described in the main text of the manuscript. For the purposes of reproducing figures, some geospatial data are included here in the directory `geo`.  All geospatial data are projected using coordinate reference system EPSG:32612 (UTM Grid 12, WGS 84).

  * `sites.shp` -- shapefile containing locations of each site.
  * `elev.tif` -- raster of elevation in meters. Original data from The National Map [3DEP](https://www.usgs.gov/core-science-systems/ngp/3dep).
  * `grass_forb.tif` -- raster of percent cover of grass and forbs. Original data from [RAP](https://rangelands.app/products/#cover).
  * `tree.tif` -- raster of percent cover of trees. Original data from [RAP](https://rangelands.app/products/#cover).
  * `NDVI_stack.tif` -- multiband (3) raster with values of NDVI on 1 April 2019, 1 May 2019, and 1 June 2019. Original data from [MODIS product MOD09Q1](https://lpdaac.usgs.gov/products/mod09q1v006/).
  * `IRG_stack.tif` -- multiband (3) raster with values of IRG on 1 April 2019, 1 May 2019, and 1 June 2019. Original data from [MODIS product MOD09Q1](https://lpdaac.usgs.gov/products/mod09q1v006/).

## R Script Structure  

Scripts in this repository are meant to be run in numerical order. To save space in this repository, intermediate data files that are created by each script are saved as `data/*.rds`, and this pattern has been added to `.gitignore`. This project is comprised of the following scripts:

  * `01_prep_data.R` -- takes the raw data (`data/muledeer_density.csv`) and formats it for analysis using NIMBLE.
  * `02_fit_model.R` -- runs the MCMC estimation algorithm using [NIMBLE](https://r-nimble.org/). Sources the model definition in `99_NIMBLE_defn.R`.
  * `03_model_diag.R` -- runs convergence diagnostics on the chains of posterior samples, including the Gelman-Rubin statistic, traceplots, and density plots for the main model parameters.
  * `04_predict.R` -- predicts state variables across the input data using the fixed effects estimated in each iteration of the MCMC. Used later to produce Figs. 3 & 4 in the main text of the manuscript.
  * `05_map_pred.R` -- predicts state variables across all pixels of the study area using the fixed effects estimated in each iteration of the MCMC. Used later to produce Fig. 5 in the main text of the manuscript.
  * `06_figures.R` -- produces final figures that appear in the manuscript and supplementary information.
  * `99_fun.R` -- contains helper functions sourced into various other scripts. Not intended to run on its own.
  * `99_NIMBLE_defn.R` -- contains the NIMBLE model definition sourced in script 2. Not intended to run on its own.
  
