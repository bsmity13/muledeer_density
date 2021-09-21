## Code for "An Eulerian Perspective on Green-Wave Surfing"

_Authors_:  

  * Tatum Del Bosco

  * Brian J. Smith <a itemprop="sameAs" content="https://orcid.org/0000-0002-0531-0492" href="https://orcid.org/0000-0002-0531-0492" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

  * Tal Avgar <a itemprop="sameAs" content="https://orcid.org/
0000-0002-8764-6976" href="https://orcid.org/
0000-0002-8764-6976" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>

Submitted to *Ecology*.

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
  
