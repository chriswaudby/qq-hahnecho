# qq-hahnecho

[![DOI](https://zenodo.org/badge/415567153.svg)](https://zenodo.org/badge/latestdoi/415567153)

Code for [Analysis of Conformational Exchange Processes using Methyl-TROSY-Based Hahn Echo Measurements of Quadruple-Quantum Relaxation](https://mr.copernicus.org/preprints/mr-2021-60/).

This repository contains code for two analyses:
* `hahnecho-analysis` - analysis of field-dependent Hahn echo measurements for single residues (Fig. 5 in publication)
* `hahnecho-cpmg-analysis` - global fitting of Hahn echo and CPMG data (Fig. 6 in publication)

Analysis scripts are written in Julia. Dependencies for both analyses are contained in the `Project.toml` file, and can be set up from the package manager by running `] instantiate` in the top level directory. Individual analyses should then be run from within the appropriate directory.


## Analysis of individual field-dependent Hahn echo measurements

* Directory: `hahnecho-analysis`
* Main file: `go.jl`
* All data are contained within `data` subdirectory.


## Global Hahn echo and CPMG analysis

* Directory: `hahnecho-cpmg-analysis`
* Main file: `go.jl`
* All data are contained within `data` subdirectory. CPMG data are loaded from parsed ChemEx format.


## Citing

If you use this software, please cite it: see `CITATION.cff` for further information.
