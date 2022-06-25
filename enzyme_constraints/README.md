# Enzyme-constrained genome-scale metabolic models for diverse yeast species
Collection of scripts for generation of ecModels for diverse yeast species and analysis of flux distributions, enzyme usage and parameterization properties. The yeast species encompassed in this study are:
- *Dekkera bruxellensis*
- *Eremothecium sinecaudum*
- *Kluyveromyces lactis*
- *Kluyveromyces marxianus*
- *Kluyveromyces dobzhanski*
- *Komagataella pastoris*
- *Lachancea fermentati*
- *Lachancea thermotolerans*
- *Naumovozyma castelli*
- *Saccharomyces eubyanus*
- *Schizosaccharomyces pombe*
- *Tetrapisispora blattae*
- *Tetrapisispora phaffii*
- *Zygosaccharomyces rouxii*

### Required Software:
* A functional Matlab installation (MATLAB 7.3 and higher)
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN) (version 2.3.0)
* For generating some figures, a functional [R installation](https://www.r-project.org/) (version 3.6.1)

### Dependencies - Recommended Software:
* The libSBML MATLAB API (version [5.17.0](https://sourceforge.net/projects/sbml/files/libsbml/5.17.0/stable/MATLAB%20interface/) is recommended)
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulation
* A [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

## Regenerating the ecModels:
The master script for generating the ecModels from the species-specific GEMs (provided in the `models/` subdirectory) is `code/update_ecModels.m`. Run this  script in MATLAB to regenerate the 14 ecModels. This pipeline generates:
- Unconstrained ecModels for incorporation of condition-dependent protein abundance data (under the name `ecModels/specific-species/ecModel.mat`)
- Constrained ecModels by a total protein pool available for metabolic enzymes (under the name `ecModels/specific-species/ecModel_batch.mat`). These models are able to predict the experimental batch growth rate for each species growing on glucose, however they are not curated for further phenotype predicitions.

## Curating the ecModels:
In order to generate biologically-meaningful phenotype predictions with the provided ecModels, these were curated according to physiology data (available at: `data/yeasts_parameters.txt`), by following an automated iterative curation of Kcat parameters based on sensitivity analysis over:
- growth rate
- ethanol production rate
- glucose uptake rate 

These curated model files are available in the `ecModels/` subdirectory, for each species, under the name `ecModel_improved.mat` and they are ready for FBA simulations and phenotype predictions.
* Now these models were compressed in zip files.
