# Cosmological coupled black holes and the power-law-plus-peak distribution (CCBH-PLPP)

Davi C. Rodrigues (UFES, Brazil & Heidelberg University, Germany)

This code was used for the paper `Constraints on cosmologically coupled black holes from gravitational wave observations and minimal formation mass` by Amendola, Rodrigues, Kumar & Quartin [[arXiv 2307.02474]](https://arxiv.org/abs/2307.02474). It was used in all PLPP-method applications, which is the method based on the power-law-plus-peak distribution.  See also the [Direct method](https://github.com/itpamendola/CCBH-direct) code.

**Quick start**: clone the repository in your machine (or download the full zip file) and run one of the notebooks. Try starting with `CCBH-PLPP-FormationPDF`.  


## Files and folders


### Notebooks
There are four independent notebooks:
* `CCBH-PLPP-FormationPDF` - Finds the modified-PLPP distribution for m1. It also includes probability evaluations with variable minimal mass. If in doubt, start here.
* `CCBH-PLPP-FormationPDFforM2` - Finds the modified-PLPP distribution for m2.
* `CCBH-PLPP-ExclusionPlots` - Generates the exclusion plots.
* `CCBH-PLPP-ParameterDependenceTest` - Considers random selection of PLPP parameters in accordance with the GWTC-3 parameter distribution. 

All the notebooks are provided in `nb` and `wl` formats. The former runs in the official Mathematica notebook, the latter can be read as plain text and executed in different environments, like Jupyter notebooks. To execute load a `wl` file in a Jupyter notebook use [Mathematica Engine](https://www.wolfram.com/engine/) and [Wolfram Language for Jupyter](https://github.com/WolframResearch/WolframLanguageForJupyter). To execute the `wl` files, first move them to the root folder of the CCBH-PLPP code.

### `codes` folder
Contains five `wl` files that are called by the notebooks when needed:
* `directories` - define functions for file manipulation.
* `ObsDataPreparationGWTC-3` - load and prepare the GWTC-3 events for analysis.
* `Constants` - defines constants that are relevat for cosmology, the delay-time distribution and other computational parameters.
* `PowerLawPlusPeakDefinition` - defines the PLPP distribution.
* `Cosmology` - defines cosmology. For the first time run only, it will compute a table that should not take more than about 10-15 min. 
* `MassFactorCorrection` - Defines the main effect of CCBH and its relation to the delay-time distribution.

### `input` folder
All the data in this folder are provided for convenience, they are not part of the CCBH-PLPP code, they are simply used by the CCBH-PLPP code. All data from GWTC-3.

* `GWTC.csv` -  GWTC-3 [events data](https://www.gw-openscience.org/eventapi/html/GWTC/).
* `GWlistPopulationExtended.csv` - classification from [[arXiv 2111.03634]](https://arxiv.org/abs/2111.03634)
* `all_samples_PLPP_GWTC3.h5` - data on the PLPP parameters distribution.

### `output` folder
All the generated plots are saved in this folder. 

### `auxiliary` folder
Some intermediary steps that take time to run are saved in this folder. This folder is provided empty, apart from a readme file.

 # Acknowledgements

This work was in part supported by CNPq (Brazil) and FAPES (Brazil).