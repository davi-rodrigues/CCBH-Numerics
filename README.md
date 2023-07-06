# Cosmological coupled black holes and the power-law-plus-peak distribution (CCBH-PLPP)

Davi C. Rodrigues (UFES, Brazil & Heidelberg University, Germany)

This code was used for the paper `Constraints on cosmologically coupled black holes from gravitational wave observations and minimal formation mass` by Amendola, Rodrigues, Kumar & Quartin (arXiv XXXXX). It was used in all PLPP-method applications, which is the method based on the power-law-plus-peak distribution.  See also the [Direct method](https://github.com/itpamendola/CCBH-direct) code.

**Quick start**: clone the repository in your local machine (or download the full zip file) and run one of the notebooks. Try starting with `CCBH-PLPP-FormationPDF`.


## Files and folders


### Notebooks
There are four independent notebooks:
* `CCBH-PLPP-FormationPDF` - Finds the "modified" PLPP distribution for m1. It also includes probability evaluations with variable minimal mass. If in doubt, start here.
* `CCBH-PLPP-FormationPDFforM2` - Finds the "modified" PLPP distribution for m2.
* `CCBH-PLPP-ExclusionPlots` - Generates the exclusion plots.
* `CCBH-PLPP-ParameterDependenceTest` - Considers random selection of PLPP parameters in according with the GWTC-3 parameter distribution. 

All the notebooks are provided in `nb` format (convenient for those with a Mathematica license) and a `wl` file that can be read as plain text and executed in different free environments, like Jupyter.

### `codes` folder
Contains five `wl` files that are called by the notebooks when needed:
* `ObsDataPreparationGWTC-3`
* `Constants`
* `PowerLawPlusPeakDefinition` 
* `Cosmology`
* `MassFactorCorrection`

### `input` folder
Contains data that are here analysed (all data from GWTC-3)

* `GWTC.csv` -  GW events data
* `GWlistPopulationExtended.csv` - classification
* `all_samples_PLPP_GWTC3.h5` - data on the PLPP parameters distribution.

### `output` folder
All the generated plots are saved in this folder. 

### `auxiliary` folder
Some intermediary steps that take time to run are saved in this folder. This folder is provided empty, apart from a readme file.

### `directories.wl`
This file has information on the directory structure and defines related functions. For convenience, it is not put in a subfolder.

 # Acknowledgements

This work was in part supported by CNPq (Brazil) and FAPES (Brazil).