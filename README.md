[![arXiv](https://img.shields.io/badge/arXiv-2307.02474-b31b1b.svg)](https://arxiv.org/abs/2307.02474)
<a href="https://ascl.net/2402.004"><img src="https://img.shields.io/badge/ascl-2402.004-blue.svg?colorB=262255" alt="ascl:2402.004" /></a>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11112984.svg)](https://doi.org/10.5281/zenodo.11112984)


# Cosmologically-coupled-black-holes formation mass: Numerics (CCBH-Numerics)

Davi C. Rodrigues <br>
Federal University of Espirito Santo (UFES), Brazil <br>
Heidelberg University, Germany

## Purpose
For given observational data of binary black holes (BBHs) from gravitational waves, `CCBH-Numerics` computes the probability of existence of a single black hole (BH) with formation mass below a threshold (e.g., 2 $M_\odot$) considering that the BHs are subject to a cosmological couplying such that their masses increases as ([Farrah et al ApJL 2023](https://doi.org/10.3847/2041-8213/acb704)) 
```math
m_f = a^k \, m_i \, , 
```
where $m_i$ is the initial (formation) mass, $m_f$ is the final (observed) mass, $a$ is the cosmological scale factor and $k$ is a constant. It is also assumed that the detected BBHs are formed from stellar evolution, they should not be primordial BHs. The code works for NSBH pairs as well, but all the explanation is focused on BBHs.

`CCBH-Numerics` was built to use the unbiased population of BBHs, as given by the power-law-plus-peak (PLPP) profile, as the observational input. It also works with individual data from BBHs. The previous version of this code is called `CCBH-PLPP` and it can be found in this repository as a branch.

This code is one of the two codes used for the paper `Constraints on cosmologically coupled black holes from gravitational wave observations and minimal formation mass` by Amendola, Rodrigues, Kumar & Quartin [MNRAS (2024)](https://academic.oup.com/mnras/advance-article-abstract/doi/10.1093/mnras/stae143/7529208?utm_source=advanceaccess&utm_campaign=mnras&utm_medium=email), [arXiv 2307.02474](https://arxiv.org/abs/2307.02474). The other code is named [CCBH-direct](https://github.com/itpamendola/CCBH-direct) and it focus on the direct method, as discussed in the paper above cited.

**Quick start**: clone the repository in your machine and run one of the notebooks. Try starting with `CCBH-PLPP-FormationPDFforM1`.  

All the notebooks are **Mathematica** notebooks. Later I will provide **Jupyter** notebook equivalents.

## Files and folders descriptions

### Notebooks

All the notebooks (in `nb` format) are independent among themselves. Each of them is focused on a specific approach that leads to a specific result in the paper. 

### Folders

* `input`. Contains data that were not generated by `CCBH-Numerics`, that are necessary and that are here provided for convenience. In particular: 
  * `GWTC.csv` -  GWTC-3 [events data](https://www.gw-openscience.org/eventapi/html/GWTC/).
  * `GWlistPopulationExtended.csv` - classification from [[arXiv 2111.03634]](https://arxiv.org/abs/2111.03634)
  * `all_samples_PLPP_GWTC3.h5` - data on the PLPP parameters distribution

* `codes` contains specific codes in `wl` format that are part of `CCBH-Numerics` and that are used repeatedly by different notebooks.

* `auxiliary` contains files that were generated by `CCBH-Numerics` but that do not constitute the main output. They contain intermediary results helpful for running quickly the notebooks. All these files can be deleted and can be regenerated by the code. They are provided for convenience.

* `output` contains the main outputs.

<!--

 All the notebooks are provided in `nb` and `wl` formats. The former runs in the official Mathematica notebook, the latter can be read as plain text and executed in different environments, like Jupyter notebooks. To load a `wl` file in a Jupyter notebook use [Mathematica Engine](https://www.wolfram.com/engine/) and [Wolfram Language for Jupyter](https://github.com/WolframResearch/WolframLanguageForJupyter). To execute the `wl` files that are in the `notebooks_in_wl_format` folder, first move them to the root folder of the CCBH-PLPP code.  

-->

### Large files - Git LFS
This repository includes a few large files with extensions `h5` or `mx` that are provided for convinience but that are not essential. The largest is about 100 MBs. The largest files depend on [Git Large File Storage (LFS)](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage). If you use github desktop, everything should be handled automatically.

## Acknowledgements
I acknowledge support from Federal University of Espirito Santo (Brazil), Heidelberg University (Germany), CNPq (Brazil) and FAPES (Brazil).
