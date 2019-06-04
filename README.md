Data and code for Integrative Structure Determination using data from point mutant epistatic miniarray profile (pE-MAP). 

# pE-Map data and processing



# pE-Map restraint

# Modeling of the H3-H4 dimer

## Data for simulations (directory modeling_histones/data)

Data used for modeling includes:

1) Script and alignments to generate comparative mdoels.

2) Comparative models of H3 and H4

3) Processed pE-MAP data files containing the pairs of residues to which the pE-MAP distance restraints will be applied. The format of these files is:

```
protein1 protein2 residue1 residue2 MIC_value distance_in_xray_structure (if known)
```

4) Python scripts to shuffle and resample the pE-MAP data files


## Running the simulations (directory modeling_histones/modeling)

-  `top_his_comp_models.dat`: Topology file containing the representation of the h3-h4 system. Each protein has presented as a rigid body. Histones tails are not considering for modeling. 

-  `mod_pemap_histones.py`: PMI modeling scripts for running the production simulations. The search for good-scoring models relied on Replica Exchange Gibbs sampling, based on the Metropolis Monte Carlo (MC) algorithm.  We recommend producing at least 2,500,000 models from 50 independent runs, each starting from a different initial conformation of H3-H4 dimer to have proper statistics.

## Analysis of simulations (directory modeling_histones/analysis)
The analysis uses the PMI_analysis module in: https://github.com/salilab/PMI_analysis


1) `run_analysis_trajectories.py`: 

2) `run_extract_models.py`:

3) `run_clustering.sh`:

4) `run_clustering.sh`:

# Modeling of RNAP II

# Modeling of Bacterial RNAP II

## Running the simulations (directory XX)

## Analysis of simulations (directory XX)


## Information

_Author(s)_: Ignacia Echeverria, Hannes Braberg

_Date_: June 3rd, 2019

_License_: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.

_Last known good IMP version_: [![build info](https://integrativemodeling.org/systems/?sysstat=25&branch=master)](https://integrativemodeling.org/systems/) [![build info](https://integrativemodeling.org/systems/?sysstat=25&branch=develop)](https://integrativemodeling.org/systems/)

_Testable_: Yes

_Parallelizeable_: Yes

_Publications_:
Xiaorong Wang, Ilan E Chemmama, Clinton Yu, Alexander Huszagh, Yue Xu, Rosa Viner, Sarah Ashley Block, Peter Cimermancic, Scott D Rychnovsky, Yihong Ye, Andrej Sali, and Lan Huang
*Quantitative genetic interaction mapping informs integrative structure determination of biomolecular assemblies***XX**