# Neuron-Growth-Morphometrics-In-Vitro
[![DOI](https://zenodo.org/badge/473702722.svg)](https://zenodo.org/badge/latestdoi/473702722)

This is a semi-automated quantitative analysis method for assessing neuron morphology imaged using bright-field microscopy as described in "Quantitative evaluation of neuron developmental morphology in vitro using the change-point test" by AS Liao, W Cui, V Webster-Wood, and YJ Zhang (submitted to Neuroinformatics 2022). This method uses outputs from neurite traces obtained using NeuronJ, a plugin in ImageJ. It evaluates the traces using seven morphometrics:
* Common Morphometrics (used previously to discriminate between the morphology of neuron cell types) 
  * Total Length (Sum of All Neurite Lengths) Per Cell
  * Number of End Points Per Cell (Degree)
  * Number of Neurites Per Cell
  * Average Tortuosity Per Cell

* Change-Point-Test-Based Morphometrics
  * Number of Change Points Per Cell
  * Average Segment Length (Distance Between Change Points) Per Cell
  * Average Absolute Turning Angle Between Change Points Per Cell


## Table of Contents
* [General Info](#general-info)
* [Directories](#Directories)

## General Info

## Directories
### analysisCode
This directory includes the R code used to automatically quantitatively assess the morphometrics of the neurites. Use runNeuronQuant.R to do the complete analysis.

### exampleDataset
This directory includes an example dataset to run the analysis code (it is a smaller subset of what was analyzed in Liao et al. (2022)). There are two sub-directories inside: Unprocessed and processed. The Unprocessed directory includes the test data files in the appropriate file tree structure that runNeuronQuant.R expects. There is also a zipped version of the dataset.

### reportedDataset
This directory hosts the reported dataset from Liao et al. (2022). There are two sub-directories inside: allFeatures_fig7 and highlightedFeatures_fig5fig6. These two directories include the Jupyter files, data spreadsheet (.xlsx), and resulting .png and .svg files used to generate the violin plots found in Figures 5,6, and 7 in Liao et al. (2022). 
The full dataset can be found here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6415474.svg)](https://doi.org/10.5281/zenodo.6415474)

