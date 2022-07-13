# Neuron-Growth-Morphometrics-In-Vitro
[![DOI](https://zenodo.org/badge/473702722.svg)](https://zenodo.org/badge/latestdoi/473702722)

This is a semi-automated quantitative analysis method for assessing neuron morphology imaged using bright-field microscopy as described in "Semi-automated quantitative evaluation of neuron developmental morphology in vitro using the change-point test" by AS Liao, W Cui, V Webster-Wood, and YJ Zhang (submitted to Neuroinformatics 2022). This method uses outputs from neurite traces obtained using NeuronJ, a plugin in ImageJ. It evaluates the traces using seven morphometrics:
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
* [Directories](#directories)

## General Info

The tool to complete the CPT and the subsequent morphological calculations can be found in the 'analysisCode' directory. This tool requires the trace coordinate outputs from NeuronJ (.txt file) and an associated .csv file for the neurite identification. In addition, metadata on the time point and the well ID that the trace data belongs to are present in the organization of directories that house the data, as seen in the 'unprocessed' directory of 'exampleDataset'. The neurite trace coordinate .txt files and the .csv files need to be the innermost directory (i.e. 'dir4_wellRowID'). To use the tool, run the runNeuronQuant.R script. The user will be prompted to enter a directory that houses the data, which should be the path for the 'dir0_experimentID' directory. Expected outputs can be found in the 'processed' directory. 

In addition, we have also included the R script (RStatTest.R) for our follow up statistical analysis to compare the morphometrics at different time points. This script was developed to run the Anderson-Darling test, the Kruskal-Wallis test, and the Dunn test with a Bonferroni correction for our organized dataset. This script has not been generalized and is specific for our datasheet (perCellMetrics(den10,rem8,remMiss).csv located in 'reportedDataset/allFeatures_fig8'), but was included for completeness in detailing the results of our manuscript. The results from this script were reported in Tables 2-16 in Liao et al. (2022).

## Directories
### 'analysisCode'
This directory includes the R code used to automatically quantitatively assess the morphometrics of the neurites. Use runNeuronQuant.R to do the complete analysis.

### 'exampleDataset'
This directory includes an example dataset to run the analysis code (it is a smaller subset of what was analyzed in Liao et al. (2022)). There are two sub-directories inside: 'unprocessed' and 'processed'. The 'unprocessed' directory includes the test data files in the appropriate file tree structure that runNeuronQuant.R expects. The 'processed' directory houses a copy of the data and directory tree structure as 'unprocessed' and includes the expected outputs after running 'runNeuronQuant.R' and inputting the path to 'dir0_experimentID'.

The directory tree structure can be found below. The outermost directory is 'dir0_experimentID', which should house all the data from a given experiment. The next directory ('dir1_div###') organizes the data into the time points in which they were imaged. After that, the 'dir2_plateID#' is used to group the data by the well plate ID, which would be useful if multiple plates were used simultaneously. Following that, the 'dir3_wellColumnID#' and 'dir4_wellRowID#' further group the data based on the well ID within a given plate. This enables multiple groups of data to be analyzed with a single run of 'runNeuronQuant.R'.

    └───dir0_experimentID
        ├───dir1_div000
        │   ├───dir2_plateID0
        │   │   ├───dir3_wellColumnID0
        │   │   │   ├───dir4_wellRowID0
        │   │   │   └───dir4__wellRowID1
        │   │   └───dir3_wellColumnID1
        │   │       ├───dir4_wellRowID0
        │   │       └───dir4_wellRowID1
        │   └───dir2_plateID1
        │       ├───dir3_wellColumnID0
        │       │   ├───dir4_wellRowID0
        │       │   └───dir4_wellRowID1
        │       └───dir3_wellColumnID1
        │           ├───dir4_wellRowID0
        │           └───dir4_wellRowID1
        └───dir1_div010
            ├───dir2_plateID0
            │   ├───dir3_wellColumnID0
            │   │   ├───dir4_wellRowID0
            │   │   └───dir4_wellRowID1
            │   └───dir3_wellColumnID1
            │       ├───dir4_wellRowID0
            │       └───dir4_wellRowID1
            └───dir2_plateID1
                ├───dir3_wellColumnID0
                │   ├───dir4_wellRowID0
                │   └───dir4_wellRowID1
                └───dir3_wellColumnID1
                    ├───dir4_wellRowID0
                    └───dir4_wellRowID1
                

### 'reportedDataset'
This directory hosts the reported dataset from Liao et al. (2022). There are two sub-directories inside: 'allFeatures_fig8' and 'highlightedFeatures_fig6fig7'. These two directories include the Jupyter files, data spreadsheet (.xlsx), and resulting .png and .svg files used to generate the violin plots found in Figures 6, 7, and 8 in Liao et al. (2022). 
The full dataset can be found here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6415474.svg)](https://doi.org/10.5281/zenodo.6415474)

