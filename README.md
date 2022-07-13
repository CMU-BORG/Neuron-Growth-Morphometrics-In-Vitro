# Neuron-Growth-Morphometrics-In-Vitro
[![DOI](https://zenodo.org/badge/473702722.svg)](https://zenodo.org/badge/latestdoi/473702722)

This is a semi-automated quantitative analysis method for assessing neuron morphology as described in "Semi-automated quantitative evaluation of neuron developmental morphology *in vitro* using the change-point test" by AS Liao, W Cui, VA Webster-Wood, and YJ Zhang (submitted to Neuroinformatics 2022). This method uses outputs from neurite traces obtained using NeuronJ, a plugin for ImageJ. This code evaluates the traces using seven morphometrics:
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
* [General Information](#general-info)
* [Directories](#directories)

## General Information

The code to conduct the Change Point Test and the subsequent morphological calculations based on neurite traces can be found in the 'analysisCode' directory. This tool requires the trace coordinate outputs (which are output from NeuronJ as a .txt file) and an associated .csv file for the neurite identification. In addition, metadata on the time point and the well ID that the trace data belongs to are present in the organization of directories that house the data, as seen in the 'unprocessed' directory of 'exampleDataset'. The neurite trace coordinate .txt files and the .csv files need to be the innermost directory (i.e. 'row_0' and 'row_1'). To use the tool, run the runNeuronQuant.R script. The user will be prompted to enter a directory that houses the data, which should be the path for the 'experiment_experimentID' directory. Expected outputs can be found in the 'processed' directory. 

In addition, we have also included the R script for our follow up statistical analysis (RStatTest.R) to compare the morphometrics at different time points. This script was developed to run the Anderson-Darling test to assess normality, the Kruskal-Wallis test, and the Dunn test with a Bonferroni correction for our organized dataset. This script has not been generalized and is specific for our datasheet (perCellMetrics(den10,rem8,remMiss).csv located in 'reportedDataset/allFeatures_fig8'), but was included for completeness in detailing the results of our manuscript. The results from this script were reported in Tables 2-16 in Liao et al. (2022).

## Directories
### 'analysisCode'
This directory includes the R code used to automatically quantitatively assess the morphometrics of the neurites. Run the script named runNeuronQuant.R to perform the complete analysis.

### 'exampleDataset'
This directory includes an example dataset to run the analysis code (it is a smaller subset of the data analyzed in Liao et al. (2022)). There are two sub-directories inside: 'unprocessed' and 'processed'. The 'unprocessed' directory includes the test data files in the appropriate file tree structure that runNeuronQuant.R expects. The 'processed' directory includes the expected outputs after running 'runNeuronQuant.R' and inputting the path to 'experiment_001'.

The directory tree structure can be found below. This directory tree has been designed to support the analysis of datasets that may include multiple experiments with individual IDs (experimentID), multiple time points (div###), multiple well plates with individual IDs (plateID), multiple columns per plate (columnID), and multiple wells within a given column (rowID). This tree structure was designed to be as flexible as possible for describing experiments using well plates. Alternative tree structures may be used but will require that the user update the code appropriately. As implemented, the outermost directory is 'experiment_{experimentID}', which should house all the data from a given experiment. Within the experiment, directory sub-folders for each time point are included ('timepoint_div{###}') to organize the data into the time points in which they were imaged. Within each time point folder, plate subfolders 'plate_{plateID}' are used to group the data by the well plate ID, which would be useful if multiple plates were used simultaneously. Each well plate folder contains folders organized by the column within the well plate, the 'column_{columnID}' and subsequently organized by the row in each column 'row_{rowID}'. This enables multiple groups of data to be analyzed with a single run of 'runNeuronQuant.R'.

As an example, the experimental directory tree for one experiment with two time points, two plates, and samples in columns 0 and 1 and in rows 0 and 1 found in the subfolers within the exampleDataset directory appears as follows:

    └───experiment_001
        ├───timepoint_div000
        │   ├───plate_0
        │   │   ├───column_0
        │   │   │   ├───row_0
        │   │   │   └───row_1
        │   │   └───column_1
        │   │       ├───row_0
        │   │       └───row_1
        │   └───plate_1
        │       ├───column_0
        │       │   ├───row_0
        │       │   └───row_1
        │       └───column_1
        │           ├───row_0
        │           └───row_1
        └───timepoint_div010
            ├───plate_0
            │   ├───column_0
            │   │   ├───row_0
            │   │   └───row_1
            │   └───column_1
            │       ├───row_0
            │       └───row_1
            └───plate_1
                ├───column_0
                │   ├───row_0
                │   └───row_1
                └───column_1
                    ├───row_0
                    └───row_1
                

### 'reportedDataset'
This directory hosts the reported dataset from Liao et al. (2022). There are two sub-directories inside: 'allFeatures_fig8' and 'highlightedFeatures_fig6fig7'. These two directories include the Jupyter files, data spreadsheet (.xlsx), and resulting .png and .svg files used to generate the violin plots found in Figures 6, 7, and 8 in Liao et al. (2022). 
The full dataset can be found here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6415474.svg)](https://doi.org/10.5281/zenodo.6415474). For a complete description of the dataset, please refer to the 'README.txt' file that is included with the full dataset. Briefly, the dataset contains the unprocessed raw images, the neurite traces from those images and the processed data that includes the results from the morphometric and statistical analyses. The unprocessed data can be used with the software tools presented in this GitHub repository to produce the results in the processedData folder. 

