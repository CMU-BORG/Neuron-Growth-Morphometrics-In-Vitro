#####
# This runs automated neuron quantification for multiple folders of data based on
# A S Liao, W Cui, V Webster-Wood, Y J Zhang
# Quantitative evaluation of neuron developmental morphology in vitro using the change-point test
# Submitted to: Neuroinformatics 2022

# This script runs the change-point test (CPT) as described by Byrne et al. (2009) (refer to CPTautoFunc.R for more details) over
# all of the text (.txt) files in a given directory (filepath) (refer to neuronQuantFunc.R) for multiple directories. It iterates 
# over all of the directories of data based on a user-input directory path once the script is initiated.

# After the runCPTauto is complete, the script will call comboNeuriteData() (combineNeuriteDataFunc.R) to combine all of the data
# into grouped .xlsx files based on the directory the data is located in. It will create a .xlsx for each of the nested directories
# (listed below). This will include both per-neurite and per-cell metrics. (For more details, refer to combineNeuriteDataFunc.R)

# This script expects the data to be stored in the following file tree structure
# >User-Entered Directory 0
#   >Directory 1 (For this study, it was associated with the time stamp of the data set group) (note, include -div### for metadata extraction purposes)
#       >Directory 2 (For this study, it was associated with the plate ID)
#           >Directory 3 (For this study, it was associated with the row ID)
#               >Directory 4 (For this study, it was associated with the column ID)
#                   >Inside Directory 4 is where all of the data files (text files) should be located
#                   This is the directory that is inputted into the runCPTauto() function (See neuronQuantFunc.R)

##### 

# Clear All ---------------------------------------------------------------


print("Start 'runNeuronQuant.R'")

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014")

# Remove all existing R objects
rm(list=ls())

#load relevant libraries and user-defined functions

library(stringr)
scriptDir = dirname(rstudioapi::getSourceEditorContext()$path) #the directory this current script is located

source(paste(scriptDir,"/neuronQuantFunc.R",sep="")) #this file MUST be in the same filepath as this script. If it is not, must be updated to source("path_to_file/neuronQuantFUnc.R")
source( paste(scriptDir,"/combineNeuriteDataFunc.R",sep="") )

#note: neuronQuantFunc.R and combineNeuriteDataFunc.R will need the "openxlsx" and "raster" libraries as well

xFN <- "neuriteMetrics18" #base name for the exported xlsx files with all of the morphometrics


# Set Up ------------------------------------------------------------------

#Asks the user to enter a directory in the console
print("Enter the filepath where the data directories are located:")
originalFilepath <- readline()
#in the console, enter a directory (ex. G:\Shared drives\CMU BORG - Neuron Culture and Quantification\Processing_Queue\Metrics\In Process)
# Directory 0

# Run Main ----------------------------------------------------------------

#Directory 0 (User-Input)
filepath <- gsub("\\\\","/",originalFilepath) #R cannot handle '\', which is defaulted in Windows thus all of the '\' must be converted to '/'
dirInPath <- list.dirs(filepath,recursive = FALSE) #list all of the directories in the file path that was entered (Directories 1)


for (i in 1:length(dirInPath)){ #iterate through all of the directories found in Directory 0 (All of the Directories 1)
  print(paste("Folder",i,"of",length(dirInPath)," Folders in Directory 0"))
  
  plateID <- list.dirs(dirInPath[i],recursive = FALSE) #List all of the directories found in Directory 1 (i) (Directories 2)
  
  pDirsName <- list.dirs(dirInPath[i],recursive = FALSE, full.names = FALSE) #all of the directory names without the full path
  div <- str_extract(dirInPath[i],"(-div[\\d]+)")
  
  
  for (k in 1:length(plateID)){ #iterate through all of the directories found in Directory 1(i) (All of the Directories 2)
    
    rowID <- list.dirs(plateID[k],recursive = FALSE) #List all of the directories found in Directory 2(k) of Directory 1(i) (Directories 3)
    
    rDirsName <- list.dirs(plateID[k],recursive = FALSE, full.names = FALSE) #all of the directory names without the full path
    
    for (j in 1:length(rowID)){ #iterate through all of the directories found in Directory 2(k) (All of the Directories 3)
      colID <- list.dirs(rowID[j],recursive = FALSE) #List all of the directories found in Directory 3(j) of Directory 2(k) of Directory 1(i) (Directories 4)
      cDirsName <- list.dirs(rowID[j],recursive = FALSE, full.names = FALSE) #all of the directory names without the full path
      
      
      for (ii in 1:length(colID)){ #iterate through all of the directories found in Directory 3(j) (All of Directories 4)
        path2data <- colID[ii] #get the full path of Directory 4(ii)
        print(path2data)
        runCPTauto(path2data) #use the runCPTauto function on Directory 4 (ii) (neuronQuantFunc.R)
        results <- comboNeuriteData(path2data) #use the comboNeuriteData function (combineNeuriteDataFunc.R)
        
        #file name and path for the .xlsx file containing all of the data for the Directory 4 (ii) group
        xFP <- paste(dirInPath[i],paste(xFN,"_",pDirsName[k],rDirsName[j],cDirsName[ii],div,".xlsx",sep=""),sep="/")
        write.xlsx(results,file=xFP)
        
        #save the results - append if this is not the first loop
        #this will be used for a larger grouping to include all of the data in the upper Directory 3 group
        if(ii==1){
          cResults <- results
        } else {
          cResults <- Map(rbind,cResults,results)
        }
        
      }
      
      #file name and path for the .xlsx file containing all of the data for the Directory 3 (j) group
      cFP <- paste(dirInPath[i],paste(xFN,"_",pDirsName[k],rDirsName[j],div,".xlsx",sep=""),sep="/")
      write.xlsx(cResults,file=cFP)
      
      #save the results - append if this is not the first loop
      #this will be used for a larger grouping to include all of the data in the upper Directory 2 group
      if(j == 1){
        rResults <- cResults
      } else{
        rResults <- Map(rbind,rResults,cResults)
      }
      
    }
    
    #file name and path for the .xlsx file containing all of the data for the Directory 2 (k) group
    rFP <- paste(dirInPath[i],paste(xFN,"_",pDirsName[k],div,".xlsx",sep=""),sep="/")
    write.xlsx(rResults,file=rFP)
    
    #save the results - append if this is not the first loop
    #this will be used for a larger grouping to include all of the data in the upper Directory 1 group
    if(k == 1){
      pResults <- rResults
    } else{
      pResults <- Map(rbind,pResults,rResults)
    }
  }
  
  #save the results - append if this is not the first loop
  #this will be used for a larger grouping to include all of the data in the upper Directory 0 group
  if(i == 1){
    tResults <- pResults
  } else{
    tResults <- Map(rbind,tResults,pResults)
  }
  
  
  
}

#file name and path for the .xlsx file containing all of the data for the Directory 0 group (all of the data)
tFP <- paste(filepath,paste(xFN,"_all.xlsx",sep=""),sep="/")
write.xlsx(tResults,file=tFP)
