#####
# This performs the automated neuron quantification based on
# A S Liao, W Cui, V Webster-Wood, Y J Zhang
# Quantitative evaluation of neuron developmental morphology in vitro using the change-point test
# Submitted to: Neuroinformatics 2022

# This script takes the results from NeuronJ (coordinates of the trace, tracing metadata) and the change-point test results
# (change-point coordinates) to calculate the morphometrics described in Liao et al. (2022):
# -total length per cell (sum of all neurite lengths)
# -average tortuosity per cell
# -number of end points per cell (degree)
# -number of neurites per cell
# -number of change points per cell
# -average segment length (distance between change points) per cell
# -average relative turning angle between change points per cell

# This script MUST be run AFTER using the neuronQuantFunc.R. It requires the output: 
#     (1) <imageName>.N<traceID>.xlsx (from neuronQuantFunc.R (runCPTauto))
#     (2) <imageName>_tracings.csv (from NeuronJ export of "Tracings")
#     (3) <imageName>.N<traceID>.txt (from NeuronJ export of trace coordinates)

# INPUTS
#   (i) dirpath - path to a directory that contains the data files described above
# 
# OUTPUTS 
#   (i) results (list) of the morphometrics calculated
#       (a) ang - contains the relative angles
#       (b) dist - contains the total and segment lengths
#       (c) numCPs - the number of change points per neurite
#       (d) perCellMetrics - the angles, lengths, number of change points, average tortuosity, degree, and number of neurites per cell
#       (e) tortPerImage - tortuosity per neurite
#       (f) noCPTraceFiles - all of the files that did not have any change points
#       


comboNeuriteData <- function(dirpath){
  
  #necessary libraries
  library(openxlsx)
  library(raster)
  
  #
  
  
  # Look for all of the files with the .csv file extension
  # Each image should only have 1 of these since they should have been directly exported from NeuronJ
  imageCSV <- list.files(dirpath,pattern="*.csv",full.names=TRUE) #full paths
  imageCSV_only <- list.files(dirpath,pattern="*.csv",full.names=FALSE) #file names only
  
  # Look for all of the trace file results from runCPTauto (should only be 1 per neurite that has change points)
  traceX <- list.files(dirpath,pattern="*(.N[:digit:]{1,}.xlsx)",full.names=TRUE) #full paths
  traceX_only <- list.files(dirpath,pattern="*(.N[:digit:]{1,}.xlsx)",full.names=FALSE) #file names only
  
  #preallocate results variables as empty data frames
  ang <- data.frame()
  dist <- data.frame()
  numCPs <- data.frame()
  perCellMetrics <- data.frame()
  tortPerImage <- data.frame()
  noCPTraceFiles <- c()
  
  counter <- 0
  
  #iterate over all of the possible images
  for(i in 1:length(imageCSV)){
    
    #for a given image, list all of the associated trace .xlsx files (output from runCPTauto)
    traceImage <- list.files(dirpath, pattern=str_c("(",str_remove(imageCSV_only[i],"_tracings.csv"),".N\\d+.xlsx)"),full.names=TRUE) #full paths
    traceImage_only <- list.files(dirpath, pattern=str_c("(",str_remove(imageCSV_only[i],"_tracings.csv"),".N\\d+.xlsx)"),full.names=FALSE) #file names only
    
    #read the *_tracings.csv for the given image - this provides metadata about all of the traces in that image
    traceInfo <- read.csv(imageCSV[i],header=TRUE)
    
    #extract the metadata for all of the tracings from the *_tracings.csv
    traceLabel <- str_extract(as.character(traceInfo$Tracing),"(N[\\d]+)")
    cellLabel <- str_extract(as.character(traceInfo$Label),"([\\d]+)")
    div <- as.numeric(str_extract(str_extract(dirpath,"(-div[\\d]+)"),"([\\d]+)"))/10
    plateID <- substr(str_extract(imageCSV_only[i],"([:upper:][:upper:][:digit:])"),1,1)
    rowID <- substr(str_extract(imageCSV_only[i],"([:upper:][:upper:][:digit:])"),2,2)
    colID <- as.numeric(substr(str_extract(imageCSV_only[i],"([:upper:][:upper:][:digit:])"),3,3))
    
    #the column ID indicates the density in which the cells were cultured at
    if (colID == 2){
      density = 10000 #cells/cm2
    } else if (colID == 4){
      density = 20000 #cells/cm2
    } else if (colID == 6){
      density = 60000 #cells/cm2
    } else {
      density = -1 #no valid colID listed, return as -1
    }
    
    #preallocate empty data frames for calculations per image
    angPerImage <- data.frame()
    distPerImage <- data.frame()
    cpPerImage <- data.frame()
    
    #list all of the text files that contain the trace coordinates (*.N<traceID>.txt)
    cFN <- list.files(dirpath,pattern=paste(str_remove_all(imageCSV_only[i],"_tracings.csv"),".N.*\\.txt$",sep=""),full.names=TRUE)
    
    #preallocate vectors for calculations
    tracing <- c(1:length(cFN))*0
    tortuosity <- c(1:length(cFN))*0
    imageID2 <- c(1:length(cFN))*0
    cellID2 <- c(1:length(cFN))*0
    div2 <- c(1:length(cFN))*0
    plateID2 <- c(1:length(cFN))*0
    rowID2 <- c(1:length(cFN))*0
    colID2 <- c(1:length(cFN))*0
    density2 <- c(1:length(cFN))*0
    
    #loop for each trace coordinate text file
    for (jj in 1:length(cFN)){
      
      print(paste("Trace Coordinate File",jj,"of",length(cFN)))
      print(cFN[jj])
      
      #read in the coordinate file
      cData <- scan(cFN[jj],list(y=0,x=0))
      #save trace coordinate data as xy
      xy <- cbind(cData$x,cData$y)
      
      #calculate the neurite length
      totPathLen <- sum(pointDistance(xy[1:(length(cData$y)-1), ],xy[-1,],lonlat = FALSE))
      #calculate the distance between the neurite's start and end points
      endPtDist <- pointDistance(xy[1, ], xy[length(cData$y),], lonlat=FALSE)
      #calculate tortuosity (neurite path length divided by the distance between the start and end points)
      tortuosity[jj] <- totPathLen/endPtDist
      
      #associate this with the relevant metadata for this particular trace - this will be saved as 1 row
      tracing[jj] <- str_extract(cFN[jj],"([N][\\d]+)")
      cellID2[jj] <- cellLabel[traceInfo$Tracing == tracing[jj]]
      imageID2[jj] <- str_remove_all(imageCSV_only[i],"_tracings.csv")
      div2[jj] <- as.numeric(str_extract(str_extract(dirpath,"(-div[\\d]+)"),"([\\d]+)"))/10
      plateID2[jj] <- substr(str_extract(cFN[jj],"([:upper:][:upper:][:digit:])"),1,1)
      rowID2[jj] <- substr(str_extract(cFN[jj],"([:upper:][:upper:][:digit:])"),2,2)
      colID2[jj] <- as.numeric(substr(str_extract(cFN[jj],"([:upper:][:upper:][:digit:])"),3,3))
      density2[jj] <- density
      
    }
    #save all of the tortuosity data calculated as one data frame
    tortPerTracing <- data.frame(imageID2, div2, plateID2, rowID2, colID2, density2, cellID2, tracing, tortuosity)
    tortPerImage <- rbind(tortPerImage,tortPerTracing)
    
    
    
    if(length(traceImage)==0){
      #do nothing - there are no files
      counter = counter+1 
      noCPTraceFiles[counter] <- imageCSV[i] #files without change points
    } else{
      for(j in 1:length(traceImage)){
        #read in the relative angles and lengths calculated by neuronQuantFunc.R (runCPTauto)
        angTemp <- read.xlsx(traceImage[j],1)
        distTemp <- read.xlsx(traceImage[j],2)
        
        cpSheet <- read.xlsx(traceImage[j],5)
        
        #metadata for the given trace
        traceID <- str_extract(traceImage_only[j],"(N[:digit:]{1,})")
        cellID <- traceInfo$Label[traceInfo$Tracing == traceID]
        
        #save the relevant data from the xlsx file into the current data frames
        angTemp$angCPsRDd180ABS <- abs(angTemp$angCPsRDd180)
        angTemp$traceID <- traceID
        angTemp$cellID <- cellID
        angTemp$div <- div
        angTemp$plateID <- plateID
        angTemp$rowID <- rowID
        angTemp$colID <- colID
        angTemp$density <- density
        angTemp$imageID <- str_remove_all(imageCSV_only[i],"_tracings.csv")
        angPerImage <- rbind(angPerImage, angTemp)
        
        distTemp$traceID <- traceID
        distTemp$cellID <- cellID
        distTemp$div <- div
        distTemp$plateID <- plateID
        distTemp$rowID <- rowID
        distTemp$colID <- colID
        distTemp$density <- density
        distTemp$imageID <- str_remove_all(imageCSV_only[i],"_tracings.csv")
        distPerImage <- rbind(distPerImage, distTemp)
        
        cpTemp <- data.frame(imageID = str_remove_all(imageCSV_only[i],"_tracings.csv"), traceID = traceID, cellID = cellID, div = div, plateID = plateID, rowID = rowID, colID = colID, density = density, numCPs = nrow(cpSheet))

        cpPerImage <- rbind(cpPerImage, cpTemp)
      }
      uniqueCells <- unique(angPerImage$cellID) #number of unique cells in the image
      
      
      perCellMetricsPerImage <- data.frame() #preallocate
      
      #the following will calculate the morphometrics as "per cell" rather than "per neurite"
      #loop through the number of unique cells per image
      for(k in 1:length(uniqueCells)){
        #average turning angle per cell
        avgRelTurnAngle <- mean((angPerImage$angCPsRDd180[angPerImage$cellID == uniqueCells[k]]),na.rm=TRUE)
        #absolute value of the average turning angle
        avgAbsRelTurnAngle <- mean(abs(angPerImage$angCPsRDd180[angPerImage$cellID == uniqueCells[k]]),na.rm=TRUE)
        
        #average segment length per cell
        avgSegLen <- mean(distPerImage$distSegUM[distPerImage$cellID == uniqueCells[k]])
        
        #total length (sum of all neurite lengths) per cell
        totLen <- sum(distPerImage$totDistUM[distPerImage$cellID == uniqueCells[k]],na.rm=TRUE)
        
        #number of end points per cell
        degree <- sum(cellLabel == uniqueCells[k])
        
        #number of neurites per cell - calculated by the number of "primary" neurites from the NeuronJ metadata
        numNeurites <- sum(traceInfo$Type[cellLabel == uniqueCells[k]] == "Primary")
        
        #average tortuosity per cell
        avgTort <- mean((tortuosity[cellID2 == uniqueCells[k]]),na.rm=TRUE)
        
        #total number of change points per cell
        totCPs <- sum(cpPerImage$numCPs[cpPerImage$cellID == uniqueCells[k]])
        
        #average number of change points per cell
        avgCPs <- mean(cpPerImage$numCPs[cpPerImage$cellID == uniqueCells[k]])
        
        #save all of the per cell metrics with relevant metadata in one data frame
        metricsTemp <- data.frame(imageCSV_only[i],div, plateID, rowID, colID,density,uniqueCells[k],avgRelTurnAngle,avgAbsRelTurnAngle,avgSegLen,totLen, degree, numNeurites, avgTort, totCPs, avgCPs)
        perCellMetricsPerImage <- rbind(perCellMetricsPerImage,metricsTemp)
      }
      
      #save all of the ang, dist, perCellMetrics, and numCPs per loop in one data frame
      ang <- rbind(ang,angPerImage)
      dist <- rbind(dist, distPerImage)
      perCellMetrics <- rbind(perCellMetrics,perCellMetricsPerImage)
      numCPs <- rbind(numCPs,cpPerImage)

    }
    
    
  }
  
  #return the morphometrics results as one list of data frames
  results <- list(ang = ang, dist = dist, numCPs = numCPs, perCellMetrics = perCellMetrics, tortPerImage = tortPerImage, noCPTraceFiles = data.frame(noCPTraceFiles))
  
  return(results)
  
}