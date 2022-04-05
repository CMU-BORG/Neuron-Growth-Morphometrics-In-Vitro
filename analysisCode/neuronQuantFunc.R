#####
# This performs the automated neuron quantification based on
# A S Liao, W Cui, V Webster-Wood, Y J Zhang
# Quantitative evaluation of neuron developmental morphology in vitro using the change-point test
# Submitted to: Neuroinformatics 2022

# This script runs the change-point test (CPT) as described by Byrne et al. (2009) (refer to CPTautoFunc.R for more details) over
# all of the text (.txt) files in a given directory (filepath). It iterates over varying q-values (1 to 10, integers) and saves the
# CPT results that had the largest number of change-points with the smallest q-value. After the CPT iterations, the script calculates
# the length of the neurite, the segment length (distance between change points), and the relative turning angle between change points.

# INPUTS
#   (i) filepath - path to a directory that contains text files (*.N#.txt) of the data (expected to be a 2 column list of coordinates)
#                  the file name must also include either "_20_" or "_40_" for the magnification - otherwise pixel values will be returned rather than real world units
# 
# OUTPUTS (none to Environment)
#   (i) in the same directory as filepath, 2 types of .xlsx files should be saved:
#   
#       (a) filesWithNoCP.xlsx: a list of data file names that did not return any change points 
#                               (it is likely these contained too few coordinate points: <= 3 points (refer to CPTautoFunc.R))
# 
#       (b) <dataFileName>.xlsx: this is the results from the change-point test (CPTautoFunc.R), 
#                                resulting CPT-based morphometrics (relative angle and segment length) and the length of the neurite
#   

#####
runCPTauto <- function(filepath){
  set.seed(0) #the CPTautoFunc introduces some randomness in its algorithm. By setting the seed to a constant value, this ensures reproducibility if the analysis is rerun
  print("Starting runCPTauto")
  
  #Load in necessary libraries and function
  library(openxlsx) #load in openxlsx library - required for saving files to a .xlsx file
  library(raster) #load in raster library - required for pointDistance function
  
  scriptDir = dirname(rstudioapi::getSourceEditorContext()$path); #the directory this current script is located
  #Load Function CPTauto from File
  source(paste(scriptDir,"/CPTautoFunc.R",sep="")) #this file MUST be in the same filepath as this script. If it is not, must be updated to source("path_to_file/CPTautoFunc.R")
  
  
  #Look for all of the .txt files in the directory "filepath"
  txtFiles <- list.files(filepath,pattern="*.txt")
  filenameMaster <- sapply(strsplit(txtFiles,split='.t',fixed=TRUE),function(x) (x[1])) #hold only the file names without the .txt extension
  fileext <- "txt"
  
  #user set q and alpha values (refer to Byrne et al. (2009) and the CPTauto function)
  q <- 1:10 #iterate over q=1 to q=10
  alpha <- 0.05
  
  
  #####
  #Run the Change Point Test (CPT)
  
  #number of q and alpha values entered
  numq <- length(q)
  numa <- length(alpha)
  numf <- length(txtFiles)
  g <- 0 #counter for number of files that have no change points
  noCPfiles <- "All files have change points" #if all files have change points, this will be the default message
  
  #loop for all of the text files found
  for (a in 1:numf){
    
    filename <- filenameMaster[a] #file name for the given iteration
    
    #print to console which file is currently being run
    print(paste("number of text files: ",a," in ",numf))
    print(paste("filename: ",filename))
    
    #looks for "N#" in the file name to identify which trace - if this is not found, the code will NOT process the text file
    if (str_detect(filename,"([N][\\d]+)")){
      saveDataHere <- paste(filepath,filename,sep="/")
      #preallocate a matrix (numCP) with zeros
      #this will hold the number of change points (CPs) for every combination of q and alpha
      
      numCP <- matrix(0,numq,numa)
      
      #preallocate vectors to hold row and column names (q = number, alpha = number)
      coln <- c(1:numa)
      rown <- c(1:numq)
      maxCPs <- 0
      maxNumCPs <- 0
      maxq <- 0
      maxa <- 0
      
      #this nested for loop will run the data through the CPTauto function with different q and alpha
      for (i in 1:numa){
        print(paste("alpha =",alpha)) #print the i value to console - for debugging
        
        coln[i] <- paste("alpha = ",alpha[i]) #alpha name
        
        for (j in 1:numq){
          print(paste("q =",q[j])) #print j value to console - for debugging
          
          rown[j] <- paste("q = ",q[j]) #q name
          
          results <- CPTauto(saveDataHere,fileext,q[j],alpha[i]) #run the data through the CPTauto function
          
          dimCPS <- dim(results$cps) #number of CPs
          
          numCP[j,i] <- dimCPS[1] #save the number of CPs to the appropriate spot in numCP
          
          #We are searching for the max number of CPs. Save the data if it is the max CP iteration to be used later
          if(dimCPS[1]==max(numCP) & (dimCPS[1]!=maxNumCPs) & (length(results$cps)>1)){
            maxCPs <- results$cps #coordinates corresponding to the maximum number of CPs
            maxNumCPs <- dimCPS[1] #maximum number of CPs
            maxq <- q[j] #q corresponding with the max number of CPs
            maxa <- alpha[i] #alpha corresponding with the max number of CPs
            Y <- results$x1 #Y COORDINATES
            X <- results$x2 #X COORDINATES
          }
          
        }
      }
      
      print("CPT iteration complete") #debugging, print to console
      
      #after iterating, if there was more than one CP found, this is considered to be a valid run so the data will be saved
      if(length(maxCPs) > 1){
        #Attach the row and column names to the numCP to show which alpha and q values each numCP value belongs to
        rownames(numCP) <- rown
        colnames(numCP) <- coln
        
        numCPplotName <- paste(saveDataHere,"CPplot.png",sep="_")
        png(numCPplotName)
        plot(q,numCP) #plot q vs numCP
        
        #dev.copy(png,numCPplotName) #for debugging - see plot to determine which q corresponds with maxCP
        dev.off()
        
        coords <- cbind(X,Y) #put the neurite coordinates into one matrix (Y,X)
        
        #####
        #Calculate the Neurite Lengths & Angles
        
        
        print("Neurite length, segment length, and relative angles calculations begin") #print to console, debugging
        
        
        
        diffCoordY <- Y[-1]-Y[1:(length(Y)-1)] #difference in the sequential Y coordinates
        diffCoordX <- X[-1]-X[1:(length(X)-1)] #difference in the sequential X coordinates
        
        angCoords <- atan2(diffCoordY,diffCoordX) #calculate the angle between each sequential coordinate
        #4-quadrant angle in radians
        
        angSeg <- (1:maxNumCPs)*0 #preallocate vector with zeroes, angSeg will hold the average angle based on each coordinate between 2 sequential CPs
        
        distCoords <- pointDistance(coords[1:(length(Y)-1), ],coords[-1,],lonlat = FALSE) #length between each sequential coordinate
        totalDist1 <- sum(distCoords) #total neurite length
        
        distSeg <- (1:maxNumCPs)*0 #preallocate vector for the sum of coordinate lengths to determine length of each segment (between 2 sequential CPs)
        oldIndex <- 1; #first index used for the for loop (Index the initial coordinate point)
        totalDist2 <- 0 #variable to sum the total distance (check with totalDist1 for debugging)
        
        #This for loop will calculate the length of each segment (between 2 sequential CPs) and the average coordinate angle
        for (k in 1:(maxNumCPs+1)){
          if (k != (maxNumCPs+1)){ #for every iteration except the last one...
            index <- which((Y == maxCPs[k,3]) & (X == maxCPs[k,4])) #find the index of the CP
            distSeg[k] <- sum(distCoords[oldIndex:(index-1)]) #sum the distance between the previous index (oldIndex) up to the current one (Index)
            angSeg[k] <- mean(angCoords[oldIndex:(index-1)]) #average coordinate angle (same indexing as above)
            oldIndex <- index #update the oldIndex for next iteration
          } else { #for the final iteration, use the final coordinate point with the last CP (held in oldIndex)
            distSeg[k] <- sum(distCoords[oldIndex:length(distCoords)])
            angSeg[k] <- mean(angCoords[oldIndex:length(angCoords)])
          }
          totalDist2 <- totalDist2 + distSeg[k] #add to the total distance
          
        }
        maxCPsXY <- cbind(maxCPs[,4],maxCPs[,3])
        startCPend <- rbind(coords[1,],maxCPsXY,coords[length(Y),]) #beginning coordinate, change points, end coordinate
        
        #distance between change points
        distCP <- pointDistance(startCPend[1:(nrow(startCPend)-1),],startCPend[-1,], lonlat = FALSE)
        totalDist3 <- sum(distCP) #total neurite distance based on distCP
        
        
        #diffCPs <- rbind(maxCPs[,3:4],coords[length(Y),])-rbind(coords[1,],maxCPs[,3:4]) #difference in sequential CPs & the initial and final coordinate points
        diffCPs <- startCPend[-1,] - startCPend[1:(nrow(startCPend)-1),]
        angCPs <- atan2(diffCPs[,2],diffCPs[,1]) #angle created by the CPs
        angCPsRelDiff <- angCPs[-1] - angCPs[1:(length(angCPs)-1)] #difference in angle relative to the previous segment angle
        angCPsRD <- c(NA,angCPsRelDiff)
        
        #Convert to Degrees
        rad2deg <- 180/pi #180 degrees per pi radians
        angSegd <- angSeg*rad2deg
        angCPsd <- angCPs*rad2deg
        angCPsRelDiffd <- angCPsRelDiff*rad2deg
        angCPsRDd <- angCPsRD*rad2deg
        
        ind360p <- which(angCPsRDd > 180)
        ind360n <- which(angCPsRDd < (-180))
        
        angCPsRDd180 <- angCPsRDd
        
        angCPsRDd180[ind360p] <- angCPsRDd180[ind360p] - 360
        angCPsRDd180[ind360n] <- angCPsRDd180[ind360n] + 360
        
        allAng <- cbind(angSegd,angCPsd,angCPsRDd,angCPsRDd180)
        
        #Look for the magnification that was used for imaging in the file name
        if (str_detect(filename,"([_][\\d][\\d][_])")){
          mag <- as.numeric(gsub("[^[:digit:].]","", str_extract(filename,"([_][\\d][\\d][_])")))
        } else if (str_detect(filename,"([-][\\d][\\d][-])")){
          mag <- as.numeric(gsub("[^[:digit:].]","", str_extract(filename,"([-][\\d][\\d][-])")))
        } else if (str_detect(filename,"([-][\\d][\\d][_])")){
          mag <- as.numeric(gsub("[^[:digit:].]","", str_extract(filename,"([-][\\d][\\d][_])")))
        } else if (str_detect(filename,"([_][\\d][\\d][-])")){
          mag <- as.numeric(gsub("[^[:digit:].]","", str_extract(filename,"([_][\\d][\\d][-])")))
        } else{
          mag = -1 #no valid magnification value could be extracted from the file name
        }
        
        #Convert from Pixels to um (this ratio is used based on default scale bars from the Echo Revolve for 20X and 40X objectives)
        if (mag==20){
          #For 20X, it is 711 pixels for 180 um
          pix2um <- 180/711
          unit <- "um"
        } else if (mag == 40){
          #For 40X, it is 711 pixels for 90 um
          pix2um <- 90/711
          unit <- "um"
        } else{
          print(paste("mag =",mag))
          print("no valid magnification value found (20 or 40), returning lengths in pixels")
          pix2um <- 1
          unit <- "pixels"
        }
        
        #convert from pixel units to real world units (unless invalid magnification - then pixel units will be maintained)
        totDistUM <- totalDist1*pix2um
        distSegUM <- distSeg*pix2um
        distCPUM <- distCP*pix2um
        totDistCPUM <- totalDist3*pix2um
        
        totDistUM_resized <- rbind(cbind(totDistUM,totDistCPUM),matrix(NA,(length(distSegUM)-1),2))
        
        
        allDist <- cbind(distSegUM,distCPUM,totDistUM_resized)
        
        inputs <- cbind(maxa,maxq,mag,pix2um)
        
        
        saveDataHereExt <- paste(saveDataHere,"xlsx",sep=".")
        
        saveData <- list("allAng"=allAng,"allDist"<-allDist,"distUnit"=unit,"coords"=coords,"CPs"=maxCPs,"inputs"=inputs)
        write.xlsx(saveData,file=saveDataHereExt)
        
        
      }else{ #there are no change points found for this file
        
        if(g == 0){ #if this is the 1st file that was identified to have no change points, begin the 'noCPfiles' list
          noCPfiles <- filename; #change the content of noCPfiles to be the filename of the file without any change points
          
        }else{
          noCPfiles <- rbind(noCPfiles,filename) #add to the list of the files with no change points
        }
        g <- g + 1 #counter for number of files without any change points
      }
    }
    
  }
  
  saveNoCP <- paste(filepath,"filesWithNoCP.xlsx",sep="/")
  write.xlsx(noCPfiles,file = saveNoCP)
  
  
  print("End of runCPTauto") #print to console, debugging
  
}