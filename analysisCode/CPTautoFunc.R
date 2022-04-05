#####
# R code for performing the "circular change point test" (CPT)
# as described in "How did they get here from there?
# Detecting changes of direction in terrestrial ranging."
# by R W Byrne, R G Noser, L A Bates & P E Jupp. 
# Animal Behaviour 77, 619-631, 2009.
# A manual for this code can be found at 
# http://www.mcs.st-andrews.ac.uk/~pej/CPTauto_Rcode_Manual.pdf.

# P. E. Jupp 26 June 2012
# Code supplied without guarantee

#####
# The CPT Rcode from Byrne et al. (2009) was modified for 
# A S Liao, W Cui, V Webster-Wood, Y J Zhang
# Quantitative evaluation of neuron developmental morphology in vitro using the change-point test
# Submitted to: Neuroinformatics 2022
#
# This is an "automated" version in which 
# (a) q  and  alpha  are set by the user;
#       (i)  alpha = significance level; alpha = 0.05 is a convenient default
#       (ii) The user should replace this by a value of q that is appropriate for the species and data 
#            considered, e.g. that obtained by running CPT_Rcode on part of the data. 
# (b) the parameters N (number of permutations) and
# tol (maximum distance between indistinguishable 
# positions) may be changed by the user manually in the code itself; 

# INPUTS
#     (i)   filetot - the filename of the file with the data (coordinates) - includes full file path
#     (ii)  fileext - the file extension of the data file (typically "txt")
#     (iii) q - q parameter (refer to Byrne et al (2009))
#     (iv)  alpha - significance level (refer to Byrne et al (2009))

# OUTPUTS:
#     (i) (to Environment) a list (in cps) of (row numbers and 
#         coordinates of) waypoints that are detected as  significant
#     (ii) (to directory) a plot of the path on which change points are indicated.

# (Code for the original version is available at http://www.mcs.st-andrews.ac.uk/~pej/CPT_Rcode.)


CPTauto <- function(filetot,fileext,q,alpha){
  filetottxt <- paste(filetot,fileext,sep = ".")
  filetotpng <- paste(filetot,"_q",q,"_a",alpha,".png",sep = "")
  
  # Read in data file
  inp<-scan(filetottxt,list(x1=0,x2=0))
  x1<-inp[[2]]
  x2<-inp[[1]]
  x1 <- (max(x1)+1)-x1
  x2 <- x2 - min(x2)
  # Inspect first few rows of the data
  xy <- cbind(x1,x2)
  head(xy)
  
  if (length(inp$x1) > 3){ #check that the number of data points is greater than 3 (not too short, check for sufficient amount of data)
    # input N 
    # N = total number of permutations (1 observed and N-1 simulated)
    # N = 10000 is a convenient number
    N <- 10000
    
    # input tol 
    # tol = tolerance 
    # = maximum distance between indistinguishable positions
    # tol = 0 is a convenient default
    # The user may like to replace this by 
    # twice the average GPS error in research area
    tol <- 0
    
    
    # PRELIMINARIES
    
    # Reverse the time-ordering,
    # so that (bx1[1], bx2[1]) refers to (final) 
    #  putative goal
    
    bx1 <- 0*x1
    bx2 <- 0*x2
    
    n <- length(x1)
    for (j in 1:n){
      bx1[j] <- x1[n-j+1]
      bx2[j] <- x2[n-j+1]
    }
    
    
    # Calculate the steps (bxdiff1, bxdiff2)
    bx1diff <- diff(bx1)
    bx2diff <- diff(bx2)
    
    
    # REMOVE POINTS AT WHICH ANIMAL STAYS STILL
    # ind is (reverse) time ordering of points
    # newp  > 0 if point differs from previous point
    tolsq <- tol^2
    ind <- c(1:n)
    newp <- c(1:n)
    for (j in 2:n){
      newp[j] <- ind[j]*( bx1diff[j-1]^2 + bx2diff[j-1]^2 > tolsq)
    }
    # bz1, bz2 are coordinates of points(in reverse time order) at which there is movement 
    bz1 <- bx1[newp >0]
    bz2 <- bx2[newp >0]
    nz <- length(bz1)
    # We shall apply the CPT to points with coordinates (bz1,bz2)
    
    # Calculate the steps (bzdiff1, bzdiff2)
    bz1diff <- diff(bz1)
    bz2diff <- diff(bz2)
    
    
    # SOME DECLARATIONS
    
    # goal.no = number (from end) of current putative goal (with goal.no = 1 for end position 
    # Thus goal.no = goal.t + 1 
    # where goal.t was used in original code to refer to 
    # time of current putative goal (with goal.t = 0 for # end position 
    # start with goal.no = 1 
    goal.no <- 1
    
    # last.no = number (backwards in time) of last position of interest
    # last.no = length(bz1) is a convenient default
    # (and includes all the positions)
    last.no <- length(bz1) 
    
    no.of.nos <- length(bz1)
    
    Rsumrand  <- rep(0, len=N)
    Rsumrandr <- rep(0, len=N)
    
    # Pr will store observed p-values in a run of r
    Pr <- rep(1, len= no.of.nos)
    
    # sig is a vector indicating whether or not 
    # waypoint k is detected as a possible change point:
    # sig[k] = 1  if change point detected at waypoint k
    # sig[k] = 0  otherwise
    sig <- rep(0,nz)
    
    
    # k is number of steps in "k-leg"
    k <- 0
    
    # LOOK (SEQUENTIALLY) FOR NEXT POSSIBLE CHANGE POINT
    
    while(goal.no < last.no - q){
      
      k <- 0 
      
      # INCREASE k UNTIL NEXT POSSIBLE CHANGE POINT IS FOUND
      
      # Re-initialise Pr 
      Pr <- 0*Pr + 1
      
      P <- 1
      
      while ((P > alpha) && (goal.no + q + k < last.no)){ 
        # B
        
        k <- k + 1
        
        R1 <- sqrt((bz1[goal.no+k] - bz1[goal.no])^2 + (bz2[goal.no+k] - bz2[goal.no])^2)
        R2 <- sqrt((bz1[goal.no+k+q] - bz1[goal.no+k])^2 + (bz2[goal.no+k+q] - bz2[goal.no+k])^2)
        
        Rsum <- R1 + R2
        
        # Rsumrand[1] = observed value of statistic R1 + R2
        Rsumrand[1] <- Rsum
        
        # Now calculate statistic R1 + R2 for a further N-1 random permutations
        # and store in Rsumrand
        for (it in 2:N){ 
          # start of it loop C
          u <- runif(k+q,0,1)
          perm <- order(u)
          bz1r <- bz1[goal.no]
          bz2r <- bz2[goal.no]
          for (j in 1:k){ # start of j loop D
            bz1r <- bz1r + bz1diff[goal.no-1+perm[j]]
            bz2r <- bz2r + bz2diff[goal.no-1+perm[j]]
          } # end of j loop D
          
          R1rand <- sqrt((bz1r - bz1[goal.no])^2 + (bz2r  - bz2[goal.no])^2)
          R2rand <- sqrt((bz1[goal.no+k+q] - bz1r)^2 + (bz2[goal.no+k+q] - bz2r)^2)
          Rsumrand[it] <- R1rand + R2rand
        } # end of it loop C
        
        P <- sum(Rsumrand >= Rsum)/N
        
        Pr[k] <- P
        
      }  # end of "while" on 
      #  ((P > alpha) && (goal.no + q + k < last.no)) B
      
      # f is the first value of k in the current run that is significant
      f <- k 
      
      while ((P <= alpha) && (goal.no + q + k < last.no)){
        # E
        
        k <- k+1
        R1 <- sqrt((bz1[goal.no+k] - bz1[goal.no])^2 + (bz2[goal.no+k] - bz2[goal.no])^2)
        R2 <- sqrt((bz1[goal.no+k+q] - bz1[goal.no+k])^2 + (bz2[goal.no+k+q] - bz2[goal.no+k])^2)
        
        
        Rsum <- R1 + R2
        
        # Rsumrand[1] = observed value of statistic R1 + R2
        Rsumrand[1] <- Rsum
        
        # Now calculate statistic R1 + R2 for a further N-1 random permutations
        # and store in Rsumrand
        for (it in 2:N){ # start of it loop F
          u <- runif(k+q,0,1)
          perm <- order(u)
          bz1r <- bz1[goal.no]
          bz2r <- bz2[goal.no]
          for (j in 1:k){ # start of j loop G
            bz1r <- bz1r + bz1diff[goal.no-1+perm[j]]
            bz2r <- bz2r + bz2diff[goal.no-1+perm[j]]
          } # end of j loop G
          
          R1rand <- sqrt((bz1r - bz1[goal.no])^2 + (bz2r  - bz2[goal.no])^2)
          R2rand <- sqrt((bz1[goal.no+k+q] - bz1r)^2 + (bz2[goal.no+k+q] - bz2r)^2)
          Rsumrand[it] <- R1rand + R2rand
        } # end of it loop F
        
        P <- sum(Rsumrand >= Rsum)/N
        
        Pr[k] <- P
        
      }  # end of "while" on 
      #  ((P < alpha) && (goal.no + q + k < last.no)) E
      
      # l is the last value of k in the current run that is significant
      l <- k - 1
      
      
      # apply "peak rule"
      Pr[f:l]
      rmin <- min(which(Pr[f:l] == min(Pr[f:l]))) 
      rmin <- if(l >=f) rmin else 0
      
      goal.no <- ifelse(rmin > 0, goal.no + f + rmin - 1, last.no)
      
      sig[goal.no] <- ifelse(goal.no == last.no,0,1)
      
    }  # end of "while" on (goal.no < last.no - q) # A
    
    #  Remove putative goal from list of CP's
    sig[1] <- 0
    
    # OUTPUTS
    
    # OUTPUT OF 
    # (first, last) (row nos. of first and last times at "significant" change points) 
    # (east, north) (their coordinates) 
    # cp lists indices of waypoints identified as change points
    cp <- which(sig==1)
    
    # bsig is sig in reverse time-order
    bsig <- rep(0,nz)
    for (j in 1:nz){
      bsig[j] <- sig[nz-j+1]
    }
    
    # "times" (row nos.) and coordinates of change points
    newpp <- newp[newp > tol]
    cp.time <- newpp[sig==1]
    cp.bx1 <- bz1[sig==1]
    cp.bx2 <- bz2[sig==1]
    cp.no <- n + 1 - cp.time
    # last is cp.no in reverse order
    last <- 0*cp.no
    cp.leng <- length(cp.no)
    for (j in 1:cp.leng){
      last[j] <- cp.no[cp.leng-j+1]
    }
    # last contains row nos. of (last times at) change points 
    # (east, north) are their coordinates 
    east <- x2[last]
    north <- x1[last]
    
    # first contains row nos. of first times at change points 
    first <- 0*last
    for (j in 1:length(last)){
      first[j] <- n+2 - min(which(newp > (n+1-last[j])))
    }
    
    # PLOT WAYPOINTS AND MARK CHANGE POINTS
    png(file=filetotpng) #save plot to png
    # plot data
    # and overlay with change points in colour 
    # For plotting purposes, reverse the order of 
    # bz1, bz2 to get z1, z2
    z1 <- 0*bz1
    z2 <- 0*bz2
    nz <- length(z1)
    for (j in 1:nz){
      z1[j] <- bz1[nz-j+1]
      z2[j] <- bz2[nz-j+1]
    }
    # s is vector to index the points
    s <- c(1:nz)
    
    # NOTE: If the data have come from a GPS then 
    # it is likely that the first column of the 
    # data file contains northings and the second 
    # column contains eastings.
    # The portion of code below assumes this ordering
    # and produces a plot with the conventional 
    # orientation (north at the top)
    
    # set limits of plot
    cxlim <- c(min(bz2),max(bz2))
    cylim <- c(min(bz1) - sd(bz1diff),max(bz1))
    plot(bz2,bz1, pch=18, xlim=cxlim, ylim=cylim, xlab="X (x2)", ylab="Y (y2)")
    title(main=paste("q = ", q, ", " , "alpha = ", alpha, ", ", "N = ", N , ", ", "tol = ", tol , sep=""), 
          sub ="Blue triangle = putative goal, red star = change pt., red no. = row of data file")
    segments(bz2[s], bz1[s], bz2[s+1], bz1[s+1])
    par(new ="T")
    plot(bz2[1],bz1[1],pch=17,col="blue", xlim=cxlim, ylim=cylim, xlab="", ylab="")
    par(new ="T")
    z1bsig <- z1[bsig==1]
    z2bsig <- z2[bsig==1]
    plot(z2bsig,z1bsig,pch=8,col="red", xlim=cxlim, ylim=cylim, xlab="", ylab="")
    cpsame <- which(first == last)
    cpdiff <- which(first != last)
    lastsame <- last[cpsame]
    firstdiff <- first[cpdiff]
    lastdiff <- last[cpdiff]
    z1bsame <- z1bsig[cpsame] 
    z2bsame <- z2bsig[cpsame]
    z1bdiff <- z1bsig[cpdiff]
    z2bdiff <- z2bsig[cpdiff]
    
    if(length(cp)!=0){
      par(new ="T")
      text(z2bsame, z1bsame - 0.6*sd(bz1diff), lastsame, col="red", cex = 0.7)
      
      if(length(firstdiff)!=0){
        par(new ="T")
        text(z2bdiff, z1bdiff - 0.6*sd(bz1diff), paste(firstdiff,lastdiff,sep="-"),col="red", cex = 0.7)
      }
      
    }
    
    # par(new ="T")
    # text(z2bsame, z1bsame - 0.6*sd(bz1diff), lastsame, col="red", cex = 0.7)
    # if(length(firstdiff)!=0){
    #   par(new ="T")
    #   text(z2bdiff, z1bdiff - 0.6*sd(bz1diff), paste(firstdiff,lastdiff,sep="-"),col="red", cex = 0.7)
    # }
    
    dev.off()
    
    # "first" and "last" are the row nos. of the first and last times at "significant" change points 
    # "north" and "east" are the coordinates of these change points
    # "cps" contains "first", "last", "north" and "east". 
    cps <- cbind(first, last, north, east)
    # cps
    #Ashlee note: north = Y (X1), east = X (X2)
    # Values set by user
    # q
    # alpha
    # N
    # tol
  } else{ #the coordinate data was 2 coordinates or less - too small to run CPTauto
    cps <- matrix(-Inf,1,1)
  }

  results <- list(x1=x1,x2=x2,cps=cps) #combine important results into list to be returned for the function
  return(results)
}