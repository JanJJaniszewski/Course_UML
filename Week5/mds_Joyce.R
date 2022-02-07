# install packages and import dataset
if (!require("pacman")) install.packages("pacman")
pacman::p_load(knitr,
               tibble,
               tidyverse,
               smacof)

#############################################################
# (b) load data and calculate dissimilarities 
#############################################################
#load("~/BDS/B3 UML/w5 MDS/basket.Rdata")
load("~/Desktop/Unsupervised Machine Learning/Course_UML/Week5/basket.Rdata")
set.seed(2022)

getDissFromSim <- function(mSim){ # right method of slide 16
  ### Add DocString ###
  
  vDimSim <- dim(mSim)
  mDiss <- matrix(0, nrow = vDimSim[1], ncol = vDimSim[2])
  for(i in 1:vDimSim[1]){
    for(j in i:vDimSim[2]){
      mDiss[i,j] <- as.numeric(log((mSim[i,i]*mSim[j,j])/(mSim[i,j]*mSim[j,i])))
      mDiss[j,i] <- mDiss[i,j]
    }  
  }
  return(mDiss)
}

#diss <- sim2diss(cor_basket, method = "corr")
mDissimilarities <- getDissFromSim(basket)

#############################################################
# (c) write own function to get euclidean distance
#############################################################

getEuclideanDist <- function(mData){
  # Calculates the Euclidean distances between the rows of an N x P matrix mData
  # mEuclideanDist has dimension N x N. 
  
  iN <- dim(mData)[1]
  mEuclideanDist <- matrix(0, nrow = iN, ncol = iN)
  for(i in 1:iN){
    for(j in i:iN){
      squaredDiff <- (mData[i,] - mData[j,])^2
      mEuclideanDist[i,j] <- sqrt(sum(squaredDiff))
      mEuclideanDist[j,i] <- mEuclideanDist[i,j]
    }
  }
  return(mEuclideanDist)
}

#euc_d_basket <- turn_corr_to_eucl(cor_basket)

#############################################################
# (d) program the SMACOF
#############################################################

getBX <- function(mDelta, mX){
  # Computes the mBX matrix
  #Inputs:
  #  mDelta: N x N matrix of dissimilarities
  #  mX: N x P matrix of initial coordinates
  #Outputs:
  #  mBX
  
  mX.EuclideanDist <- getEuclideanDist(mX)
  mF <- mDelta * (1/mX.EuclideanDist)
  mF[mX.EuclideanDist == 0] <- 0 # This values are NA since denominator is 0
  mBX <- diag(rowSums(mF)) - mF
  return(mBX)
}

getBY <- function(mDelta, mY){
  # Computes the mBY matrix
  #Inputs:
  #  mDelta: N x N matrix of dissimilarities
  #  mX: N x P matrix of initial coordinates
  #Outputs:
  #  mBY
  
  mBY <- matrix(0, nrow = iN, ncol = iN)
  mY.EuclideanDist <- getEuclideanDist(mY)
  iN <- dim(mY)[1]
  for(i in 1:(iN-1)){
    for(j in (i+1):iN){
      if(mY.EuclideanDist[i,j] > 0){
        dB <- mDelta[i,j]/mY.EuclideanDist[i,j]
        vE_i = rep(0,iN)
        vE_j = rep(0,iN)
        vE_i[i] = 1
        vE_j[j] = 1
        mA <- (vE_i - vE_j) %*% t(vE_i - vE_j)
        mBY <- mBY + dB*mA
      }
    }
  }
  return(mBY)
}

getStress <- function(mDelta, mX, mV, dEta2Delta){
  # Computes the raw stress value of mX
  #Inputs:
  #  mDelta: N x N matrix of dissimilarities
  #  mX: N x P matrix of initial coordinates
  #Outputs:
  #  dRawStress
  
  # Compute eta^2(X)
  dEta2 <- sum(diag( t(mX) %*% mV %*% mX ))
  
  # Compute rho(X)
  mBX <- getBX(mDelta, mX)
  dRhoX <- sum(diag(t(mX) %*% mBX %*% mX))
  
  # Compute stress
  dRawStress <- dEta2Delta + dEta2 - 2*dRhoX

  return(dRawStress)
}

getMDS <- function(mDelta, mX, eps = 10e-6, debug = FALSE){
  # Returns the MDS approximation of the dissimilarity matrix mDelta, starting
  #  at the provided mX. 
  #Inputs:
  #  mDelta: N x N matrix of dissimilarities
  #  mX: N x P matrix of initial coordinates
  #Outputs:
  #  dNormStress 
  #  mX: 
  
  iN <- dim(mDelta)[1]
  mV <- iN*(diag(iN) - (1/iN)*matrix(1, nrow = iN, ncol = iN))
  mOnes <- matrix(1, nrow = iN, ncol = iN)
  mV.MPinv <- solve(mV + (1/iN)*mOnes) - (1/iN)*mOnes
  
  # Initializing Normalized Stress
  dEta2Delta <- sum(mDelta[upper.tri(mDelta)]^2)
  dNormStress <- getStress(mDelta, mX, mV, dEta2Delta)/dEta2Delta 
  
  iK <- 1
  while(iK == 1 || dNormStressDiff > eps){
    iK <- iK + 1
    
    # Updating mY and mBY
    mY <- mX
    mY.EuclideanDist <- getEuclideanDist(mY)
    mBY <- getBY(mDelta, mY)
    
    # Updating mX
    mX = mV.MPinv %*% mBY %*% mY # Moore Penrose approach
    
    # Printing stress values
    dPrevNormStress <- dNormStress
    dNormStress <- getStress(mDelta, mX, mV, dEta2Delta)/dEta2Delta
    dNormStressDiff <- dPrevNormStress - dNormStress # should always be >= 0
    if(debug){
      if(iK == 2){
        print("Running Own MDS:")
      }
      print(sprintf("Iteration: %i, Normalized stress: %f, Difference: %f",
                    iK, dNormStress, dNormStressDiff))
    }
  }
  
  output <- list(dNormStress, mX)
  names(output) <- c("dNormStress", "mX")
  return(output)
}

# Own function results
iN <- dim(mDissimilarities)[1]
iP <- 2
mX0 <- matrix(seq(1:(iN*iP)), nrow = iN, ncol = iP)
mX.ownMDS <- getMDS(mDissimilarities, mX0, debug = TRUE)

# Configuration plot
vDim1.ownMDS <- mX.ownMDS$mX[,1]
vDim2.ownMDS <- mX.ownMDS$mX[,2]
plot(vDim1.ownMDS, vDim2.ownMDS, main = "Configuration plot: Own SMACOF results", 
     xlab = "Dim 1", ylab = "Dim 2", pch = 18, col = "blue", xlim = c(-3,3), 
     ylim = c(-3,3))
text(vDim1.ownMDS, vDim2.ownMDS, colnames(basket), cex=0.6, pos=3, col="red")

#############################################################
# (e) compare own function with package function
#############################################################



#############################################################
# (f) save initial configuration
#############################################################

# Implement package
# Do ratio MDS with ndim = 2
mX.packageMDS <- mds(mDissimilarities, type="ratio", ndim = 2)

# Own MDS results with package initialization 
init <- mds(mDissimilarities, type="ratio", ndim = 2, itmax = 1)
mX.ownMDS.packageInit <- getMDS(mDelta = mDissimilarities, mX = init$conf, debug = TRUE)

#############################################################
# (g) compare plots
#############################################################

# Configuration plot, own MDS with package init
vDim1.ownMDS.packageInit <- mX.ownMDS.packageInit$mX[,1]
vDim2.ownMDS.packageInit <- mX.ownMDS.packageInit$mX[,2]
plot(vDim1.ownMDS.packageInit, vDim2.ownMDS.packageInit, 
     main = "Configuration plot: Own ratio MDS results", 
     xlab = "Dim 1", ylab = "Dim 2", pch = 18, col = "blue", xlim = c(-4,4), 
     ylim = c(-4,4))
text(vDim1.ownMDS.packageInit, vDim2.ownMDS.packageInit, colnames(basket),
     cex=0.6, pos=3, col="red")

# Shepard plot, own MDS with package init
vOwnConfigDistances <- as.vector(getEuclideanDist(mX.ownMDS.packageInit$mX))
vDissimilarities <- as.vector(mDissimilarities)
plot(vDissimilarities, vOwnConfigDistances, 
     main = "Shepard Plot: Own ratio MDS results", 
     xlab = "Dissimilarities", ylab = "Configuration Distances", pch = 18, 
     col = "blue", xlim = c(0,7), ylim = c(0,6))

# Configuration plot, package MDS
vDim1.packageMDS <- mX.packageMDS$conf[,1]
vDim2.packageMDS <- mX.packageMDS$conf[,2]
plot(vDim1.packageMDS, vDim2.packageMDS, 
     main = "Configuration plot: Package ratio MDS results", 
     xlab = "Dim 1", ylab = "Dim 2", pch = 18, col = "blue", xlim = c(-1,1), 
     ylim = c(-1,1))
text(vDim1.packageMDS, vDim2.packageMDS, colnames(basket),
     cex=0.6, pos=3, col="red")

# Shepard plot, package MDS 
vPackageConfigDistances <- as.vector(getEuclideanDist(mX.packageMDS$conf))
plot(vDissimilarities, vPackageConfigDistances,
     main = "Shepard Plot: Package ratio MDS results", 
     xlab = "Dissimilarities", ylab = "Configuration Distances", pch = 18, 
     col = "blue", xlim = c(0,6), ylim = c(0,3))

#############################################################
# (h) MDS with ordinal transformation
#############################################################

mX.packageOrdinalMDS <- mds(mDissimilarities, type="ordinal", ndim = 2)

# Configuration plot, package ordinal MDS
vDim1.packageOrdinalMDS <- mX.packageOrdinalMDS$conf[,1]
vDim2.packageOrdinalMDS <- mX.packageOrdinalMDS$conf[,2]
plot(vDim1.packageOrdinalMDS, vDim2.packageOrdinalMDS, 
     main = "Configuration plot: Package ordinal MDS results", 
     xlab = "Dim 1", ylab = "Dim 2", pch = 18, col = "blue", xlim = c(-1,1), 
     ylim = c(-1,1))
text(vDim1.packageOrdinalMDS, vDim2.packageOrdinalMDS, colnames(basket),
     cex=0.6, pos=3, col="red")

# Kruskal's rule of thumb (only applicable to ordinal MDS) and elbow rule
vStressVals <- rep(0,iN-1)
for(i in 1:length(vStressVals)){
  vStressVals[i] <- mds(mDissimilarities, type="ordinal", ndim = i)$stress
}

plot(seq(1,iN-1), vStressVals, xlab = "Dimensionality", ylab = "Stress", 
     pch = 18, col = "blue", xlim = c(1,iN-1), ylim = c(0,0.25))

# According to Kruskal's rule of thumb 2 dim is fair, 3 is approx good and 4 is 
# already excellent. By the elbow rule 3 dim seems optimal. 

# Shepard plot, own MDS with package init
vPackageOrdConfigDistances <- as.vector(getEuclideanDist(mX.packageOrdinalMDS$conf))
plot(vDissimilarities, vPackageOrdConfigDistances,
     main = "Shepard Plot: Package ordinal MDS results", 
     xlab = "Dissimilarities", ylab = "Configuration Distances", pch = 18, 
     col = "blue", xlim = c(0,3), ylim = c(0,6))



