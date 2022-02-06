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
mDiss <- getDissFromSim(basket)

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
#Inputs:
#  mDelta: N x N matrix of dissimilarities
#  mX: N x P matrix of initial coordinates
#Outputs:
#  (normalized) Stress value of the final configuration
#############################################################

getB <- function(mDelta, mX, updating = FALSE){
  # Computes the B matrix of mX
  # updating denotes if B matrix of mY or mX is being computed
  
  mX.EuclideanDist <- getEuclideanDist(mX)
  mF <- mDelta * (1/mX.EuclideanDist)
  mF[mX.EuclideanDist == 0] <- 0
  mB <- diag(rowSums(mF)) - mF
  if(updating){
    mB[mX.EuclideanDist <= 0] <- 0
  }
  
  return(mB)
}

getNormStress <- function(mDelta, mX, mV){
  # Computes the normalized stress value of mX
  
  # Step1: eta^2_Delta
  dEta2Delta <- sum(mDelta[upper.tri(mDelta)]^2)
  
  # Step2: eta^2(X)
  dEta2 <- sum(diag(t(mX) %*% mV %*% mX))
  
  # Step3: rho(X)
  mBX <- getB(mDelta, mX)
  dRhoX <- sum(diag(t(mX) %*% mBX %*% mX))
  
  # compute raw stress
  dRawStress <- dEta2Delta + dEta2 - 2*dRhoX
  
  # compute normalised stress
  dNormStress <- dRawStress / dEta2Delta
  
  return(iNormStress)
}

getMDS <- function(mDelta, mX, eps = 1e-5){
  iN <- dim(mDelta)[1]
  mV <- iN*(diag(iN) - (1/iN)*matrix(1, nrow = iN, ncol = iN))
  dNormStress <- getNormStress(mDelta, mX, mV) 
  
  iK <- 1
  while(iK == 1 || dNormStressDiff > eps){
    iK <- iK + 1
    
    # Updating mY and mBY
    mY <- mX
    mY.EuclideanDist <- getEuclideanDist(mY)
    mBY <- getB(mDelta, mY, updating = TRUE)
    
    # Computing Moore Penrose inverse of V
    mOnes <- matrix(1, nrow = iN, ncol = iN)
    mV.MPinv <- solve(mV + (1/iN)*mOnes) - (1/iN)*mOnes   
    
    # Updating mX
    mX <- mV.MPinv %*% mBY %*% mY 
    
    # Printing stress values
    dNormStress.Prev <- dNormStress
    dNormStress <- getNormStress(mDelta, mX, mV) # stress not being updated
    dNormStressDiff <- dNormStress.Prev - dNormStress
    print(dNormStress, dNormStressDiff)
  }
  
  return(mX)
}

#############################################################
# (e) compare results from own function and package
#############################################################
# Own function results
mX0 <- matrix(1, nrow = dim(mDiss)[1], ncol = 2)
mX.ownMDS <- getMDS(mDiss, mX0)

# Implement package
# Do ratio MDS with ndim = 2
mds_ratio <- mds(diss,type="ratio")
mds_ratio$stress

#############################################################
# (f) save initial configuration
#############################################################
init <- mds(diss, type="ratio", ndim = 2, itmax = 1)
my <- my_mds(diss, X = init$conf)

#############################################################
# (g) compare plots
#############################################################
plot(mds_ratio)
plot(my)
