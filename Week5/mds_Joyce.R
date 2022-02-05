# install packages and import dataset
if (!require("pacman")) install.packages("pacman")
pacman::p_load(knitr,
               tibble,
               tidyverse,
               smacof)

#############################################################
# (b) load data and calculate dissimilarities 
#############################################################
load("~/BDS/B3 UML/w5 MDS/basket.Rdata")
set.seed(2022)
basket <- as_tibble(basket)
# Change matrix into dissimilarity
cor_basket <- cor(basket)
diss <- sim2diss(cor_basket, method = "corr")

#############################################################
# (c) write own function to get euclidean distance
#############################################################
# calculate euclidean distance measure
turn_corr_to_eucl <- function(m_corr){
  D <- sqrt(2*(1 - m_corr))
  return(D)
}
euc_d_basket <- turn_corr_to_eucl(cor_basket)

#############################################################
# (d) program the SMACOF
#Inputs:
#  Delta: n*n matrix of dissimilarities
#  X: n*p matrix
#Outputs:
#  (normalized) Stress value of the final configuration
#############################################################
# compute without weights
compute_stress <- function(Delta, X){
  
  # compute sum of squared dissimilarities eta2_delta
  eta2_delta <- sum(Delta[upper.tri(Delta)]^2)
  
  # Step1: compute squared distances
  eta2 <- nrow(X) * sum(diag(t(X) %*% X))
  
  # Step2: compute crossproduct: sum of distances times dissimilarites
  # compute Bx
  euclidean <- dist(X, method = "euclidean")
  mF <- Delta/euclidean
  # take care of division by zero issues
  mF[is.na(mF)] <- 0
  Bx <- diag(rowSums(mF)) - mF
  # slide 30: construct rho
  rho <- sum(diag(t(X) %*% Bx %*% X))
  
  # compute raw stress
  sigma_r <- eta2_delta + eta2 - 2*rho
  # compute normalised stress
  sigma_n <- sigma_r / eta2_delta
  
  return(list(stress = sigma_n, Bx = Bx))
}

my_mds <- function(diss, X = NULL, eps = 1e-6){
  
    
}

#############################################################
# (e) compare results from own function and package
#############################################################
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
