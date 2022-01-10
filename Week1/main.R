################################################################# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, SVMMaj, vtreat, kernlab)

######################################################### Load and prepare data 
load('Week1/Data/bank.RData')

# dropping weird columns
bank <- bank %>% select(-c(emp.var.rate, euribor3m))

# Transform months to numbers
bank$month <- factor(bank$month, levels = c("jan", "feb", 'mar', 'apr', 'may',
                                            'jun', 'jul', 'aug', 'sep', 'oct',
                                            'nov', 'dec'))
bank$month <- bank$month %>% as.numeric

# Transform education to numeric
bank$education <- factor(bank$education, levels = c('illiterate', 'basic.4y',
                                                    'basic.6y', 'basic.9y',
                                                    'high.school',
                                                    'professional.course',
                                                    'university.degree'))
bank$education <- bank$education %>% as.numeric

# Transform data for modeling
treatment <- designTreatmentsC(bank,
                               colnames(bank),
                               'y',
                               'yes', 
                               codeRestriction = c('clean', 'lev'))
bank_prepared <- prepare(treatment, bank)   

# set seed for sampling 
set.seed(42)

# Define train and test sets
n <- 1000 # number of obs in train
trainIndeces <- sample(nrow(bank_prepared), n)
trainSet <- bank_prepared[trainIndeces,]
testSet <- bank_prepared[-trainIndeces, ]

# Scale X and add a column with 1's to X
X.train <- trainSet %>% select(-y) %>% scale
X.train <- cbind(1, X.train)

X.test <- testSet %>% select(-y) %>% scale
X.test <- cbind(1, X.test)

# Set y to -1 and 1 for negative and positive cases
y.train <- trainSet %>% pull(y) %>% as.numeric
y.train[y.train == 1] <- -1
y.train[y.train == 2] <- 1

y.test <- testSet %>% pull(y) %>% as.numeric
y.test[y.test == 1] <- -1
y.test[y.test == 2] <- 1

##################################################### Defining own SVM function

# Initializing loss, update and majorizing parameter functions based on hinge
quadSVMLoss <- function(y, q, lambda, w){
  return( sum(pmax(0, 1 - y * q)^2) + lambda * crossprod(w) )
}

quadSVMParams <- function(y, q){
  a <- 1
  b <- ifelse(y == -1, pmin(q, -1), pmax(q, 1))
  params <- list(a,b)
  names(params) <- c("a","b")
  return(params)
}
  
quadSVMUpdate <- function(X, A, lambda, P, b, Z){
  return( Z %*% t(X) %*% b )
}

absSVMLoss <- function(y, q, lambda, w){
  return( sum(pmax(0, 1 - y * q)) + lambda * crossprod(w) )
}

absSVMParams <- function(y, q){
  a <- pmax(4*abs(1 - y * q), 1e-4)^(-1)
  b <- y*(a + 1/4)
  params <- list(a,b)
  names(params) <- c("a","b")
  return(params)
}

absSVMUpdate <- function(X, A, lambda, P, b, Z){
  return( solve(t(X) %*% A %*% X + lambda * P, t(X) %*% b) )
}

ownSVMMAJ <- function(X, y, lambda, eps, hinge, debug = FALSE){
  ####
  # Purpose :
  #   Fit an SVM with absolute or quadratic error for predicting y based on X.
  # Inputs  : 
  #   X           smth
  #   y           vector of 1's and (-1)'s, treatment variable
  #   lambda      smth
  #   eps         smth
  #   hinge       smth
  #   debug       smth
  # Returns :
  #   constant    smth
  #   w           smth
  ####
  
  # Initializing values
  m <- ncol(X) - 1 # Number of columns (excluding the constant column)
  w <- matrix(1, m, 1) # Initial weights
  constant <- 1 # Initial c
  v <- t(cbind(constant, t(w))) # [c, wT]
  v.prevprev <- v
  v.prev <- v
  P <- diag(1, m+1)
  P[1,1] <- 0
  q <- X %*% v
  Z <- solve(t(X) %*% X + lambda * P)
  
  # Defining appropriate functions based on hinge
  if(hinge == 'absolute'){
    SVMLoss <- absSVMLoss
    SVMUpdate <- absSVMUpdate
    SVMParams <- absSVMParams
  } else { # quadratic hinge
    SVMLoss <- quadSVMLoss
    SVMUpdate <- quadSVMUpdate
    SVMParams <- quadSVMParams
  }
  
  # Entering while function
  k <- 0
  l_svm <- SVMLoss(y, q, lambda, w)
  
  while( k == 0 || ((l_svm_old - l_svm)/l_svm_old) > eps ){
    # BUG: This terminates because l_svm increases after some point which 
    # shouldn't happen.
    
    # Update number of iterations
    k = k+1
    
    # Assign previous run loss to old loss
    l_svm_old <- l_svm
    
    # Predict q
    q <- X %*% v
    
    # Compute a, b, c, A
    params <- SVMParams(y, q)
    a <- params$a
    b <- params$b
    A <- diag(x = a %>% as.vector, n, n)
    
    # Update v based on hinge
    v.prevprev <- v.prev
    v.prev <- v
    v <- SVMUpdate(X, A, lambda, P, b, Z)
    w <- v[2:length(v)]
    
    # compute the new loss
    l_svm <- SVMLoss(y, q, lambda, w)
    
    # Print update information for debugging
    if(debug){
      print('---------------')
      print(k) # number of iterations
      print(l_svm[1,1])
      print( ((l_svm_old - l_svm)/l_svm_old)[1,1]) # loss function improvement, must be positive
    }
  }
  
  output <- list(w, constant)
  names(output) <- c("w","c")
  return(v)
}

########################################################## Defining Diagnostics

F1Score <- function(y, yhat){
  TP <- sum(yhat[y == 1] == 1)
  P <- sum(y[y == 1])
  Pstar <- sum(yhat[yhat == 1])
  recall <- TP/P
  precision <- TP/Pstar
  f1 <- 2*precision*recall/(precision + recall)
  return(f1)
}
  
hitRate <- function(y, yhat){
  TP <- sum(yhat[y == 1] == 1)
  TN <- sum(yhat[y == -1] == -1)
  n <- length(y)
  hitRate <- (TP + TN)/n
  return(hitRate)
}

############################################################# Computing results

# Initializing params
eps <- 1e-3 # stopping criterion for improvement of function
lambda <- 15 # lambda parameter

# Initializing results table 1
resTable <- matrix(nrow = 4, ncol = 2)
colnames(resTable) <- c("Hit rate", "F1 score")
rownames(resTable) <- c("Own, abs hinge", "Package, abs hinge", 
                        "Own, quad hinge", "Package, quad hinge")

# Fitting abs hinge
resOwnSVMMAJabs <- ownSVMMAJ(X.train, y.train, lambda = lambda, eps,
                          hinge = 'absolute')
resTable[1,1] <- hitRate(y.test, ifelse(X.test%*%resOwnSVMMAJabs > 0, 1, -1)) 
resTable[1,2] <- F1Score(y.test, ifelse(X.test%*%resOwnSVMMAJabs > 0, 1, -1))

resSVMMAJabs <- svmmaj(X.train, y.train, lambda = lambda,
                              hinge = 'absolute')
resTable[2,1] <- hitRate(y.test, ifelse(X.test%*%resSVMMAJabs$beta > 0, 1, -1)) 
resTable[2,2] <- F1Score(y.test, ifelse(X.test%*%resSVMMAJabs$beta > 0, 1, -1))

# Fitting quad hinge
resOwnSVMMAJquad <- ownSVMMAJ(X.train, y.train, lambda = lambda, eps,
                           hinge = 'quadratic')
resTable[3,1] <- hitRate(y.test, ifelse(X.test%*%resOwnSVMMAJquad > 0, 1, -1)) 
resTable[3,2] <- F1Score(y.test, ifelse(X.test%*%resOwnSVMMAJquad > 0, 1, -1))

resSVMMAJquad <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'quadratic')
resTable[4,1] <- hitRate(y.test, ifelse(X.test%*%resSVMMAJquad$beta > 0, 1, -1)) 
resTable[4,2] <- F1Score(y.test, ifelse(X.test%*%resSVMMAJquad$beta > 0, 1, -1))

# Printing table with results
resTable

# Initializing results table 2
resTable2 <- matrix(nrow = 4, ncol = 1)
colnames(resTable2) <- c("Loss")
rownames(resTable2) <- c("RBF abs hinge", "RBF quad hinge", 
                        "Poly abs hinge", "Poly quad hinge")

# Fitting nonlinear SVM's
resSVMRBFabs <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'absolute', 
                              kernel = rbfdot, kernel.sigma = 1)
resTable2[1] <- resSVMRBFabs$loss 

resSVMRBFquad <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'quadratic', 
                              kernel = rbfdot, kernel.sigma = 1)
resTable2[2] <- resSVMRBFquad$loss

resSVMPolyabs <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'absolute', 
                       kernel = polydot, kernel.scale = 1, kernel.degree = 1, 
                       kernel.offset = 1)
resTable2[3] <- resSVMPolyabs$loss

resSVMPolyquad <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'quadratic', 
                        kernel = polydot, kernel.scale = 1, kernel.degree = 1, 
                        kernel.offset = 1)
resTable2[4] <- resSVMPolyquad$loss

# Printing table with results
resTable2


