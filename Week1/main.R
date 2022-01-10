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
      print(l_svm_old - l_svm) # loss function improvement, must be positive
    }
    
    # Trying to fix the overshoot 
    if(l_svm_old - l_svm < 0){
      v <- (1/2)*(v.prevprev + v.prev)
      w <- v[2:length(v)]
      l_svm <- SVMLoss(y, q, lambda, w)
      print(l_svm_old - l_svm)
    }
  }
  
  output <- list(w, constant)
  names(output) <- c("w","c")
  return(v)
}

########################################################## Defining Diagnostics

F1Score <- function(y, yhat){
  TP <- sum(y[y == 1] == yhat[yhat == 1])
  P <- sum(y[y == 1])
  Pstar <- sum(yhat[yhat == 1])
  recall <- TP/P
  precision <- TP/Pstar
  f1 <- 2*precision*recall/(precision + recall)
  return(f1)
}
  
hitRate <- function(y, yhat){
  TP <- sum(y[y == 1] == yhat[yhat == 1])
  TN <- sum(y[y == -1] == yhat[yhat == -1])
  n <- length(y)
  hitRate <- (TP + TN)/n
  return(hitRate)
}

############################################################# Computing results

# Initializing params
eps <- 1e-25 # stopping criterion for improvement of function
lambda <- 15 # lambda parameter

X <- X.train
y <- y.train
hinge <- 'absolute'
debug <- TRUE

# Fitting abs hinge
resSVMMAJabs <- ownSVMMAJ(X.train, y.train, lambda = lambda, eps,
                          hinge = 'absolute', debug = TRUE)
# These dont work because we need to define a way of predicting 1 or -1 
ownSVMAbsHitRate <- hitRate(ytest, Xtest%*%resSVMMAJabs) 
ownSVMAbsF1Score <- F1Score(ytest, Xtest%*%resSVMMAJabs)

resPackageSVMMAJabs <- svmmaj(X.train, y.train, lambda = lambda,
                              hinge = 'absolute')

# Fitting quad hinge
resSVMMAJquad <- ownSVMMAJ(X.train, y.train, lambda = lambda, eps,
                           hinge = 'quadratic', debug = TRUE)

resPackageSVMMAJquad <- svmmaj(X.train, y.train, lambda = lambda,
                               hinge = 'quadratic')


# Fitting nonlinear SVM's
resSVMRBFabs <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'absolute', 
                              kernel = rbfdot, kernel.sigma = 1)

resSVMRBFquad <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'quadratic', 
                              kernel = rbfdot, kernel.sigma = 1)

resSVMRBFabs <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'absolute', 
                       kernel = polydot, kernel.scale = 1, kernel.degree = 1, 
                       kernel.offset = 1)

resSVMRBFquad <- svmmaj(X.train, y.train, lambda = lambda, hinge = 'quadratic', 
                        kernel = polydot, kernel.scale = 1, kernel.degree = 1, 
                        kernel.offset = 1)




