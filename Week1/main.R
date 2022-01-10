################################################################# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, SVMMaj, vtreat, kernlab)

######################################################### Load and prepare data 
load('Week1/Data/bank.RData')

# Take only a sample of n
n <- 1000 
bank <- bank[sample(nrow(bank), n), ]

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

# Scale X and add a column with 1's to X
X <- bank_prepared %>% select(-y) %>% scale
X <- cbind(1, X) 

# Set y to -1 and 1 for negative and positive cases
y <- bank_prepared %>% pull(y) %>% as.numeric
y[y == 1] <- -1
y[y == 2] <- 1

##################################################### Defining own SVM function

# Initializing loss, update and majorizing parameter functions based on hinge
quadSVMLoss <- function(y, q, lambda, w){
  return( sum(pmax(0, 1 - y * q)^2) + lambda * crossprod(w) )
}

quadSVMParams <- function(y, q){
  a <- 1
  b <- ifelse(y == -1, min(q, -1), max(q, 1))
  c <- ifelse(y == -1, 
              ifelse(q <= -1, a - 2*(1 + q) + (1 + q)^2, 1),
              ifelse(q > 1, a - 2*(1 - q) + (1 - q)^2, 1))
  params <- list(a,b,c)
  names(params) <- c("a","b","c")
  return(params)
}
  
quadSVMUpdate <- function(X, A, lambda, P, b){
  return( solve(t(X) %*% X + lambda * P) %*% t(X) %*% b )
}

absSVMLoss <- function(y, q, lambda, w){
  return( sum(pmax(0, 1 - y * q)) + lambda * crossprod(w) )
}

absSVMParams <- function(y, q){
  a <- pmax((4 * abs(1 - y * q)), 1e-4) ^ (-1)
  b <- y*(a + 1/4)
  c <- a + 1/2 + abs(1 - y * q)/4
  params <- list(a,b,c)
  names(params) <- c("a","b","c")
  return(params)
}

absSVMUpdate <- function(X, A, lambda, P, b){
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
  w <- matrix(0.1, m, 1) # Initial weights
  constant <- 0 # Initial c
  v <- t(cbind(constant, t(w))) # [c, wT]
  P <- diag(1, m+1)
  P[1,1] <- 0
  q <- X %*% v
  
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
  l_svm_old <- l_svm + l_svm * eps + 1
  
  while(((l_svm_old - l_svm) / l_svm_old) > eps){
    # Update number of iterations
    k = k+1
    
    # Assign previous run loss to old loss
    l_svm_old <- l_svm
    
    # Predict q
    q <- X %*% v
    
    # Compute a, b, c, A
    params <- SVMParams(y,q)
    a <- params$a
    b <- params$b
    c <- params$c
    A <- diag(x = a %>% as.vector, n, n)
    
    # compute the new loss
    l_svm <- SVMLoss(y, q, lambda, w)
    
    # Update v based on hinge
    v <- SVMUpdate(X, A, lambda, P, b)
    constant <- v[1]
    w <- v[2:length(v)]
    
    # Print update information for debugging
    if(debug){
      print('---------------')
      print(k) # number of iterations
      print(l_svm[1,1]) # loss function value
    }
  }
  
  output <- list(w, constant)
  names(output) <- c("w","c")
  return(v)
}

############################################################# Computing results

# Initializing params
eps <- 1e-15 # Stopping criterion for improvement of function
lambda <- 15 # Lambda (parameter)

# Fitting abs hinge
set.seed(42)
resSVMMAJabs <- ownSVMMAJ(X, y, lambda = lambda, eps, hinge = 'absolute')
set.seed(42)
resPackageSVMMAJabs <- svmmaj(X, y, lambda = lambda, hinge = 'absolute')
resPackageSVMMAJabs$beta - resSVMMAJabs

# Fitting quad hinge
set.seed(42)
resSVMMAJquad <- ownSVMMAJ(X, y, lambda = lambda, eps, hinge = 'quadratic')
set.seed(42)
resPackageSVMMAJquad <- svmmaj(X, y, lambda = lambda, hinge = 'quadratic')


# Fitting nonlinear SVM's
resSVMRBFabs <- svmmaj(X, y, lambda = lambda, hinge = 'absolute', 
                              kernel = rbfdot, kernel.sigma = 1)

resSVMRBFquad <- svmmaj(X, y, lambda = lambda, hinge = 'quadratic', 
                              kernel = rbfdot, kernel.sigma = 1)

resSVMRBFabs <- svmmaj(X, y, lambda = lambda, hinge = 'absolute', 
                       kernel = polydot, kernel.scale = 1, kernel.degree = 1, 
                       kernel.offset = 1)

resSVMRBFquad <- svmmaj(X, y, lambda = lambda, hinge = 'quadratic', 
                        kernel = polydot, kernel.scale = 1, kernel.degree = 1, 
                        kernel.offset = 1)

# Comparing results


