#### Config ####
set.seed(1)

#### Load packages ####
library(pacman)
p_load(tidyverse)
p_load(SVMMaj)
p_load(vtreat)

#### Load and prepare data ####
load('Week1/Data/bank.RData')

# Take only a sample of n
n <- 1000 
bank <- bank[sample(nrow(bank), n), ]

# dropping weird columns
bank <- bank %>% select(-c(emp.var.rate, euribor3m))

# Transform months to numbers
bank$month <- factor(bank$month, levels = c("jan", "feb", 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'))
bank$month <- bank$month %>% as.numeric

# Transform education to numberic
bank$education <- factor(bank$education, levels = c('illiterate', 'basic.4y', 'baic.6y', 'basic.9y', 'high.school', 'professional.course', 'university.degree'))
bank$education <- bank$education %>% as.numeric

# Transform data for modeling
treatment <- designTreatmentsC(bank,
                               colnames(bank),
                               'y',
                               'yes', 
                               codeRestriction = c('clean', 'lev'))
bank_prepared <-  prepare(treatment, bank)   

# Scale X and add a column with 1's to X
X <- bank_prepared %>% select(-y) %>% scale
X <- cbind(1, X) 

# Set y to -1 and 1 for negative and positive cases
y <- bank_prepared %>% pull(y) %>% as.numeric
y[y == 1] <- -1
y[y == 2] <- 1

# Execute SVMMaj function on data
# svmmaj(X = X, y=y, hinge = 'absolute')

#### Own SVM function ####
# Initialize
error <- 0.0001 # Stopping criterion for improvement of function
m <- ncol(X) - 1 # Number of columns (excluding the constant column)
w <- matrix(0.1, m, 1) # Initial weights
lambda <- 4 # Lambda (parameter)
constant <- 0 # Initial c
v <- t(cbind(constant, t(w))) # [c, wT]

# P matrix
P <- diag(1, m+1)
P[1,1] <- 0


# Entering while function
k <- 1
l_svm_old <- 0
l_svm <- 0

while((k <= 2) | ((l_svm_old - l_svm) / l_svm_old) > error){
  k = k+1
  # Assign previous run loss to old loss
  l_svm_old <- l_svm
  
  # Predict q
  q <- X %*% v
  
  q %>% head
  y %>% head
  
  # Compute a, b, c, A using absolute hinge
  a <- (4 * abs(1 - y * q)) ^ (-1)
  b <- y*(a + 1/4)
  c <- a + 1/2 + abs(1 - y * q)/4
  A <- diag(x=a %>% as.vector, n, n)
  
  # compute the new loss
  l_svm <- sum(pmax(0, 1 - y * q)) + lambda * crossprod(w)
  
  # Update v
  v <- solve(t(X) %*% A %*% X + lambda * P, t(X) %*% b)
  
  constant <- v[1]
  w <- v[2:length(v)]
  
  # Print update information
  print('---------------')
  print(k)
  print(l_svm[1,1])
}
