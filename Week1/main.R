# Config
set.seed(1)

# Load packages
library(pacman)
p_load(tidyverse)
p_load(SVMMaj)
p_load(vtreat)

# Load and prepare data
load('Week1/Data/bank.RData')
n <- 1000
bank <- bank[sample(nrow(bank), n), ]

bank <- bank %>% select(-c(emp.var.rate, euribor3m))

bank$month <- factor(bank$month, levels = c("jan", "feb", 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'))
bank$month <- bank$month %>% as.numeric

bank$education <- factor(bank$education, levels = c('illiterate', 'basic.4y', 'baic.6y', 'basic.9y', 'high.school', 'professional.course', 'university.degree'))
bank$education <- bank$education %>% as.numeric



treatment <- designTreatmentsC(bank,
                               colnames(bank),
                               'y',
                               'yes', 
                               codeRestriction = c('clean', 'lev'))
                               
bank_prepared <-  prepare(treatment, bank)   

X <- bank_prepared %>% select(-y) %>% scale
X <- cbind(1, X)
y <- bank_prepared %>% pull(y) %>% as.numeric
y[y == 1] <- -1
y[y == 2] <- 1

# Execute SVMMaj function on data
svmmaj(X = X, y=y, hinge = 'absolute')

# Own SVM function
error <- 0.1
m <- ncol(X) - 1
w <- matrix(0.1, m, 1)
lambda <- 1
constant <- 0
v_T <- cbind(constant, t(w))
P <- diag(1, m+1)
P[1,1] <- 0

# Loss
l_svm <- 1
l_svm_old <- -Inf

# Majorization for later
Z <- (crossprod(X) + lambda * P)^(-1) %*% t(X)

while(((l_svm - l_svm_old) / l_svm) > error){
  l_svm_old <- l_svm
  
  # Predictions
  q <- X %*% t(v_T)
  
  # Absolute hinge computation
  a <- (4 * abs(1 - y * q)) ^ (-1)
  b <- y*(a + 1/4)
  c <- a + 1/2 + abs(1 - y * q)/4
  A <- diag(x=a %>% as.vector, n, n)
  
  # Loss computation
  l_svm <- sum(a * q^2) - 2 * sum(b * q) + sum(c) + lambda * crossprod(w)
  
  # Updating v
  v <- Z %*% b
  constant <- v[1]
  w <- v[2:length(v)]
  
  print(l_svm)
}