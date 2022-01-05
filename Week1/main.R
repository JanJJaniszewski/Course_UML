# Config
set.seed(1)

# Load packages
library(pacman)
p_load(tidyverse)
p_load(SVMMaj)
p_load(vtreat)

# Load and prepare data
load('Week1/Data/bank.RData')
bank <- bank[sample(nrow(bank), 1000), ]

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
m <- ncol(X)
w <- matrix(1.1, m, 1)
lambda <- 1
c <- 1
v_T <- cbind(c, t(w))

q <- X %*% t(v_T)

l_svm <- sum(max(0, 1 - y*q)) + (lambda * crossprod(w))

