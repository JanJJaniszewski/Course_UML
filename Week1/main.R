# Config
set.seed(1)

# Load packages
library(pacman)
p_load(tidyverse)
p_load(SVMMaj)

# Load and prepare data
load('Week1/Data/bank.RData')
bank <- bank[sample(nrow(bank), 1000), ]

X <- bank %>% select(-y)
y <- bank %>% pull(y)

# Execute SVMMaj function on data
svmmaj(X = X, y=y, hinge = 'absolute')
