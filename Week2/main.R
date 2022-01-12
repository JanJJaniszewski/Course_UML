#### Load packages and init ####################################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, vtreat, caret)
set.seed(42)


#### Data #################################################################
##### Load #####
load('Week2/bookings.RData')

##### Analyze #####
bookings %>% group_by(is_cancelled) %>% summarise_all(mean) %>% head()
dim(bookings)


##### Prepare #####
# Define train and test sets
n <- (0.75 * nrow(bookings)) %>% round()
trainIndeces <- sample(nrow(bookings), n)
trainSet <- bookings[trainIndeces,]
testSet <- bookings[-trainIndeces, ]

treatment <- designTreatmentsC(trainSet,
                               colnames(trainSet),
                               'is_cancelled',
                               'yes', 
                               codeRestriction = c('clean', 'lev'))
trainSet <- prepare(treatment, trainSet)  
testSet <- prepare(treatment, testSet)

# Scale X and add a column with 1's to X
X.train <- trainSet %>% select(-is_cancelled) %>% scale
X.test <- testSet %>% select(-is_cancelled) %>% scale

# Set y to -1 and 1 for negative and positive cases
y.train <- trainSet %>% pull(is_cancelled) %>% as.numeric %>% `-`(1)
y.test <- testSet %>% pull(is_cancelled) %>% as.numeric %>% `-`(1)

