
################################################################# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, gbm, vtreat, kernlab, caret)

######################################################### Load and prepare data 
load("~/GitHub/Course_UML/Week2/bookings.RData")
adaboost_own <- function(data,){
  formula <- 
}
