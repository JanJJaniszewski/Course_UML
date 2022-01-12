
################################################################# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, gbm, vtreat, kernlab, caret, rpart)

######################################################### Load and prepare data 
load("~/GitHub/Course_UML/Week2/bookings.RData")
adaboost_own <- function(df, X, y, formula,iteration, maxdepth=4){
  #set initial value of the weight
  n <- nrow(df)
  weight <<- rep(1/n,n)
  
  #make an array to storage the alpha 
  alpha_list <- c()
  
  #make a list to storage the classifier fitting model
  tree_list <- list()
  
  # create a list to storage the prediction
  y_predict_list <- list()
  
  #into interations
  for (ite in 1:iteration){
    # 1 fit classifier with tree
    tree <- rpart(formula =formula,
                  data = df,
                  weights = weight,
                  method = "class",
                  control =rpart.control(maxdepth = maxdepth, cp = 0,xval = 0, 
                                         minsplit = 2, minbucket = 1))
    #append new tree to trees
    trees <- c(tree_list, tree)
    
    #collect prediction results for calculating error rate
    y_predict <- predict(tree, data = bookings_prepared, type = "class")
    y_predict_list <- c(y_predict_list, y_predict)
    
    # 2 calculate weighted error rate
    error=0
    index <- as.numeric(y != y_predict)
    err <- sum(weight * index)
    
    # 3 calculate alpha
    alpha <- log((1 - err)/err)
    alpha_list <- c(alpha_list, alpha)
                  
    # 4 set weight
    weight <- weight * exp(alpha * index)
    
    # 5 normalize the weights to sum to 1
    weight <- weight / sum(weight)
  }
  #end the loop and return output
  result <- list(alpha * y_predict)
  return(result)
}

adaboost_own(df, X, y, is_cancelled ~ ., 20)
