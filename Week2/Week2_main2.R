#### Load packages and init ################################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, vtreat, caret, rpart, ohenery, gbm, pROC)
set.seed(42)

##### Own AdaBoost #########################################################

fitAdaBoost <- function(data.train, response.var, n.trees = 100, max.depth = 5,
                        debug = 0){
  ####################################################################
  # Purpose :
  #   Fit AdaBoost classification algorithm on response.var using data in 
  #     data.train, using (rpart) trees as base learner. 
  # Inputs  : 
  #   data.train    [...]
  #   response.var  [...]
  #   n.trees       [...]
  #   max.depth     [...]
  #   debug         [...]
  # Returns :
  #   output        [...]
  ####################################################################
  
  # Initializing formula, response, predictors and outputs
  formula <- as.formula(paste(response.var, "~ .")) # defining formula for rpart
  response.index <- which(colnames(data.train) == response.var)
  X <- data.train[,-response.index]
  y <- data.train[,response.index]
  N <- dim(data.train)[1] # number of observations
  w <- rep(1/N, N) # initializing observation weights
  G <- list() # initializing list of trees to be fit
  alpha <- list() # initializing list of tree weights
  
  # Iterating over all trees
  for (m in 1:n.trees){
    G[[m]] <- rpart(formula, data.train, weights = w, 
                    control = rpart.control(maxdepth = max.depth, cp = 0, 
                                            xval = 0, minsplit = 2,
                                            minbucket = 1)) # fitting m-th tree
    preds_Gm <- predict(G[[m]], X) # computing in-sample predictions
    preds_Gm <- sign(preds_Gm) # converting to 0 or 1 predictions
    err_m <- sum(w * (y != preds_Gm)) / sum(w)
    alpha[[m]] <- log((1 - err_m) / err_m)
    w <- (w * exp(alpha[[m]] * (y != preds_Gm))) %>% normalize
    if(debug){
      print('------------')
      print(m)
    }
    if(debug >= 2){
        print('SD(w)')
        print(sd(w))
        print('Alpha')
        print(alpha[[m]])
        print('Error')
        print(err_m)
    }
  }
  
  # Initializing output
  output <- list()
  output[['alpha']] <- alpha
  output[['G']] <- G
  return(output)
}

predictAdaBoost <- function(data.test, response.var, fitAdaBoost.res){
  ####################################################################
  # Purpose :
  #   [...]
  # Inputs  : 
  #   data.test         [...]
  #   response.var      [...]
  #   fitAdaBoost.res   [...]
  # Returns :
  #   preds_final       [...]
  ####################################################################
  
  M <- fitAdaBoost.res$alpha %>% length
  alpha <- fitAdaBoost.res$alpha
  G <- fitAdaBoost.res$G
  response.index <- which(colnames(data.test) == response.var)
  X <- data.test[,-response.index]
  y <- data.test[,response.index]
  preds <- data.frame(matrix(nrow = nrow(X), ncol = M))
  
  for(m in 1:M){
    preds[,m] <- alpha[[m]] * sign(predict(G[[m]], X)) 
  }
  preds_final <- sign(apply(preds, MARGIN = 1, FUN = sum))
  return(preds_final)
}

#### Data #################################################################

##### Loading
load('Week2/bookings.RData')

##### Defining response variable
response.var <- "is_cancelled"
response.index <- which(colnames(bookings) == response.var)

##### Processing
# Creating dummies for categorical variables to avoid different test and train sample structures
bookingsData <- data.frame(cbind(bookings[,response.index],
                                 model.matrix(is_cancelled ~ . + 0,
                                              data = bookings))) 

# Defining train and test sets
numObsTrainSet <- (0.75 * nrow(bookings)) %>% round() # using 3/4 of obs to train
trainIndeces <- sample(nrow(bookings), numObsTrainSet) # extracting train data indeces

# Defining train and test set for own implementation
trainSet <- bookingsData[trainIndeces,] %>% data.frame
trainSet[,response.index] <- trainSet[,response.index] %>% as.numeric %>%
  replace(. == 1, -1) %>% replace(. == 2, 1) # Set y to 0 and 1 for negative and positive cases

testSet <- bookingsData[-trainIndeces,] %>% data.frame 
testSet[,response.index] <- testSet[,response.index] %>% as.numeric %>%
  replace(. == 1, -1) %>% replace(. == 2, 1)

# Defining train and test set for GBP implementation
trainSet.GBM <- bookingsData[trainIndeces,] %>% data.frame
trainSet.GBM[,response.index] <- trainSet.GBM[,response.index] %>% as.numeric %>%
  replace(. == 1, 0) %>% replace(. == 2, 1) # Set y to -1 and 1 for negative and positive cases

testSet.GBM <- bookingsData[-trainIndeces,] %>% data.frame 
testSet.GBM[,response.index] <- testSet.GBM[,response.index] %>% as.numeric %>%
  replace(. == 1, 0) %>% replace(. == 2, 1)

#### Execution ##############################################################

######  Initializing parameters
n_trees <- 100
shrinkage <- 1
maxTreeDepth <- 5

######  Executing own version
fitAdaBoost.res <- fitAdaBoost(trainSet, response.var, n.trees = n_trees,
                               max.depth = maxTreeDepth, debug = 2)
predsAdaBoost <- predictAdaBoost(testSet, response.var, fitAdaBoost.res)
confusionMatrix(predsAdaBoost %>% as.factor, testSet[,response.index] %>%
                  as.factor)

######  Executing package version
fit_gbm <- gbm(is_cancelled ~ ., data = trainSet.GBM, n.trees = n_trees, 
               shrinkage = shrinkage, interaction.depth = maxTreeDepth, 
               distribution = "adaboost")
pred_gbm <- ifelse(predict(fit_gbm, newdata = testSet.GBM, n.trees = n_trees) > 0.5, 1, 0)
confusionMatrix(pred_gbm %>% as.factor, testSet.GBM[,response.index] %>%
                  as.factor)
ROC_own <- roc(testSet[,response.index],predsAdaBoost,
           smoothed = TRUE,
           # arguments for ci
           ci=TRUE, ci.alpha=0.9, stratified=FALSE,
           # arguments for plot
           plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
           print.auc=TRUE, show.thres=TRUE)
sens.ci <- ci.se(ROC)
plot(sens.ci, type="shape", col="lightblue")
plot(sens.ci, type="bars")
