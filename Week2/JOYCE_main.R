
################################################################# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, gbm, vtreat, kernlab, caret, rpart,xgboost,ohenery)

######################################################### Load and prepare data 
load("~/GitHub/Course_UML/Week2/bookings.RData")

##### Analyze #####
bookings %>% group_by(is_cancelled) %>% summarise_all(mean) %>% head()

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
X.train <- trainSet %>% select(-is_cancelled) %>% scale %>% as_tibble
X.test <- testSet %>% select(-is_cancelled) %>% scale%>% as_tibble

# Set y to -1 and 1 for negative and positive cases
y.train <- trainSet %>% pull(is_cancelled) %>% as.numeric %>% replace(.==1, -1) %>% replace(.==2, 1)
y.test <- testSet %>% pull(is_cancelled) %>% as.numeric %>% replace(.==1, -1) %>% replace(.==2, 1)

######################################################## Build adaboost
adaboost_own <- function(X, y, formula, iteration, maxdepth=4){
  
  #set df 
  df <- cbind(y,X)
  df <- as.data.frame(df)
  #give formula 
  formula <- lm(y~., data = df)
  #set initial value of the weight
  n <- nrow(X)
  vWeight <- rep(1/n,n)
  
  #make an array to storage the alpha 
  vAlpha_list <- c()
  
  #make a list to storage the classifier fitting model
  vTree_list <- list()
  
  # create a list to storage the prediction
  vY_predict_list <- list()
  
  #into interations
  for (m in 1:iteration){
    # 1 fit classifier with tree
    tree_m <- rpart(formula =formula,
                  data = X.train,
                  weights = vWeight,
                  method = "class",
                  control =rpart.control(maxdepth = maxdepth, cp = 0,xval = 0, 
                                         minsplit = 2, minbucket = 1))

    #append new tree to trees collection
    vTree_list[[m]] <- tree_m
    
    #collect prediction results for calculating error rate
    vY_predict_m <- predict(tree_m, data = X.train, type = "class")
    vY_predict_list[[m]] <- vY_predict_m
    
    # 2 calculate weighted error rate
    error=0
    index_m <- as.numeric(y != vY_predict_m)
    err_m <- sum(vWeight * index_m)
    
    # 3 calculate alpha
    alpha_m <- log((1 - err_m)/err_m)
    vAlpha_list[m] <- alpha_m
                  
    # 4 set weight
    vWeight <- vWeight * exp(alpha_m * index_m)
    
    # 5 normalize the weights to sum to 1
    vWeight <- vWeight %>% normalize
  }
  #end the loop and return output
  vResult <- data.frame(row.names = seq(1, nrow(X)))
  vResult[m] <- vAlpha_list[[m]] * vY_predict_list[[m]]
  vResult_final <- apply(X = vResult, MARGIN = 1, FUN = sum) %>% sign
  return(vResult_final)
}
ada_own <-adaboost_own(X=X.train, y=y.train, formula=formula, iteration=30)

