#### Load packages and init ####################################################
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, vtreat, caret, rpart, ohenery)
set.seed(42)


#### Data #################################################################
##### Load #####
load('Week2/bookings.RData')

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



#### AdaBoost ####

##### Init #####
give_ada_booster <- function(X, y, M=10, classifier_method='rpart', verbose=0){
  w <- rep(1/nrow(X.train), nrow(X.train))
  G <- list()
  M <- seq(1, M)
  alpha <- list()
  
  if(verbose >= 1) print('Ada gets boosted')
  for (m in M){
    G[[m]] <- train(x = X, y=y, method = classifier_method, weights = w)
    preds_Gm <- sign(predict(G[[m]], X))
    err_m <- sum(w * (y != preds_Gm)) / sum(w)
    alpha[[m]] <- log((1 - err_m) / err_m)
    w <-  (w * exp(alpha[[m]] * (y != preds_Gm))) %>% normalize
    
    if (verbose>=2){
      print('------------')
      print('Run')
      print(m)
      print('SD(w)')
      print(sd(w))
      print('Alpha')
      print(alpha[[m]])
      print('Error')
      print(err_m)
    }
  }
  ada <- list()
  ada[['alpha']] <- alpha
  ada[['G']] <- G
  return(ada)
}

ada <- give_ada_booster(X.train, y.train, M=30, classifier_method = 'rpart', verbose=2)

ask_boosted_ada_to_predict <- function(X, ada, verbose=0){
  M <- ada$alpha %>% length
  M <- seq(1, M)
  alpha <- ada$alpha
  G <- ada$G
  preds <- data.frame(row.names = seq(1, nrow(X)))
  
  if(verbose >= 1) print('Ada uses her incredible booster powers to predict outcomes')
  for(m in M){
    preds[m] <- alpha[[m]] * predict(G[[m]], X)
    
  }
  preds_final <- apply(X = preds, MARGIN = 1, FUN = sum) %>% sign
  
  return(preds_final)
}

preds <- ask_boosted_ada_to_predict(X=X.test, ada=ada)
confusionMatrix(preds %>% as.factor, y.test %>% as.factor)