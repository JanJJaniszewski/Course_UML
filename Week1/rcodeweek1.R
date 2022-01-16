## ---- assignment1
#library(dsmle)
#library(MASS)
#library(kernlab)

set.seed(1234)

# construct data
load("bank.RData")

data_idx <- sample(nrow(bank), size=1000)
data <- bank[data_idx, -c(16, 19)]

# make train test
y <- subset(data, select=y)
X <- model.matrix(y~.,data=data)

dim_data <- dim(X)
n <- dim_data[1]
p <- dim_data[2]

n.train <- 0.7 * n
idx <- sample(nrow(X), size=n.train)

X.train <- X[idx,]
X.train[, -c(1, 20, 25)] <- scale(X.train[, -c(1, 20, 25)])
y.train <- y[idx,]

X.test <- X[-idx,]
y.test <- y[-idx,]

levels(y.test) <- c(-1, 1)
levels(y.train) <- c(-1,1)
y.test <- as.numeric(levels(y.test))[y.test]
y.train <- as.numeric(levels(y.train))[y.train]

#start svmmj
SVMloss <- function(X, y, v, lambda){
  q <- X %*% v
  w <- v[-1]
  loss <- sum(pmax(0, 1- y*q)) + lambda * t(w) %*% w
  return(loss)
}

calc_a_i <- function(y_i, q_i, eps){
  return(ifelse(y_i == 1, 1/(4*pmax(abs(1-q_i), eps)), 
  1/(4*pmax(abs(q_i + 1), eps)) ))
}

calc_b_i <- function(y_i, a_i){
  return(ifelse(y_i == 1, a_i + 0.25, -(a_i + 0.25) ))
}

SVMmaj_own <- function(X, y, eps, lambda=20){
  k <- 1
  vP_old <- rep(1, p)
  vP_new <- vP_old
  loss_old <- SVMloss(X, y, vP_old, lambda)
  loss_new <- loss_old
  b <- mapply(calc_b_i, y, A_asvector)
  lambdaP <- lambda * diag(nrow=ncol(X))
  P[1,1] <- 0
  
  while(k==1 | (loss_old - loss_new)/loss_old > eps){
    #print((loss_old - loss_new)/loss_old)
    k <- k + 1
    loss_old <- loss_new
    vP_old <- vP_new

    A_asvector <- mapply(calc_a_i, y, X %*% vP_old, eps)
    A <- diag(A_asvector)


    vP_new <- solve(t(X)%*%A%*%X + lambdaP, t(X)%*%b)
    loss_new <- SVMloss(X, y, vP_new, lambda)
  }
  print(k)
  print(loss_new)
}

eps <- 0.0000001
SVMmaj_own(X.train, y.train, eps)
#model_compare <- svmmaj(X.train, y.train, lambda = 20, scale = 'zscore', 
#hinge = 'absolute')
#print(model_compare$loss)

#y.train <- as.factor(y.train)
#y.test <- as.factor(y.test)

#weights.obs <- list(no = length(which(y.train == "yes"))/length(y.train), 
#yes = length(which(y.train == "no"))/length(y.train))
#model <- svmmaj(X.train, y.train, lambda = 2, scale = 'interval', 
#hinge = 'absolute', kernel = rbfdot, kernel.sigma = 1,  
#weights.obs = weights.obs)
#print(model$loss)
#predict(model, X.new = X.test, y = y.test, show.plot = TRUE)

#model_linear <- svmmaj(X.train, y.train, lambda = 2, scale = 'interval', 
#hinge = 'absolute', weights.obs = weights.obs)
#print(model_linear$loss)
#predict(model_linear, X.new = X.test, y = y.test, show.plot = TRUE)
#plotWeights(model, plotdim = c(2, 4))


