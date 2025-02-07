library(pacman)
p_load(PMA, tidyverse, fastDummies, softImpute, vtreat, bcv, ggfortify, factoextra, ggplot2, psych,dplyr)
set.seed(1)

# PLEASE DONT CONFIG YOUR WD BUT USE THE R-PROJECT 
# setwd('C:\\Users\\Rohan Chatterjee\\Desktop\\Unsupervised ML\\Week 3')

#### Data ####
load('./Week3/FIFA2017_NL.rdata')

fifa$id <- 1:nrow(fifa)
fifa_scaled <- fifa %>% select(-c(name, Position, club)) %>% scale
imputer <- softImpute(fifa_scaled, rank.max = 33, maxit = 1000)
fifa_imputed <- complete(fifa_scaled, imputer) %>% as_tibble

## Summary Statistics
FW <- subset(fifa, fifa$Position == "FW")
Mid <- subset(fifa, fifa$Position == "Mid")
Def <- subset(fifa, fifa$Position == "Def")
Gk <- subset(fifa, fifa$Position == "Gk")

summaryFifa <- summary(fifa)
summaryFW <- summary(FW)
summaryMid <- summary(Mid)
summaryDef <- summary(Def)
summaryGk <- summary(Gk)

summaryFW
summaryMid
summaryDef
summaryGk
#### Standard code ####
model_spc <- SPC(fifa_imputed %>% as.matrix, sumabsv = sqrt(ncol(fifa_imputed)), K=2, center =F , trace = F)
model_pca <- prcomp(fifa_imputed, center = TRUE, scale. = TRUE, rank. = 2)

# Comparison of models
abs(model_pca$rotation - model_spc$v) %>% ggplot(aes(x=PC1, y=PC2, label=rownames(model_pca$rotation))) + geom_point() + xlab('Delta PC1')+ ylab('Delta PC2')

model_pca %>% autoplot(., data = fifa, colour = 'Position', loadings = F)
model_pca %>% autoplot(., data = fifa, colour = 'club', loadings = F)



fviz_eig(model_pca, addlabels = TRUE)
print(model_pca$rotation[,1])
print(model_spc$v[,1])

model_pca$rotation %>% ggplot(aes)

var <- get_pca_var(model_pca)
fviz_pca_var(model_pca, col.var = "black")

fifa_simple <- predict(model_pca) %*% t(model_pca$rotation) %>% as_tibble
fifa_simple$y <- fifa$Position == 'Gk'

model_pca$


#### Own function ####
soft_threshold <- function(vX_t_u,c){
  
  vS <- vX_t_u
  
  
  vv <- vS / ( sqrt( t( vS ) %*% vS)[1,1] ) 
  
  
  if(sum(vv) > c){
   message('sum(vv) was big')
   colnames <- 'V1'
   see <- vS%>% as.data.frame()%>%setNames(colnames)%>% 
     mutate(sign = sign(V1)) %>% group_by(sign) %>% summarise(count = n(), addition = sum(V1))
    
   lambda <- (c + see$addition[2] - see$addition[1]) / (see$count[1] - see$count[2] )
   
   vS <- (vX_t_u) * ( abs(vX_t_u) - rep(lambda, nrow(vX_t_u)) )

   vv <- vS / ( sqrt( t( vS ) %*% vS)[1,1] ) 
   
  } 
  
  return(vv)
  
  
}


sparse_pca <- function(mX,imax=20,vv=NULL){
  c <- sqrt(ncol(mX))
  
  if(is.null(vv)){
    vv <- eigen( t(mX) %*% mX )$vectors[,1]
  }
  
  i=1
  
  while(i<=imax){
    
    i <- i+1 
    
    vXv <- mX %*% vv
    
    
    
    vu <- vXv  / (sqrt(t( vXv ) %*% vXv )[1,1])
    
    vX_t_u <- t(mX) %*% vu
    
    vv <- soft_threshold(vX_t_u=vX_t_u,c=c)
  }
  
  sigma <- t(vu) %*% mX %*% vv
  
  return(list(u=vu,v=vv,sigma=sigma))
  
}

pca_own <- sparse_pca(fifa_imputed %>% as.matrix)
pca_own$v
model_spc$v

sum(check$u - see$u) #very close
check$d - see$sigma #very cl0se
sum(check$v - see$v) #very close


