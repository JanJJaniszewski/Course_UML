library(pacman)
p_load(PMA, tidyverse, fastDummies, softImpute, vtreat, bcv, ggfortify, factoextra, ggplot2, psych,dplyr)
set.seed(1)

# PLEASE DONT CONFIG YOUR WD BUT USE THE R-PROJECT 
# setwd('C:\\Users\\Rohan Chatterjee\\Desktop\\Unsupervised ML\\Week 3')

#### Data ####
load('./Week3/FIFA2017_NL.rdata')

fifa_imputed <- fifa %>% select(-c(name, Position, club, contains('eur_'))) %>% scale %>% as_tibble

#### Standard code ####
# SPCA predictions
model_spc <- SPC(fifa_imputed %>% as.matrix, sumabsv = 2, K=2, center =F , trace = F)
predictions <- tibble(.rows = nrow(fifa_imputed))
Phi <- model_spc$v
predictions$pc1 <- as.matrix(fifa_imputed) %*% Phi[, 1] 
predictions$pc2 <- as.matrix(fifa_imputed) %*% Phi[, 2] 
predictions$position <- fifa$Position
predictions %>% ggplot(aes(x=pc1, y=pc2, color=position)) + geom_point() + xlim(-3, 5.5) + ylim(-3,5.5)

# Classical PCA predictions
model_pca <- prcomp(fifa_imputed, center = TRUE, scale. = TRUE, rank. = 2)


fviz_eig(model_pca, addlabels = TRUE)

# Comparison of models
abs(model_pca$rotation - model_spc$v) %>% ggplot(aes(x=PC1, y=PC2, label=rownames(model_pca$rotation))) + geom_point() + xlab('Delta PC1')+ ylab('Delta PC2')

# Loading analysis
model_pca %>% autoplot(., data = fifa, colour = 'Position', loadings = F, xlim=c(-0.1, 0.14), ylim=c(-0.1, 0.14))

var <- get_pca_var(model_pca)
fviz_pca_var(model_pca, col.var = "black")
model_pca$rotation

# Plots
fifa_output <- fifa_imputed
fifa_output$PC1 <- as.matrix(fifa_imputed) %*% model_pca$rotation[,1]
fifa_output$PC2 <- as.matrix(fifa_imputed) %*% model_pca$rotation[,2]
fifa_output$goalkeeper <- fifa$Position == 'Gk'
fifa_output$position <- fifa$Position

fifa_output %>% group_by(position) %>% summarise(PC1=mean(PC1), PC2=mean(PC2))

t1 <- fifa_output %>% group_by(goalkeeper) %>% summarise_all(mean) %>% select(-c(goalkeeper, PC1, PC2)) %>% t 
rn <- rownames(t1)
colnames(t1) <- c('Other', 'Goalkeeper')
t1 <- t1 %>% as_tibble 
t1$Skill <- rn
t1 <- t1 %>% gather(key='Position', value='Value', 1:2)
t1 %>% ggplot(aes(x=Skill, y=Value, color=Position)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


t1 <- fifa_output %>% group_by(position) %>% summarise_all(mean) %>% select(-c(goalkeeper, PC1, PC2, position)) %>% t 
rn <- rownames(t1)
colnames(t1) <- c('FW', 'Mid' , 'Def','Gk')
t1 <- t1 %>% as_tibble 
t1$Skill <- rn
t1 <- t1 %>% gather(key='Position', value='Value', 1:3)
t1 %>% ggplot(aes(x=Skill, y=Value, color=Position)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


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


