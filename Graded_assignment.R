library(PMA)
library(dplyr)
library(fastDummies)
library(softImpute)
setwd('C:\\Users\\Rohan Chatterjee\\Desktop\\Unsupervised ML\\Week 3')


load('FIFA2017_NL.rdata')


mX <- matrix(c(1,2,3,4,5,6,7,8,9,10),nrow = 2,ncol=5)

fifa_dummy <- dummy_cols(fifa,remove_first_dummy = T)%>% select(!name&!club&!Position) %>% as.matrix()

imputed_fifa <- softImpute(fifa_dummy, rank.max =  min(nrow(fifa_dummy),ncol(fifa_dummy)),
                           lambda = 1.5,type='svd')

imputed_fifa <- complete(fifa_dummy,imputed_fifa)



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

see <- sparse_pca(imputed_fifa)

check <- SPC(imputed_fifa,sumabsv = sqrt(ncol(imputed_fifa)), K=1,center =F , trace = F)

sum(check$u - see$u) #very close
check$d - see$sigma #very cl0se
sum(check$v - see$v) #very close
