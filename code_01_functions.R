#############################################
########   MNAR-k-means functions  ##########
#############################################

library(kpodclustr)
library(MASS)
library(Rfast)
library(Rcpp)
library(RcppArmadillo) 
library(abind)
library(RSKC)
library(mice)
library(impute)
library(dplyr)
library(tidyr)
library(ggplot2)
library(umap)

sourceCpp("initcenters_oncepkpod_functions.cpp")

#---------------------------------------------

pkpod_1218<-function(x,k,lambda,itermax=100,nstart=100){
  n=dim(x)[1]
  p=dim(x)[2]
  
  
  x_imputed <- x
  mu_mat <- matrix(colMeans(x,na.rm = TRUE), nrow = n, ncol = p, byrow = TRUE)
  x_imputed[is.na(x)==1] <- mu_mat[is.na(x)==1] 
  #check NA columns in x_imputed 
  NAcols <- which( apply(x_imputed, 2, function(z) sum(is.na(z)) ) == n )
  if (length(NAcols)>0){x_imputed[,NAcols] <- rep(0,n)}
  
  obj_old=Inf
  obj_nstart=rep(0,nstart)
  for (start in 1:nstart){
    #initialization
    tmp_centers <- kmpp_init_centers_cpp(x_imputed,k)$centers
    init_centers <- kmeans(x_imputed, centers = tmp_centers, nstart = 1, iter.max = 100)$centers #for fixed init centers, nstart must =1
    res=pkpod_once_1218cpp(x=x,k=k,lambda=lambda,init_centers=init_centers, itermax = itermax)
    obj=res$obj_val
    obj_nstart[start]=obj
    #compare
    if(obj<obj_old){
      obj_old=obj
      result=res
    }
  }
  result$obj_nstart=obj_nstart
  return(result)
}

instability_1218<-function(x,tuningpara,k,clusteringfun,nstart=1,perm.max=10){
  f <- function(x, y) sum(  (x - y)^2, na.rm = TRUE )
  n=dim(x)[1]
  p=dim(x)[2]
  
  Svec=rep(0,perm.max)
  for (perm in 1:perm.max){
    
    #split datasets
    setindex=sample(1:3,n,replace = TRUE)
    train1=x[which(setindex==1),]
    train2=x[which(setindex==2),]
    testing=x[which(setindex==3),]
    ntesting=nrow(testing)
    
    #run clustering on train1
    res1=clusteringfun(x=train1,k=k,lambda=tuningpara,nstart=nstart)
    estimated_mu1=res1$centers
    
    #run clustering on train2
    res2=clusteringfun(x=train2,k=k,lambda=tuningpara,nstart=nstart)
    estimated_mu2=res2$centers
    
    #predict testing
    missing_testing=is.na(testing)
    dist_crossdis1 <- proxy::dist(testing, estimated_mu1, f)  #squared distance in observed dimensions
    dist_mat1 <- `dim<-`(c(dist_crossdis1), dim(dist_crossdis1))
    Q_umat1 <- dist_mat1 + tuningpara*(missing_testing%*%t(estimated_mu1^2)) #penalized distance
    cluster1 <- Rfast::rowMins(Q_umat1)
    umat1 <- matrix(0, ntesting , k)
    for(kk in 1:k){
      umat1[which(cluster1==kk),kk] <- 1
    }
    dist_crossdis2 <- proxy::dist(testing, estimated_mu2, f)  #squared distance in observed dimensions
    dist_mat2 <- `dim<-`(c(dist_crossdis2), dim(dist_crossdis2))
    Q_umat2 <- dist_mat2 + tuningpara*(missing_testing%*%t(estimated_mu2^2)) #penalized distance
    cluster2 <- Rfast::rowMins(Q_umat2)
    umat2 <- matrix(0, ntesting, k)
    for(kk in 1:k){
      umat2[which(cluster2==kk),kk] <- 1
    }
    
    
    #calculate instability
    Vmat=( umat1 %*% t(umat1) ) + ( umat2 %*% t(umat2) )
    Svec[perm]=sum(Vmat[upper.tri(Vmat)]==1)/((ntesting-1)*(ntesting-2)/2)
  }
  
  return(list(Svalue=mean(Svec),Svalue_sd=sd(Svec),Svec=Svec,estimated_mu1=estimated_mu1,estimated_mu2=estimated_mu2))
  
}

fun_mice <- function(x, k=2){
  miceimputeres <- mice(data=x, m = 5, maxit = 10,method = "pmm", printFlag=FALSE)
  micedata <- matrix(0,nrow(x),ncol(x))
  for (mm in 1:5){
    micedata <- micedata + mice::complete(miceimputeres, action = mm)
  }
  micedata <- micedata/5
  miceres <- kmeans(x=micedata,centers=k,iter.max = 100,nstart=100)
  return(miceres)
}

diffcode<-function(codebook0,codebook){
  K=dim(codebook0)[1]
  p=dim(codebook0)[2]
  
  mse=rep(0,K)
  for (k in 1:K){
    dist=apply(codebook,1,function(z) sum((codebook0[k,]-z)^2)  )
    mse[k]=min(dist)
  }
  return(sum(mse))
}
