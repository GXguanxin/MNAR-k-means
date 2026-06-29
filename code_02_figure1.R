###################################
########   Figure 1   #############
###################################


# generate data

generate_data_intro <- function(n,p,k,centers=matrix(0,k,p),variance_vec=rep(1,p),
                                   missmechanism=c("MCAR","MAR","MNAR1","MNAR2","MNAR0"),
                                   MCARmissprob=rep(0,p),
                                   MARpsi=c(0.1,0.1),
                                   MNAR1phi=c(3,2),
                                   MNAR2percent=rep(0.1,p),
                                   MNAR0lambda_star=0.1,
                                   nomissing=FALSE                    ){
  ##original data
  x <- matrix(0, nrow = n, ncol = p)
  clusterlabel <- sample(1:k, n, replace = TRUE)
  for (kk in 1:k){
    kk_cluster_size <- sum(clusterlabel==kk)
    x_kk_cluster <- MASS::mvrnorm(kk_cluster_size, mu=centers[kk,], Sigma=diag(variance_vec) ) 
    x[which(clusterlabel==kk),] <- x_kk_cluster
  }
  
  if(nomissing==TRUE){
    return(list( Orig=x, Missing=NULL, CompleteCase=x, 
                 Origlabel=clusterlabel, Misslabel=NULL, Complabel=clusterlabel  ))
  }
  
  ##missingness
  
  if (missmechanism[1] == "MCAR"){
    miss_mat <- matrix(0,n,p)
    for (j in 1:p){
      miss_mat[,j] <- rbinom(n, 1, prob=MCARmissprob[j])
    }
  }
  
  if (missmechanism[1] == "MAR"){
    missprob=rep(0,n)
    miss_mat=matrix(0,n,p)
    for (i in 1:n){
      missprob[i]=1/(1+exp(- (MARpsi[1])*( x[i,1] - (MARpsi[2]) ) ) )
      for (j in 2:p){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i])
      }
    }
  }
  
  if (missmechanism[1] == "MNAR1"){
    missprob=1/(1+exp(- (MNAR1phi[1])*( x^2 - (MNAR1phi[2]) ) ) )
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      for (i in 1:n){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i,j])
      }
    }
  }
  
  if (missmechanism[1] == "MNAR2"){
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      threshold=quantile(abs(x[,j]),probs=MNAR2percent[j])
      miss_mat[,j]=( abs(x[,j])< threshold )
    }
  }
  
  if (missmechanism[1] == "MNAR0"){
    missprob=exp(-MNAR0lambda_star*(x^2) )
    miss_mat=matrix(0,n,p)
    for (j in 1:p){
      for (i in 1:n){
        miss_mat[i,j]=rbinom(n=1,size=1,prob=missprob[i,j])
      }
    }
  }
  
  ##missing data and complete data
  x_tmp <- x
  x_tmp[miss_mat==1] <- NA
  nullrows <- which( (Rfast::rowsums(miss_mat)) == p )
  comprows <- which( (Rfast::rowsums(miss_mat)) == 0 )
  is_null_row <- ( (Rfast::rowsums(miss_mat)) == p )
  
  if(length(nullrows)==0){
    x_miss <- x_tmp
    clusterlabel_miss <- clusterlabel
  }else{
    x_miss <- x_tmp[-nullrows,]
    clusterlabel_miss <- clusterlabel[-nullrows]
  }
  
  if(length(comprows)==0){
    x_comp <- NULL
    clusterlabel_comp <- NULL
  }else{
    x_comp <- x_tmp[comprows,]
    clusterlabel_comp <- clusterlabel[comprows]
  }
  
  return(list(Orig=x, Missing=x_miss, CompleteCase=x_comp, 
              Origlabel=clusterlabel, Misslabel=clusterlabel_miss, Complabel=clusterlabel_comp,
              is_null_row=is_null_row  ))
}

set.seed(2026*8)
mustar=rbind( c(0,0) , c(3,0), c(3/2,3*sqrt(3)/2) )
data <- generate_data_intro(n=100*3, p=2, k=3, 
                               centers=mustar, 
                               variance_vec = c(rep(1,2),rep(2,0)), 
                               missmechanism="MNAR0", MNAR0lambda_star = 1, nomissing=FALSE )

# run different methods

kmres <- kmeans(data$Orig[-which(data$is_null_row==1),],centers=3,iter.max = 100,nstart=100)
CER(kmres$cluster,data$Misslabel)
diffcode(mustar,kmres$centers)

kpodres <- pkpod_1218(x=data$Missing,k=3,lambda=0,nstart=100)  
CER(kpodres$cluster,data$Misslabel)
diffcode(mustar,kpodres$centers)

pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=1,nstart=100) 
CER(pkpodres$cluster,data$Misslabel)
diffcode(mustar,pkpodres$centers)


# draw figures

col.vecs <- matrix(0,nrow(data$Missing),3)
col.vecs[,1] <- kmres$cluster
map <- apply(table(kpodres$cluster,kmres$cluster), 2, function(z) names(which.max(z)))
col.vecs[,2] <- map[kpodres$cluster]
map <- apply(table(pkpodres$cluster,kmres$cluster), 2, function(z) names(which.max(z)))
col.vecs[,3] <- map[pkpodres$cluster]

col_map <- c("1"="red", "2"="blue", "3"="green4")
pch_map <- c("1"=0, "2"=1, "3"=2)

plot(data$Orig[-which(data$is_null_row==1),],
     col=col_map[col.vecs[,1]], pch=pch_map[data$Misslabel], cex=0.8, xlab="X1",ylab="X2"  )
points(kmres$centers,col="black",pch=c(15,17,16),cex=1.5 )

plot(data$Orig[-which(data$is_null_row==1),],
     col=col_map[col.vecs[,2]], pch=pch_map[data$Misslabel], cex=0.8, xlab="X1",ylab="X2"  )
points(kpodres$centers,col="black",pch= c(17,16,15),cex=1.5 )

plot(data$Orig[-which(data$is_null_row==1),],
     col=col_map[col.vecs[,3]], pch=pch_map[data$Misslabel], cex=0.8, xlab="X1",ylab="X2"  )
points(pkpodres$centers,col="black",pch= c(15,17,16),cex=1.5 )
