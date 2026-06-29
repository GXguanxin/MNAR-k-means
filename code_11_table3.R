#################################
########   Table 3  #############
#################################

generate_data_real_data_based_missing <- function(x, clusterlabel,  
                                                  missmechanism=c("MCAR","MAR","MNAR1","MNAR2","MNAR0"),
                                                  MCARmissprob=rep(0,p),
                                                  MARpsi=c(0.1,0.1),
                                                  MNAR1phi=c(3,2),
                                                  MNAR2percent=rep(0.1,p),
                                                  MNAR0lambda_star=0.1,
                                                  nomissing=FALSE){
  n=dim(x)[1]
  p=dim(x)[2]
  k=length(unique(clusterlabel))
  
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


realdata <- apply(t(as.matrix(lymphoma.x)), 2 , function(z) (z-mean(z,na.rm=TRUE))/sd(z,na.rm=TRUE)    )

experiment_lymphoma <- function(realdata, clusterlabel=lymphoma.y$V1, k=3,
                                missmechanism=c("MNAR1","MNAR2","MNAR0"),
                                missproportion=c(0.1,0.3,0.5),
                                tuninglambdalist=c(0.1,seq(1,10,by=1)) ){
  #generate missing data
  if (missmechanism=="MNAR0"){
    if (missproportion==0.1){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR0",
                                                    MNAR0lambda_star = 50)
    }
    if (missproportion==0.3){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR0",
                                                    MNAR0lambda_star = 5)
    }
    if (missproportion==0.5){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR0",
                                                    MNAR0lambda_star = 1.5)
    }
  }
  
  if (missmechanism=="MNAR1"){
    if (missproportion==0.1){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR1",MNAR1phi = c(-15,0))
    }
    if (missproportion==0.3){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR1",MNAR1phi = c(-1.5,0))
    }
    if (missproportion==0.5){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR1",MNAR1phi = c(-0.1,0))
    }
  }
  
  if (missmechanism=="MNAR2"){
    if (missproportion==0.1){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR2",
                                                    MNAR2percent = rep(0.1,ncol(realdata)) )
    }
    if (missproportion==0.3){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR2",
                                                    MNAR2percent = rep(0.3,ncol(realdata)) )
    }
    if (missproportion==0.5){
      data <- generate_data_real_data_based_missing(x=realdata, clusterlabel = clusterlabel,
                                                    missmechanism = "MNAR2",
                                                    MNAR2percent = rep(0.5,ncol(realdata)) )
    }
  }
  
  #run different methods
  missdata_n <- dim(data$Missing)[1]
  missdata_proportion <- sum(is.na(data$Missing))/(nrow(data$Missing)*ncol(data$Missing))
  
  data$meanimpt <- data$Missing
  data$meanimpt[is.na(data$Missing)] <- matrix(colMeans(data$Missing,na.rm=TRUE),62,4026,byrow = TRUE)[is.na(data$Missing)]
  meanimptres <- kmeans(data$meanimpt,centers=3,iter.max = 100,nstart=100)
  cer_meanimpt <- CER(meanimptres$cluster,data$Misslabel)
  
  data$knnimpt <- t(impute.knn(t(as.matrix(data$Missing)), k = 10, rng.seed = 2026*88)$data)  #10-nearest-neighbors
  knnimptres <- kmeans(data$knnimpt,centers=k,iter.max = 100,nstart=100)
  cer_knnimpt <- CER(knnimptres$cluster,data$Misslabel)
  
  kpodres <- pkpod_1218(x=data$Missing,k=k,lambda=0.000001,nstart=100)  
  cer_kpod <- CER(kpodres$cluster,data$Misslabel)
  
  instability_values <- rep(0,length(tuninglambdalist))
  instability_values_sd <- rep(0,length(tuninglambdalist))
  for (t in 1:length(tuninglambdalist)){
    instability_res <- instability_1218(x=data$Missing,k=k,tuningpara = tuninglambdalist[t], 
                                        clusteringfun = pkpod_1218, nstart=10, perm.max = 10  ) 
    instability_values[t] <- instability_res$Svalue
    instability_values_sd[t] <- instability_res$Svalue_sd
  }
  selectedlambda <- tuninglambdalist[which.min(instability_values)]
  if ( selectedlambda==tail(tuninglambdalist,1) ){
    lowerbound <- tail(instability_values,1) + tail(instability_values_sd,1)
    selectedlambda <- tuninglambdalist[min(which( instability_values < lowerbound ))]
  }
  
  pkpodres <- pkpod_1218(x=data$Missing,k=k,lambda=selectedlambda,nstart=100) 
  cer_pkpod <- CER(pkpodres$cluster,data$Misslabel)
  
  return(list( missdata_n=missdata_n,  missdata_proportion=missdata_proportion, selectedlambda=selectedlambda,
               cer_meanimpt=cer_meanimpt, cer_knnimpt=cer_knnimpt, cer_kpod=cer_kpod, cer_pkpod=cer_pkpod,
               instability_values=instability_values, instability_values_sd=instability_values_sd
  ))
}


### experiments -------------

ExperimentResults_realdata_lymphoma <- list(
  allsettings=list(missmechanism = c("MNAR0","MNAR1","MNAR2"),
                   missproportion = c(0.1,0.3,0.5),
                   tuninglambdalist=seq(2,10,by=2),
                   realdata = realdata, clusterlabel = lymphoma.y$V1, k=3),
  allresults=list()
)
setting=1
for (missmechanism in c("MNAR0","MNAR1","MNAR2")){
  for (missproportion in c(0.1,0.3,0.5) ){
    res_mat=matrix(0,10,10)
    for (rept in 1:10){
      set.seed(2026+rept)
      res_once <- experiment_lymphoma(realdata = realdata, clusterlabel = lymphoma.y$V1, k=3,
                                      missmechanism = missmechanism, 
                                      missproportion = missproportion,
                                      tuninglambdalist=seq(2,10,by=2) )
      res_mat[rept,1:4] <- c(res_once$cer_meanimpt,res_once$cer_knnimpt,res_once$cer_kpod,res_once$cer_pkpod)
      res_mat[rept,5:10] <- c(res_once$selectedlambda,res_once$instability_values)
    }
    ExperimentResults_realdata_lymphoma$allresults[[setting]] <- res_mat
    setting = setting + 1
    print(c(missmechanism,missproportion,res_once$missdata_proportion))
  }
}


### summarize results -------

for (setting in 1:9){
  print(apply(ExperimentResults_realdata_lymphoma$allresults[[setting]][,1:4],2,mean))
}

for (setting in 1:9){
  print(apply(ExperimentResults_realdata_lymphoma$allresults[[setting]][,1:4],2,sd))
}

