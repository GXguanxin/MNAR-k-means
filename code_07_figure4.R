###############################
########   Figure 4  ##########
###############################

##### Computation time -----------

experiment_computational_time <- function(data,k=3,
                                          mustar=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                          tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])  ){
  time_start <- Sys.time()
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
  time_end <- Sys.time()
  time_tuningpara <- time_end - time_start
  
  time_start <- Sys.time()
  pkpodres <- pkpod_1218(x=data$Missing,k=k,lambda=selectedlambda,nstart=100) #大约0.5s
  time_end <- Sys.time()
  time_pkpod <- time_end - time_start
  
  cer_pkpod <- CER(pkpodres$cluster,data$Misslabel)
  mse_pkpod <- diffcode(mustar,pkpodres$centers)
  
  return(list( selectedlambda=selectedlambda, instability_values=instability_values, instability_values_sd=instability_values_sd,
               cer_pkpod=cer_pkpod, mse_pkpod=mse_pkpod, 
               time_pkpod=time_pkpod, time_tuningpara=time_tuningpara
  ))
}

generate_data <- function(n,p,k,centers=matrix(0,k,p),variance_vec=rep(1,p),
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


# run experiments 

ExperimentResults_computation_time <-list(
  different_p = list( setting = list( n=100*3, k=3, missmechanism="MNAR0", MNAR0lambda_star = 1 ),
                      p_list=c(10, 20, 30, 40, 50), once_time=matrix(0,30,5), tuning_time=matrix(0,30,5) ),
  different_n = list( setting = list( p=10, k=3, missmechanism="MNAR0", MNAR0lambda_star = 1 ),
                      n_list=c(100, 200, 300, 400, 500)*3, once_time=matrix(0,30,5), tuning_time=matrix(0,30,5) )
)
p_t=1
for ( p in c(10, 20, 30, 40, 50) ){
  for (rept in 1:30){
    set.seed(2026+rept*10)
    data <- generate_data(n=100*3, p=p, k=3, 
                                   centers=rbind( c( 0,0 ,rep(0,p-2)) , c( 3,0 ,rep(0,p-2)), c( 1.5, sqrt(6.75), rep(0,p-2) ) ), 
                                   variance_vec = c(rep(1,2),rep(2,p-2)), 
                                   missmechanism="MNAR0", MNAR0lambda_star = 1, nomissing=FALSE )
    res_once <- experiment_computational_time( data=data, k=3,
                                               mustar=rbind( c( 0,0 ,rep(0,p-2)) , c( 3,0 ,rep(0,p-2)), c( 1.5, sqrt(6.75), rep(0,p-2) ) ),
                                               tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])  )
    ExperimentResults_computation_time$different_p$once_time[rept,p_t] <- as.numeric(res_once$time_pkpod)
    ExperimentResults_computation_time$different_p$tuning_time[rept,p_t] <- as.numeric(res_once$time_tuningpara)
    if (rept%%10 == 0){print(rept)}
  }
  p_t <- p_t+1
  print(p)
}


n_t=1
for ( n in c(100, 200, 300, 400, 500)*3 ){
  for (rept in 1:30){
    set.seed(2026+rept*10)
    data <- generate_data(n=n, p=10, k=3, 
                                   centers=rbind( c( 0,0 ,rep(0,8)) , c( 3,0 ,rep(0,8)), c( 1.5, sqrt(6.75), rep(0,8) ) ), 
                                   variance_vec = c(rep(1,2),rep(2,8)), 
                                   missmechanism="MNAR0", MNAR0lambda_star = 1, nomissing=FALSE )
    res_once <- experiment_computational_time( data=data, k=3,
                                               mustar=rbind( c( 0,0 ,rep(0,8)) , c( 3,0 ,rep(0,8)), c( 1.5, sqrt(6.75), rep(0,8) ) ),
                                               tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])  )
    ExperimentResults_computation_time$different_n$once_time[rept,n_t] <- as.numeric(res_once$time_pkpod)
    ExperimentResults_computation_time$different_n$tuning_time[rept,n_t] <- as.numeric(res_once$time_tuningpara)
    if (rept%%10 == 0){print(rept)}
  }
  n_t <- n_t+1
  print(n/3)
}



# summarize data

ExperimentResults_computation_time$different_p$summarizedresults <- cbind(
  "p" = c(10, 20, 30, 40, 50),
  "once_time_mean" = apply(ExperimentResults_computation_time$different_p$once_time,2,mean),
  "once_time_sd" = apply(ExperimentResults_computation_time$different_p$once_time,2,sd),
  "tuning_time_mean" = apply(ExperimentResults_computation_time$different_p$tuning_time,2,mean),
  "tuning_time_sd" = apply(ExperimentResults_computation_time$different_p$tuning_time,2,sd)
)

ExperimentResults_computation_time$different_n$summarizedresults <- cbind(
  "n" = c(100, 200, 300, 400, 500)*3,
  "once_time_mean" = apply(ExperimentResults_computation_time$different_n$once_time,2,mean),
  "once_time_sd" = apply(ExperimentResults_computation_time$different_n$once_time,2,sd),
  "tuning_time_mean" = apply(ExperimentResults_computation_time$different_n$tuning_time,2,mean),
  "tuning_time_sd" = apply(ExperimentResults_computation_time$different_n$tuning_time,2,sd)
)

# draw figures

ExperimentResults_computation_time$different_p$figure <- ggplot(as.data.frame(ExperimentResults_computation_time$different_p$summarizedresults), 
                                                                     aes(x = p ) ) + 
  geom_line(aes(y = once_time_mean)) +
  geom_point(aes(y = once_time_mean), shape = 16, size=2) +
  geom_errorbar( aes(ymin = once_time_mean - once_time_sd, ymax = once_time_mean + once_time_sd), 
                 width = 1 ) +
  theme_bw() + labs(x = "p",y = "Time (s)") + ylim(0,1)


ExperimentResults_computation_time$different_n$figure <- ggplot(as.data.frame(ExperimentResults_computation_time$different_n$summarizedresults), 
                                                                     aes(x = n ) ) + 
  geom_line(aes(y = once_time_mean)) +
  geom_point(aes(y = once_time_mean), shape = 16, size=2) +
  geom_errorbar( aes(ymin = once_time_mean - once_time_sd, ymax = once_time_mean + once_time_sd), 
                 width = 30 ) +
  scale_x_continuous(breaks = c(100, 200, 300, 400, 500)*3 ) + 
  theme_bw() + labs(x = "n",y = "Time (s)") + ylim(0,1.5)


ExperimentResults_computation_time$different_p$figure
ExperimentResults_computation_time$different_n$figure

