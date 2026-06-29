###################################
########   Figure 2   #############
###################################

### define functions ------------
generate_data_MNAR0_p2 <- function(n,p,k,centers=matrix(0,k,p),variance_vec=rep(1,p),
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

experiment_MNAR0_p2 <- function(a=3,k=3,n=100*3,
                                mustar=rbind( c( 0,0 ,rep(0,0)) , c( 3,0 ,rep(0,0)), c( 1.5, sqrt(6.75), rep(0,0) ) ),
                                MNAR0lambda_star=1,
                                tuninglambdalist=c(0.01, seq(0,2,by=0.2)[-1])  ){
  data <- generate_data_MNAR0_p2(n=n, p=2, k=k, 
                                 centers=mustar, 
                                 variance_vec = c(rep(1,2),rep(2,0)), 
                                 missmechanism="MNAR0", MNAR0lambda_star = MNAR0lambda_star, nomissing=FALSE )
  missdata_n <- dim(data$Missing)[1]
  missdata_proportion <- sum(is.na(data$Missing))/(nrow(data$Missing)*ncol(data$Missing))
  
  kmres <- kmeans(data$Orig,centers=k,iter.max = 100,nstart=100)
  cer_km <- CER(kmres$cluster,data$Origlabel)
  mse_km <- diffcode(mustar,kmres$centers)
  
  kpodres <- pkpod_1218(x=data$Missing,k=k,lambda=0,nstart=100)  
  cer_kpod <- CER(kpodres$cluster,data$Misslabel)
  mse_kpod <- diffcode(mustar,kpodres$centers)
  
  
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
  mse_pkpod <- diffcode(mustar,pkpodres$centers)
  
  return(list( missdata_n=missdata_n,  missdata_proportion=missdata_proportion, selectedlambda=selectedlambda,
               cer_km=cer_km, cer_kpod=cer_kpod, cer_pkpod=cer_pkpod,
               mse_km=mse_km, mse_kpod=mse_kpod, mse_pkpod=mse_pkpod,
               instability_values=instability_values, instability_values_sd=instability_values_sd
  ))
}


### experiments -------------------

ExperimentResults_MNAR0_p2 <- list(
  allsettings=list(
    alist=c(3,4), klist=c(2,3,4),
    MNAR0lambda_star_list=c(0.1, 0.5, 1, 2, 10),
    samplesize_eachcluster_list=c(200, 100, 100, 100, 100),
    tuninglambdalist=c(0.01, seq(0,2,by=0.2)[-1])
  ),
  allresults=list()
)

tuninglambdalist=c(0.01, seq(0,2,by=0.2)[-1])
setting=1
for (a in c(2,3,4)){
  for (k in c(2,3,4) ){
    if (k==2){mustar=rbind( c(0,0) , c(a,0) ) }
    if (k==3){mustar=rbind( c(0,0) , c(a,0), c(a/2,a*sqrt(3)/2) ) }
    if (k==4){mustar=rbind( c(0,0) , c(a,0), c(0,a) , c(a,a) ) }
    for (MNAR0lambda_star in c(0.1, 0.5, 1, 2, 10) ){
      if (MNAR0lambda_star==0.1){samplesize_eachcluster=200}else{samplesize_eachcluster=100}
      setting_mat=matrix( c(a,k,MNAR0lambda_star), nrow=30, ncol=3, byrow = TRUE) 
      colnames(setting_mat) <- c("a","k","lambda_star")
      res_mat=matrix(0,30,3+3+3+length(tuninglambdalist))
      colnames(res_mat) <- c("cer_kmeans","cer_kpod","cer_pkpod",
                             "mse_kmeans","mse_kpod","mse_pkpod",
                             "n", "miss%", "selectedlambda",
                             paste0("lambda",tuninglambdalist) )
      for (rept in 1:30){
        res_once <- try( experiment_MNAR0_p2(a=a,k=k,n=samplesize_eachcluster*k, 
                                             mustar=mustar, MNAR0lambda_star=MNAR0lambda_star, 
                                             tuninglambdalist=tuninglambdalist)  )
        if (inherits(res_once,"try-error")){next}
        res_mat[rept,] <- c( res_once$cer_km, res_once$cer_kpod, res_once$cer_pkpod,
                             res_once$mse_km, res_once$mse_kpod, res_once$mse_pkpod,
                             res_once$missdata_n,res_once$missdata_proportion,res_once$selectedlambda,
                             res_once$instability_values )
      }
      ExperimentResults_MNAR0_p2$allresults[[setting]] <- cbind(setting_mat, res_mat)
      setting <- setting + 1
      
      print(c( "a=",a ,"k=", k, "lambda_star=",MNAR0lambda_star ))
    }
  }
}


### summarize results -------------------

ExperimentResults_MNAR0_p2$summarizedresults$cer=as.data.frame(do.call(rbind, ExperimentResults_MNAR0_p2$allresults)) %>%
  select(a, k, lambda_star, cer_kmeans, cer_kpod, cer_pkpod) %>%
  pivot_longer(
    cols = c(cer_kmeans, cer_kpod, cer_pkpod),
    names_to = "method",
    values_to = "cer_value"
  ) %>%
  group_by(a, k, lambda_star, method) %>%
  summarise(
    cer_mean = mean(cer_value),
    cer_sd   = sd(cer_value),
    .groups = "drop"
  )

ExperimentResults_MNAR0_p2$summarizedresults$mse=as.data.frame(do.call(rbind, ExperimentResults_MNAR0_p2$allresults)) %>%
  select(a, k, lambda_star, mse_kmeans, mse_kpod, mse_pkpod) %>%
  pivot_longer(
    cols = c(mse_kmeans, mse_kpod, mse_pkpod),
    names_to = "method",
    values_to = "mse_value"
  ) %>%
  group_by(a, k, lambda_star, method) %>%
  summarise(
    mse_mean = mean(mse_value),
    mse_sd   = sd(mse_value),
    .groups = "drop"
  )

### draw figures --------------------

ExperimentResults_MNAR0_p2$figures$cer <- ggplot(ExperimentResults_MNAR0_p2$summarizedresults$cer,
                                                      aes(x = log(lambda_star), y = cer_mean,
                                                          color = method, linetype = method, shape = method, group = method)) + 
  scale_color_manual( values = c("cer_kmeans" = "black", "cer_kpod"="blue", "cer_pkpod"="red"),
                      labels = c("k-means", "kPOD", "proposed")) +
  scale_linetype_manual( values = c("cer_kmeans" = "dotted", "cer_kpod"="dashed", "cer_pkpod"="solid"),
                         labels = c("k-means", "kPOD", "proposed")) +
  scale_shape_manual( values = c("cer_kmeans" = 15, "cer_kpod"=17, "cer_pkpod"=16),
                      labels = c("k-means", "kPOD", "proposed")) +
  geom_line() + geom_point(size=1) +
  geom_errorbar( aes(ymin = cer_mean - cer_sd, ymax = cer_mean + cer_sd), width = 0.15 , linetype="solid" ) +
  facet_grid(k ~ a, 
             labeller = labeller( a = function(z) paste0("a = ", z), k = function(z) paste0("k = ", z))
  ) +
  theme_bw() + theme( legend.position = "bottom",legend.margin = margin(t = -8)) + 
  labs(x = expression(log(lambda^{"* "})),y = "CER", color = NULL,linetype = NULL,shape = NULL)

ExperimentResults_MNAR0_p2$figures$mse <- ggplot(ExperimentResults_MNAR0_p2$summarizedresults$mse,
                                                      aes(x = log(lambda_star), y = mse_mean,
                                                          color = method, linetype = method, shape = method, group = method)) + 
  scale_color_manual( values = c("mse_kmeans" = "black", "mse_kpod"="blue", "mse_pkpod"="red"),
                      labels = c("k-means", "kPOD", "proposed")) +
  scale_linetype_manual( values = c("mse_kmeans" = "dotted", "mse_kpod"="dashed", "mse_pkpod"="solid"),
                         labels = c("k-means", "kPOD", "proposed")) +
  scale_shape_manual( values = c("mse_kmeans" = 15, "mse_kpod"=17, "mse_pkpod"=16),
                      labels = c("k-means", "kPOD", "proposed")) +
  geom_line() + geom_point(size=1) +
  geom_errorbar( aes(ymin = mse_mean - mse_sd, ymax = mse_mean + mse_sd), width = 0.15 , linetype="solid" ) +
  facet_grid(k ~ a, scales = "free_y",
             labeller = labeller( a = function(z) paste0("a = ", z), k = function(z) paste0("k = ", z))
  ) +
  theme_bw() + theme( legend.position = "bottom",legend.margin = margin(t = -8)) + 
  labs(x = expression(log(lambda^{"* "})),y = "MSE", color = NULL,linetype = NULL,shape = NULL)

ExperimentResults_MNAR0_p2$figures$cer
ExperimentResults_MNAR0_p2$figures$mse


save(ExperimentResults_MNAR0_p2, file="ExperimentResults_MNAR0_p2.Rdata")

