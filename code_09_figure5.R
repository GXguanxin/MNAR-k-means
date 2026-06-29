###############################
########   Figure 5  ##########
###############################

# a=3, k=3, p=2, MNAR0 lambdastar=1

experiment_tuningpara_effect <- function(data,k=3,
                                         mustar=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                         tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])  ){
  instability_values <- rep(0,length(tuninglambdalist))
  instability_values_sd <- rep(0,length(tuninglambdalist))
  cers_vec <- rep(0,length(tuninglambdalist))
  mses_vec <- rep(0,length(tuninglambdalist))
  activefeatures_vec <- rep(0,length(tuninglambdalist))
  for (t in 1:length(tuninglambdalist)){
    instability_res <- instability_1218(x=data$Missing,k=k,tuningpara = tuninglambdalist[t], 
                                        clusteringfun = pkpod_1218, nstart=10, perm.max = 10  ) #大约2s
    instability_values[t] <- instability_res$Svalue
    instability_values_sd[t] <- instability_res$Svalue_sd
    
    pkpodres <- pkpod_1218(x=data$Missing,k=k,lambda=tuninglambdalist[t],nstart=100)
    cers_vec[t] <- CER(pkpodres$cluster,data$Misslabel)
    mses_vec[t] <- diffcode(mustar,pkpodres$centers)
    activefeatures_vec[t] <- sum( apply( pkpodres$centers,2,function(z) sum(z^2) ) > 1e-3 )
  }
  selectedlambda <- tuninglambdalist[which.min(instability_values)]
  if ( selectedlambda==tail(tuninglambdalist,1) ){
    lowerbound <- tail(instability_values,1) + tail(instability_values_sd,1)
    selectedlambda <- tuninglambdalist[min(which( instability_values < lowerbound ))]
  }
  
  return(list( instability_values=instability_values, instability_values_sd=instability_values_sd,
               cers_vec=cers_vec, mses_vec=mses_vec, activefeatures_vec=activefeatures_vec,
               selectedlambda=selectedlambda) )
  
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

tuninglambdalist=c(0.01, seq(0,2,by=0.2)[-1], 4, 6, 8, 10)
ExperimentResults_tuningpara_effect <- list(
  setting = list(missmechanism="MNAR0", MNAR0lambdastar=1, 
                 tuninglambdalist=c(0.01, seq(0,2,by=0.2)[-1], 4, 6, 8, 10)  ),
  results = list( instability_values = matrix(0,30,length(tuninglambdalist)),
                  cers = matrix(0,30,length(tuninglambdalist)),
                  mses = matrix(0,30,length(tuninglambdalist)),
                  activefeatures = matrix(0,30,length(tuninglambdalist)),
                  selectedlambdas = rep(0,30)
  ) )
for (rept in 1:30){
  set.seed(2026+rept)
  data <- generate_data(n=100*3, p=2, k=3, 
                                 centers=rbind( c( 0,0 ,rep(0,0)) , c( 3,0 ,rep(0,0)), c( 1.5, sqrt(6.75), rep(0,0) ) ), 
                                 variance_vec = c(rep(1,2),rep(2,0)), 
                                 missmechanism="MNAR0", MNAR0lambda_star = 1, nomissing=FALSE )
  res_once <- experiment_tuningpara_effect(data=data,k=3,
                                           mustar=rbind( c( 0,0 ,rep(0,0)) , c( 3,0 ,rep(0,0)), c( 1.5, sqrt(6.75), rep(0,0) ) ), 
                                           tuninglambdalist=tuninglambdalist)
  
  ExperimentResults_tuningpara_effect$results$instability_values[rept,] <- res_once$instability_values
  ExperimentResults_tuningpara_effect$results$cers[rept,] <- res_once$cers_vec
  ExperimentResults_tuningpara_effect$results$mses[rept,] <- res_once$mses_vec
  ExperimentResults_tuningpara_effect$results$activefeatures[rept,] <- res_once$activefeatures_vec
  ExperimentResults_tuningpara_effect$results$selectedlambdas[rept] <- res_once$selectedlambda
  
  if(rept%%10==0){print(rept)}
}


# summarize data
ExperimentResults_tuningpara_effect$summarizedresults <- cbind(
  "lambda"=tuninglambdalist, 
  "instability_mean"=apply(ExperimentResults_tuningpara_effect$results$instability_values, 2, mean),
  "instability_sd"=apply(ExperimentResults_tuningpara_effect$results$instability_values, 2, sd),
  "cer_mean"=apply(ExperimentResults_tuningpara_effect$results$cers, 2, mean),
  "cer_sd"=apply(ExperimentResults_tuningpara_effect$results$cers, 2, sd),
  "mse_mean"=apply(ExperimentResults_tuningpara_effect$results$mses, 2, mean),
  "mse_sd"=apply(ExperimentResults_tuningpara_effect$results$mses, 2, sd)
)


# draw figures
ExperimentResults_tuningpara_effect$figures$instabilityVScer <- ggplot(as.data.frame(ExperimentResults_tuningpara_effect$summarizedresults)[2:15,], aes(x = log(lambda,10)) ) + 
  geom_line(aes(y = instability_mean, color="Instability", linetype="Instability") ) +
  geom_point(aes(y = instability_mean, shape="Instability", color = "Instability"), size=2) +
  geom_errorbar( aes(ymin = instability_mean - instability_sd, ymax = instability_mean + instability_sd,
                     linetype="Instability", color="Instability"), 
                 width = 0.05 ) +
  geom_line(aes(y = cer_mean, linetype="CER", color = "CER") ) +
  geom_point(aes(y = cer_mean, shape="CER", color = "CER"), size=2 ) +
  geom_errorbar( aes(ymin = cer_mean - cer_sd, ymax = cer_mean + cer_sd,
                     linetype="CER", color = "CER"), 
                 width = 0.05) +
  scale_y_continuous(limits = c(0, 0.3), name = "Instability", sec.axis = sec_axis( ~ ., name = "CER")  ) + 
  scale_color_manual( values = c("Instability" = "red", "CER" = "black") ) +
  scale_linetype_manual( values = c("Instability" = "dashed", "CER" = "solid") ) +
  scale_shape_manual( values = c("Instability" = 17, "CER" = 16) ) +
  theme_bw() +
  theme( legend.position = c(0.8,0.85),
         axis.title.y.left = element_text(),
         axis.title.y.right  = element_text() ) + 
  labs(x = expression(log(lambda)),color = NULL, linetype = NULL, shape = NULL) +
  geom_vline( xintercept = 0, linetype = "dotted", size=1, color = "grey60") +
  annotate("text", x = 0, y = 0.03,
           label = expression(lambda == 1), color = "grey60", hjust = -0.1 )

ExperimentResults_tuningpara_effect$figures$instabilityVScer

