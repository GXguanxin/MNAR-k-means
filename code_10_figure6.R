##############################
########  Figure 6  ##########
##############################



### run experiments for figure 6-1 ----------

ExperimentResults_consistency_to_truth <- list(
  settings = list(
    samplesize_each = c(0.3,0.5,1,3,5,7,10,30,50,100)*100, 
    p=2, k=3, 
    sigma2_list = (  (1/log( c(0.3,0.5,1,3,5,7,10,30,50,100)*100*3 ,base = exp(1) ))^2  )*30, 
    mustar=rbind( c(-sqrt(6)/2, -sqrt(6)/2 ), 
                  c( (sqrt(6)+3*sqrt(2))/4, (sqrt(6)-3*sqrt(2))/4 ),
                  c( (sqrt(6)-3*sqrt(2))/4, (sqrt(6)+3*sqrt(2))/4 )  ),
    missmechanism="MNAR0", MNAR0lambda_star = 1,
    lambda_opt_list = 1-1/( 2*1*(  (1/log( c(0.3,0.5,1,3,5,7,10,30,50,100)*100*3 ,base = exp(1) ))^2  )*30 + 1 )
  ),
  allresults=list(mse=list(), cer=list())
)


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


mustar=rbind( c(-sqrt(6)/2, -sqrt(6)/2 ), 
              c( (sqrt(6)+3*sqrt(2))/4, (sqrt(6)-3*sqrt(2))/4 ),
              c( (sqrt(6)-3*sqrt(2))/4, (sqrt(6)+3*sqrt(2))/4 )  )

setting=1
for (samplesize_each in c(0.3,0.5,1,3,5,7,10,30,50,100)){
  n <- samplesize_each*100*3
  sigma2 <- ((1/log(n,base = exp(1) ))^2)*30
  lambda_opt <- 1-1/( 2*1*sigma2 +1 )
  res_mat <- matrix(0,30,6)
  for (rept in 1:30){
    set.seed(2026+samplesize_each*rept)
    data <- generate_data(n=n, p=2, k=3, 
                                   centers=mustar,
                                   variance_vec = c(rep(sigma2,2)), 
                                   missmechanism="MNAR0", MNAR0lambda_star = 1, nomissing=FALSE )
    kmres <- kmeans(data$Orig,centers=3,iter.max = 100,nstart=100)
    res_mat[rept,1] <- diffcode(mustar,kmres$centers)
    res_mat[rept,4] <- CER(kmres$cluster,data$Origlabel)
    kpodres <- pkpod_1218(x=data$Missing,k=3,lambda=0,nstart=100)
    res_mat[rept,2] <- diffcode(mustar,kpodres$centers)
    res_mat[rept,5] <- CER(kpodres$cluster,data$Misslabel)
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt, nstart=100)
    res_mat[rept,3] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,6] <- CER(pkpodres$cluster,data$Misslabel)
    
    if (rept%%10==0){print(rept)}
  }
  ExperimentResults_consistency_to_truth$allresults$mse[[setting]] <- res_mat[,1:3]
  ExperimentResults_consistency_to_truth$allresults$cer[[setting]] <- res_mat[,4:6]
  setting = setting +1
  
  print( c(samplesize_each,sigma2,lambda_opt,
           nrow(data$Missing), sum(is.na(data$Missing))/(nrow(data$Missing)*2),
           diffcode(mustar,kmres$centers),diffcode(mustar,kpodres$centers), diffcode(mustar,pkpodres$centers)) )
}



# summarize data
tmp <- do.call(rbind, ExperimentResults_consistency_to_truth$allresults$mse)
tmp <- cbind( rep( c(0.3,0.5,1,3,5,7,10,30,50,100)*100*3, each=30 ) ,tmp)
colnames(tmp) <- c("n","kmeans","kpod","pkpod")

ExperimentResults_consistency_to_truth$summarizedresults$mse <- as.data.frame(tmp) %>%
  pivot_longer(
    cols = c(kmeans,kpod,pkpod),
    names_to = "method",
    values_to = "mse_value"
  ) %>%
  group_by(n,method) %>%
  summarise(
    mse_mean = mean(mse_value),
    mse_sd   = sd(mse_value),
    .groups = "drop"
  )


# draw figures
ExperimentResults_consistency_to_truth$figures$mse <- ggplot(ExperimentResults_consistency_to_truth$summarizedresults$mse,
                                                                  aes(x = n, y = mse_mean,
                                                                      color = method, linetype = method, shape = method, group = method)) + 
  scale_color_manual( values = c("kmeans" = "black", "kpod"="blue", "pkpod"="red"),
                      labels = c("k-means", "kPOD", "proposed")) +
  scale_linetype_manual( values = c("kmeans" = "dotted", "kpod"="dashed", "pkpod"="solid"),
                         labels = c("k-means", "kPOD", "proposed")) +
  scale_shape_manual( values = c("kmeans" = 15, "kpod"=17, "pkpod"=16),
                      labels = c("k-means", "kPOD", "proposed")) +
  geom_line() + geom_point(size=1.5) +
  geom_errorbar( aes(ymin = mse_mean - mse_sd, ymax = mse_mean + mse_sd), width = 0.05 , linetype="solid" ) +
  theme_bw() + theme( legend.position = c(0.95, 0.95),legend.justification = c(1, 1)) + 
  coord_cartesian(ylim = c(0, 4)) + scale_x_log10()+
  labs(x = "n",y = "MSE", color = NULL,linetype = NULL,shape = NULL)


ExperimentResults_consistency_to_truth$figures$mse



### run experiments for figure 6-2 ----------

mustar=rbind( c(-sqrt(6)/2, -sqrt(6)/2 ), 
              c( (sqrt(6)+3*sqrt(2))/4, (sqrt(6)-3*sqrt(2))/4 ),
              c( (sqrt(6)-3*sqrt(2))/4, (sqrt(6)+3*sqrt(2))/4 )  )

ExperimentResults_consistency_to_truth$allresults_otherlambda <- list(mse=list(), cer=list())

setting=1
for (samplesize_each in c(0.3,0.5,1,3,5,7,10,30,50,100)){
  n <- samplesize_each*100*3
  sigma2 <- ((1/log(n,base = exp(1) ))^2)*30
  lambda_opt <- 1-1/( 2*1*sigma2 +1 )
  res_mat <- matrix(0,30,18)
  for (rept in 1:30){
    set.seed(2026+samplesize_each*rept)
    data <- generate_data(n=n, p=2, k=3, 
                                   centers=mustar,
                                   variance_vec = c(rep(sigma2,2)), 
                                   missmechanism="MNAR0", MNAR0lambda_star = 1, nomissing=FALSE )
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt, nstart=100)
    res_mat[rept,1] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,10] <- CER(pkpodres$cluster,data$Misslabel)
    
    #smaller than lambda_opt
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*0.1, nstart=100)
    res_mat[rept,2] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,11] <- CER(pkpodres$cluster,data$Misslabel)
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*0.3, nstart=100)
    res_mat[rept,3] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,12] <- CER(pkpodres$cluster,data$Misslabel)
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*0.5, nstart=100)
    res_mat[rept,4] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,13] <- CER(pkpodres$cluster,data$Misslabel)
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*0.7, nstart=100)
    res_mat[rept,5] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,14] <- CER(pkpodres$cluster,data$Misslabel)
    
    #larger than lambda_opt
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*1.5, nstart=100)
    res_mat[rept,6] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,15] <- CER(pkpodres$cluster,data$Misslabel)
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*2, nstart=100)
    res_mat[rept,7] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,16] <- CER(pkpodres$cluster,data$Misslabel)
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*2.5, nstart=100)
    res_mat[rept,8] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,17] <- CER(pkpodres$cluster,data$Misslabel)
    pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=lambda_opt*3, nstart=100)
    res_mat[rept,9] <- diffcode(mustar,pkpodres$centers)
    res_mat[rept,18] <- CER(pkpodres$cluster,data$Misslabel)
    
    
    if (rept%%10==0){print(rept)}
  }
  ExperimentResults_consistency_to_truth$allresults_otherlambda$mse[[setting]] <- res_mat[,1:9]
  ExperimentResults_consistency_to_truth$allresults_otherlambda$cer[[setting]] <- res_mat[,10:18]
  setting = setting +1
  
  print(setting)
}


# summarize data
tmp <- do.call(rbind, ExperimentResults_consistency_to_truth$allresults_otherlambda$mse)
tmp <- cbind( rep( c(0.3,0.5,1,3,5,7,10,30,50,100)*100*3, each=30 ) ,tmp)
colnames(tmp) <- c("n","opt","opt01","opt03","opt05","opt07",
                   "opt15","opt2","opt25","opt3")

ExperimentResults_consistency_to_truth$summarizedresults_otherlambda$mse <- as.data.frame(tmp) %>%
  pivot_longer(
    cols = c(opt,opt01,opt03,opt05,opt07,opt15,opt2,opt25,opt3),
    names_to = "lambda_value",
    values_to = "mse_value"
  ) %>%
  group_by(n,lambda_value) %>%
  summarise(
    mse_mean = mean(mse_value),
    mse_sd   = sd(mse_value),
    .groups = "drop"
  )


# draw figures
library(latex2exp)

ExperimentResults_consistency_to_truth$figures_otherlambda$mse <- ggplot(ExperimentResults_consistency_to_truth$summarizedresults_otherlambda$mse,
                                                                              aes(x = n, y = mse_mean,
                                                                                  color = lambda_value, linetype = lambda_value, shape = lambda_value, group = lambda_value)) + 
  scale_color_manual( values = c( "opt"   = "red",
                                  "opt01" = "#1f77b4", "opt03" = "#4f94cd", "opt05" = "#63b8ff", "opt07" = "#00bfff", 
                                  "opt15" = "#ff7f0e", "opt2"  = "#ffa500", "opt25" = "#ff8c00", "opt3"  = "#d2691e"  ),
                      labels = c( "opt"   = TeX("$\\lambda_{opt}$"),
                                  "opt01" = TeX("$0.1\\lambda_{opt}$"), "opt03" = TeX("$0.3\\lambda_{opt}$"), "opt05" = TeX("$0.5\\lambda_{opt}$"), "opt07" = TeX("$0.7\\lambda_{opt}$"),
                                  "opt15" = TeX("$1.5\\lambda_{opt}$"), "opt2"  = TeX("$2\\lambda_{opt}$"), "opt25" = TeX("$2.5\\lambda_{opt}$"), "opt3"  = TeX("$3\\lambda_{opt}$")     )
  ) +
  scale_linetype_manual( values = c( "opt"   = "solid",
                                     "opt01" = "dashed", "opt03" = "dashed", "opt05" = "dashed", "opt07" = "dashed",
                                     "opt15" = "dotted", "opt2"  = "dotted", "opt25" = "dotted", "opt3"  = "dotted" ),
                         labels = c( "opt"   = TeX("$\\lambda_{opt}$"),
                                     "opt01" = TeX("$0.1\\lambda_{opt}$"), "opt03" = TeX("$0.3\\lambda_{opt}$"), "opt05" = TeX("$0.5\\lambda_{opt}$"), "opt07" = TeX("$0.7\\lambda_{opt}$"),
                                     "opt15" = TeX("$1.5\\lambda_{opt}$"), "opt2"  = TeX("$2\\lambda_{opt}$"), "opt25" = TeX("$2.5\\lambda_{opt}$"), "opt3"  = TeX("$3\\lambda_{opt}$")     )
  ) +
  scale_shape_manual( values = c( "opt"   = 16, 
                                  "opt01" = 17, "opt03" = 15, "opt05" = 3,  "opt07" = 7,
                                  "opt15" = 8,  "opt2"  = 4,  "opt25" = 9,  "opt3"  = 10  ),
                      labels = c( "opt"   = TeX("$\\lambda_{opt}$"),
                                  "opt01" = TeX("$0.1\\lambda_{opt}$"), "opt03" = TeX("$0.3\\lambda_{opt}$"), "opt05" = TeX("$0.5\\lambda_{opt}$"), "opt07" = TeX("$0.7\\lambda_{opt}$"),
                                  "opt15" = TeX("$1.5\\lambda_{opt}$"), "opt2"  = TeX("$2\\lambda_{opt}$"), "opt25" = TeX("$2.5\\lambda_{opt}$"), "opt3"  = TeX("$3\\lambda_{opt}$")     )
  ) +
  geom_line() + geom_point(size=1.5) +
  geom_errorbar( aes(ymin = mse_mean - mse_sd, ymax = mse_mean + mse_sd), width = 0.05 , linetype="solid" ) +
  theme_bw() + theme( legend.position = c(0.99, 0.99),legend.justification = c(1, 1),
                      legend.text = element_text(size = 8),legend.spacing.y = unit(0, "cm") ) + 
  scale_y_continuous(limits = c(0, 3)) + scale_x_log10()+
  labs(x = "n",y = "MSE", color = NULL,linetype = NULL,shape = NULL)


ExperimentResults_consistency_to_truth$figures_otherlambda$mse

