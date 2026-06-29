########################################
########   Figure 3 (k=3)  #############
########################################

generate_data_different_MNAR_k3_p50 <- function(missmechanism,missproportion){
  if (missmechanism == "MNAR0"){
    if (missproportion == 0.1){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR0",MNAR0lambda_star=50, nomissing=FALSE )
    }
    if (missproportion == 0.3){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR0",MNAR0lambda_star=5, nomissing=FALSE )
    }
    if (missproportion == 0.5){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR0",MNAR0lambda_star=1.5, nomissing=FALSE )
    }
  }
  
  if (missmechanism == "MNAR1"){
    if (missproportion == 0.1){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR1",MNAR1phi = c(-15,0), nomissing=FALSE )
    }
    if (missproportion == 0.3){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR1",MNAR1phi = c(-1.5,0), nomissing=FALSE )
    }
    if (missproportion == 0.5){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR1",MNAR1phi = c(-0.1,0), nomissing=FALSE )
    }
  }
  
  if (missmechanism == "MNAR2"){
    if (missproportion == 0.1){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR2",MNAR2percent = rep(0.1,50), nomissing=FALSE )
    }
    if (missproportion == 0.3){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR2",MNAR2percent = rep(0.3,50), nomissing=FALSE )
    }
    if (missproportion == 0.5){
      data <- generate_data_20251216(n=100*3, p=50, k=3, 
                                     centers=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                     variance_vec = c(rep(1,2),rep(1,48)), 
                                     missmechanism="MNAR2",MNAR2percent = rep(0.5,50), nomissing=FALSE )
    }
  }
  
  return(data)
}

experiment_different_MNAR_k3_p50 <- function(data,k=3,
                                             mustar=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                             tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])  ){
  missdata_n <- dim(data$Missing)[1]
  missdata_proportion <- sum(is.na(data$Missing))/(nrow(data$Missing)*ncol(data$Missing))
  
  kmres <- kmeans(data$Orig,centers=k,iter.max = 100,nstart=100)
  cer_km <- CER(kmres$cluster,data$Origlabel)
  mse_km <- diffcode(mustar,kmres$centers)
  
  data$meanimpt <- data$Missing
  data$meanimpt[is.na(data$Missing)] <- matrix(colMeans(data$Missing,na.rm=TRUE),missdata_n,50,byrow = TRUE)[is.na(data$Missing)]
  meanimptres <- kmeans(data$meanimpt,centers=k,iter.max = 100,nstart=100)
  cer_meanimpt <- CER(meanimptres$cluster,data$Misslabel)
  mse_meanimpt <- diffcode(mustar,meanimptres$centers)
  
  miceres <- fun_mice(x=data$Missing, k=k)  
  cer_mice <- CER(miceres$cluster,data$Misslabel)
  mse_mice <- diffcode(mustar,miceres$centers)
  
  data$knnimpt <- impute.knn(as.matrix(data$Missing), k = 10, rng.seed = 2026*88)$data  #10-nearest-neighbors
  knnimptres <- kmeans(data$knnimpt,centers=k,iter.max = 100,nstart=100)
  cer_knnimpt <- CER(knnimptres$cluster,data$Misslabel)
  mse_knnimpt <- diffcode(mustar,knnimptres$centers)
  
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
               cer_km=cer_km, cer_meanimpt=cer_meanimpt,cer_mice=cer_mice,cer_knnimpt=cer_knnimpt, cer_kpod=cer_kpod, cer_pkpod=cer_pkpod, 
               mse_km=mse_km, mse_meanimpt=mse_meanimpt, mse_mice=mse_mice, mse_knnimpt=mse_knnimpt, mse_kpod=mse_kpod, mse_pkpod=mse_pkpod,
               instability_values=instability_values, instability_values_sd=instability_values_sd
  ))
}


### experiments ----------------

ExperimentResults_different_MNAR_k3_p50=list(
  allsettings=list(
    missmechanism_list=c("MNAR0","MNAR1","MNAR2"),
    MNAR0lambda_star_list=c(50, 5, 1.5),
    MNAR1phi_list=rbind( c(-15,0), c(-1.5,0), c(-0.1,0) ),
    MNAR2percent_list=rbind( rep(0.1,50), rep(0.3,50), rep(0.5,50) ),
    tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])
  ),
  allresults=list()
)


tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])
setting=1
for (missmechanism in c("MNAR0","MNAR1","MNAR2")){
  for (missproportion in c(0.1,0.3,0.5)){
    setting_mat=matrix( c(missmechanism,missproportion), nrow=30, ncol=2, byrow = TRUE) 
    colnames(setting_mat) <- c("mechanism","missproportion")
    res_mat=matrix(0,30,6+6+2+length(tuninglambdalist))
    colnames(res_mat) <- c("cer_kmeans","cer_meanimpt","cer_mice","cer_knnimpt","cer_kpod","cer_pkpod",
                           "mse_kmeans","mse_meanimpt","mse_mice","mse_knnimpt","mse_kpod","mse_pkpod",
                           "miss%", "selectedlambda",
                           paste0("lambda",tuninglambdalist) )
    for (rept in 1:30){
      set.seed(2026+rept*10)
      data <- generate_data_different_MNAR_k3_p50(missmechanism=missmechanism, missproportion=missproportion)
      res_once <- try( experiment_different_MNAR_k3_p50(data,k=3,
                                                        mustar=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                                        tuninglambdalist=tuninglambdalist ) )  
      if (inherits(res_once,"try-error")){next}
      res_mat[rept,] <- c( res_once$cer_km, res_once$cer_meanimpt, res_once$cer_mice, res_once$cer_knnimpt, res_once$cer_kpod, res_once$cer_pkpod,
                           res_once$mse_km, res_once$mse_meanimpt, res_once$mse_mice, res_once$mse_knnimpt, res_once$mse_kpod, res_once$mse_pkpod,
                           res_once$missdata_proportion, res_once$selectedlambda,
                           res_once$instability_values )
      if(rept%%10==0){print(rept)}
    }
    ExperimentResults_different_MNAR_k3_p50$allresults[[setting]] <- cbind(setting_mat, res_mat)
    setting <- setting + 1
    
    print(c(missmechanism, missproportion))
  }
}


### summarize results ----------------

ExperimentResults_different_MNAR_k3_p50$summarizedresults$cer <- as.data.frame(do.call(rbind, ExperimentResults_different_MNAR_k3_p50$allresults)) %>%
  select(mechanism, missproportion, cer_kmeans, cer_meanimpt, cer_mice, cer_knnimpt, cer_kpod, cer_pkpod) %>%
  pivot_longer(
    cols = c(cer_kmeans, cer_meanimpt, cer_mice, cer_knnimpt, cer_kpod, cer_pkpod),
    names_to = "method",
    values_to = "cer_value"
  ) %>%
  mutate( cer_value = as.numeric(cer_value) ) %>%
  mutate( 
    method = factor( method,
                     levels = c("cer_kmeans","cer_meanimpt", "cer_mice","cer_knnimpt", "cer_kpod", "cer_pkpod"),
                     labels = c("k-means","meanimpt", "MICE", "knnimpt" , "k-POD", "proposed")  )
  )

ExperimentResults_different_MNAR_k3_p50$summarizedresults$mse <- as.data.frame(do.call(rbind, ExperimentResults_different_MNAR_k3_p50$allresults)) %>%
  select(mechanism, missproportion, mse_kmeans,mse_meanimpt, mse_mice, mse_knnimpt, mse_kpod, mse_pkpod) %>%
  pivot_longer(
    cols = c(mse_kmeans,mse_meanimpt, mse_mice, mse_knnimpt, mse_kpod, mse_pkpod),
    names_to = "method",
    values_to = "mse_value"
  ) %>%
  mutate( mse_value = as.numeric(mse_value) ) %>%
  mutate( 
    method = factor( method,
                     levels = c("mse_kmeans","mse_meanimpt", "mse_mice","mse_knnimpt", "mse_kpod", "mse_pkpod"),
                     labels = c("k-means","meanimpt", "MICE","knnimpt", "k-POD", "proposed")  )
  )


### draw figures ----------------

ExperimentResults_different_MNAR_k3_p50$figures$cer <- ggplot(ExperimentResults_different_MNAR_k3_p50$summarizedresults$cer,
                                                                   aes(x = method, y = cer_value, color = method)) + 
  geom_boxplot(outlier.size=0.8) + 
  scale_color_manual( values = c( "k-means"="black","meanimpt"="purple" , "MICE"="green4", "knnimpt"="orange" , "k-POD"="blue","proposed"="red") ) + 
  facet_grid(missproportion ~ mechanism, scales = "free_y",
             labeller = labeller( missproportion = function(z) paste0(as.numeric(z) * 100, "%") )
  ) + 
  labs(x="Method", y="MSE") + theme_bw() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none" )+
  scale_y_continuous( breaks = scales::pretty_breaks(n = 4) ) +
  scale_x_discrete(limits = c("k-means","meanimpt", "MICE", "knnimpt", "k-POD", "proposed"),
                   labels = c("k-means", "mean" ,"MICE", "KNN", "k-POD", "proposed") )



ExperimentResults_different_MNAR_k3_p50$figures$mse <- ggplot(ExperimentResults_different_MNAR_k3_p50$summarizedresults$mse,
                                                                   aes(x = method, y = mse_value, color = method)) + 
  geom_boxplot(outlier.size=0.8) + 
  scale_color_manual( values = c( "k-means"="black","meanimpt"="purple" , "MICE"="green4", "knnimpt"="orange" , "k-POD"="blue","proposed"="red") ) + 
  facet_grid(missproportion ~ mechanism, scales = "free_y",
             labeller = labeller( missproportion = function(z) paste0(as.numeric(z) * 100, "%") )
  ) + 
  labs(x="Method", y="MSE") + theme_bw() +
  theme( axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none" )+
  scale_y_continuous( breaks = scales::pretty_breaks(n = 4) ) +
  scale_x_discrete(limits = c("k-means","meanimpt", "MICE", "knnimpt", "k-POD", "proposed"),
                   labels = c("k-means", "mean" ,"MICE", "KNN", "k-POD", "proposed") )


ExperimentResults_different_MNAR_k3_p50$figures$cer
ExperimentResults_different_MNAR_k3_p50$figures$mse
