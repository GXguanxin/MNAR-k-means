##############################
########   Table 1  ##########
##############################

##### loss convergence -----------

experiment_loss_convergence <- function(data,k=3,
                                        mustar=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                        tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])  ){
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
  
  pkpodres <- pkpod_1218(x=data$Missing,k=k,lambda=selectedlambda,nstart=100) #大约0.5s
  cer_pkpod <- CER(pkpodres$cluster,data$Misslabel)
  mse_pkpod <- diffcode(mustar,pkpodres$centers)
  loss_value_vec <- pkpodres$obj_list/nrow(data$Missing)
  
  return(list( selectedlambda=selectedlambda, instability_values=instability_values, instability_values_sd=instability_values_sd,
               cer_pkpod=cer_pkpod, mse_pkpod=mse_pkpod, 
               loss_value_vec=loss_value_vec
  ))
}

ExperimentResults_loss_convergence <- list( 
  allsettings = list( missmechanism = c("MNAR0","MNAR1","MNAR2"), missproportion = c(0.1,0.3,0.5),
                      k=3, p=50, tuninglambdalist=c(0.1, seq(0,20,by=2)[-1]) ),
  allresults = matrix(0,9,30)
)
setting=1
for (missmechanism in c("MNAR0","MNAR1","MNAR2")){
  for (missproportion in c(0.1,0.3,0.5)){
    for (rept in 1:30){
      set.seed(2026+rept*10)
      data <- generate_data_different_MNAR_k3_p50(missmechanism=missmechanism, missproportion=missproportion)
      res_once <- experiment_loss_convergence( data,k=3,
                                               mustar=rbind( c( 0,0 ,rep(0,48)) , c( 3,0 ,rep(0,48)), c( 1.5, sqrt(6.75), rep(0,48) ) ), 
                                               tuninglambdalist=c(0.1, seq(0,20,by=2)[-1])  )
      ExperimentResults_loss_convergence$allresults[setting,rept] <- length(res_once$loss_value_vec)
      if(rept%%10==0){print(rept)}
    }
    setting <- setting + 1
    print(c(missmechanism, missproportion))
  }
}


print(apply(ExperimentResults_loss_convergence_1228$allresults,1,mean))
print(apply(ExperimentResults_loss_convergence$allresults,1,sd))


