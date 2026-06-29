#############################
########   Table 2 ##########
#############################


###### initialization fluctuation -------------------------

ExperimentResults_initialization_fluctuation <- list( allresults = list() )

#MNAR0 0.1 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR0", missproportion=0.1)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=12,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[1]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))


#MNAR0 0.3 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR0", missproportion=0.3)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=18,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[2]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))


#MNAR0 0.5 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR0", missproportion=0.5)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=16,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[3]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))




#MNAR1 0.1 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR1", missproportion=0.1)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=18,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[4]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))



#MNAR1 0.3 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR1", missproportion=0.3)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=16,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[5]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))



#MNAR1 0.5 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR1", missproportion=0.5)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=4,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[6]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))


#MNAR2 0.1 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR2", missproportion=0.1)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=18,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[7]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))


#MNAR2 0.3 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR2", missproportion=0.3)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=14,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[8]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))



#MNAR2 0.5 ---
obj_nstart_summary <- matrix(0,6,30)
obj_nstart_sd <- rep(0,30)
obj_decresing <- rep(0,30)
for (rept in 1:30){
  set.seed(2026+rept*10)
  data <- generate_data_different_MNAR_k3_p50(missmechanism="MNAR2", missproportion=0.5)
  pkpodres <- pkpod_1218(x=data$Missing,k=3,lambda=6,nstart=100) 
  obj_nstart_summary[,rept] <- as.numeric(summary(pkpodres$obj_nstart)/300)
  obj_nstart_sd[rept] <- sd(pkpodres$obj_nstart/300)
  obj_decresing[rept] <- (pkpodres$obj_list[1] - tail(pkpodres$obj_list,1))/300
}
ExperimentResults_initialization_fluctuation$allresults[[9]] <- rbind(
  obj_nstart_summary, obj_nstart_sd, obj_decresing )

print(mean(obj_nstart_sd))
print(mean(apply(obj_nstart_summary,2,function(z) z[6]-z[1] )))
print(mean(apply(obj_nstart_summary,2,function(z) z[5]-z[2] )))
print(mean(obj_decresing))

