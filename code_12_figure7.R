#################################
########   Figure 7  ############
#################################

experiment_realdata <- function(data,k=4, tuninglambdalist=c(0.001, 0.01, 0.1, 1, 10)  ){
  missdata_n <- dim(data$Missing)[1]
  missdata_proportion <- sum(is.na(data$Missing))/(nrow(data$Missing)*ncol(data$Missing))
  
  data$zeroimpt <- data$Missing
  data$zeroimpt[is.na(data$Missing)] <- 0
  kmres <- kmeans(data$zeroimpt,centers=k,iter.max = 100,nstart=100)
  cer_km <- CER(kmres$cluster,data$Misslabel)
  
  #miceres <- fun_mice(x=data$Missing, k=k)  
  #cer_mice <- CER(miceres$cluster,data$Misslabel)
  
  scimputeres <- kmeans(data$scImpute_data, centers=k,iter.max = 100,nstart=100)
  cer_scimpute <- CER(scimputeres$cluster,data$Misslabel)
  
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
               cer_km=cer_km, cer_scimpute=cer_scimpute, cer_kpod=cer_kpod, cer_pkpod=cer_pkpod,
               instability_values=instability_values, instability_values_sd=instability_values_sd
  ))
}

### load data -----------
realdata <- read.table("C:/Rstudio/kpod_MNAR/real_incomplete_data/GSE59739_cleardata_significant.txt", header = TRUE, sep = "\t")
data <- list(Missing=NULL, Misslabel=realdata[,2])
realdata <- as.matrix(realdata[,-c(1:2)])
realdata[realdata==0] <- NA
#delete fetures whose missingness>600
manynarows <- which(apply(realdata,2,function(z) sum(is.na(z)) )>600)
realdata <- realdata[,-manynarows]
dim(realdata)
data$Missing <- apply(realdata, 2 , function(z) (z-mean(z,na.rm=TRUE))/sd(z,na.rm=TRUE)    )

### impute by scImpute -----------
scImpute::scimpute(
  count_path = "C:/Rstudio/kpod_MNAR/scimpute/GSE59739_scimpute_input.csv", infile = "csv",
  outfile = "csv", out_dir = "C:/Rstudio/kpod_MNAR/GSE59739_scimpute_output_",
  labeled = FALSE,  drop_thre = 0.5,  Kcluster = 4, ncores = 1)
scImpute_data<- read.table("C:/Rstudio/kpod_MNAR/scimpute/GSE59739_scimpute_output_scimpute_count.csv", header = TRUE, sep = ",")
scImpute_data <- t(scImpute_data[,-c(1)])
scImpute_data <- apply(scImpute_data, 2, function(z) (z-mean(z,na.rm=TRUE))/sd(z,na.rm=TRUE) )
data$scImpute_data <- scImpute_data[,-manynarows]


### experiments -------------

set.seed(59739)
tuninglambdalist=c(0.001,0.01,0.1,1,10)
instability_values <- rep(0,length(tuninglambdalist))
instability_values_sd <- rep(0,length(tuninglambdalist))
for (t in 1:length(tuninglambdalist)){
  instability_res <- instability_1218(x=data$Missing,k=4,tuningpara = tuninglambdalist[t], 
                                      clusteringfun = pkpod_1218, nstart=10, perm.max = 10  ) 
  instability_values[t] <- instability_res$Svalue
  instability_values_sd[t] <- instability_res$Svalue_sd
}
selectedlambda <- tuninglambdalist[which.min(instability_values)]
if ( selectedlambda==tail(tuninglambdalist,1) ){
  lowerbound <- tail(instability_values,1) + tail(instability_values_sd,1)
  selectedlambda <- tuninglambdalist[min(which( instability_values < lowerbound ))]
}

ExperimentResults_realdata_GSE59739$tuning <- list(
  tuninglambdalist=c(0.001,0.01,0.1,1,10), 
  instability_values=instability_values, instability_values_sd=instability_values_sd,
  lowerbound=lowerbound, selectedlambda=selectedlambda
)

ExperimentResults_realdata_GSE59739$allresults <- matrix(0,30,4)
data$zeroimpt <- data$Missing
data$zeroimpt[is.na(data$Missing)] <- 0
for (rept in 1:30){
  set.seed(2026+rept*9)
  kmres <- kmeans(data$zeroimpt,centers=4,iter.max = 100,nstart=100)
  ExperimentResults_realdata_GSE59739$allresults[rept,1] <- CER(kmres$cluster,data$Misslabel)
  
  scimputeres <- kmeans(data$scImpute_data, centers=4,iter.max = 100,nstart=100)
  ExperimentResults_realdata_GSE59739$allresults[rept,2] <- CER(scimputeres$cluster,data$Misslabel)
  
  kpodres <- pkpod_1218(x=data$Missing,k=4,lambda=0.000001,nstart=100)  
  ExperimentResults_realdata_GSE59739$allresults[rept,3] <- CER(kpodres$cluster,data$Misslabel)
  
  pkpodres <- pkpod_1218(x=data$Missing,k=4,lambda=ExperimentResults_realdata_GSE59739$tuning$selectedlambda,nstart=100) 
  ExperimentResults_realdata_GSE59739$allresults[rept,4] <- CER(pkpodres$cluster,data$Misslabel)
  
  if(rept%%10==0){print(rept)}
}


### summarize results -----------
print(apply(ExperimentResults_realdata_GSE59739$allresults,2,mean))
print(apply(ExperimentResults_realdata_GSE59739$allresults,2,sd))


### draw figures ------------ 
# using the last repetition of "ExperimentResults_realdata_GSE59739$allresults"

ExperimentResults_realdata_GSE59739$labels_for_figures = list(
  original_labels=cbind( data$Misslabel, kmres$cluster, scimputeres$cluster, kpodres$cluster, pkpodres$cluster ),
  aligned_labels=matrix(0,length(data$Misslabel),5)
)
ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,1] <- data$Misslabel


## align labels for color

#kmeans
tab <- table(data$Misslabel, kmres$cluster)
map <- apply(tab, 2, function(z) names(which.max(z)))
ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,2] <- map[as.character(kmres$cluster)]
#scimpute
tab <- table(data$Misslabel, scimputeres$cluster)
map <- apply(tab, 2, function(z) names(which.max(z)))
ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,3] <- map[as.character(scimputeres$cluster)]
#kpod
tab <- table(data$Misslabel, kpodres$cluster)
map <- apply(tab, 2, function(z) names(which.max(z)))
ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,4] <- map[as.character(kpodres$cluster)]
#pkpod
tab <- table(data$Misslabel, pkpodres$cluster)
map <- apply(tab, 2, function(z) names(which.max(z)))
ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,5] <- map[as.character(pkpodres$cluster)]


## make mis-clustered point to be black

ExperimentResults_realdata_GSE59739$labels_for_figures$mis_clustered <- ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels
for (j in 2:5){
  for (i in 1:622){
    if ( (ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[i,j]) != (ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[i,1])  ){
      ExperimentResults_realdata_GSE59739$labels_for_figures$mis_clustered[i,j] <- "black"
    }
  }
}


## run UMAP 

umap_res <- umap(data$zeroimpt)
col_map <- c(
  "NF"  = "orange",
  "NP"  = "red",
  "PEP" = "green4",
  "TH"  = "blue",
  "black" = "black"
)
pch_map <- c(
  "NF"  = 0,  #square
  "NP"  = 1,  #circle
  "PEP" = 2,  #triangle
  "TH"  = 5   #diamond
)
par(mar=c(4,4,1,1))
plot(umap_res$layout,
     pch=pch_map[ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,1]], 
     col=col_map[ExperimentResults_realdata_GSE59739$labels_for_figures$mis_clustered[,2]],
     lwd=1.2, cex=0.8, xlab="UMAP 1",ylab="UMAP 2")
legend("topleft", legend = names(pch_map), pch= pch_map, col= col_map, cex=0.8)
plot(umap_res$layout,
     pch=pch_map[ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,1]], 
     col=col_map[ExperimentResults_realdata_GSE59739$labels_for_figures$mis_clustered[,3]],
     lwd=1.2, cex=0.8, xlab="UMAP 1",ylab="UMAP 2")
plot(umap_res$layout,
     pch=pch_map[ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,1]], 
     col=col_map[ExperimentResults_realdata_GSE59739$labels_for_figures$mis_clustered[,4]],
     lwd=1.2, cex=0.8, xlab="UMAP 1",ylab="UMAP 2")
plot(umap_res$layout,
     pch=pch_map[ExperimentResults_realdata_GSE59739$labels_for_figures$aligned_labels[,1]], 
     col=col_map[ExperimentResults_realdata_GSE59739$labels_for_figures$mis_clustered[,5]],
     lwd=1.2, cex=0.8, xlab="UMAP 1",ylab="UMAP 2")



