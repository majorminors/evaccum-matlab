set.seed(123)#for reproducibility
# save the root directory
wd <- here::here();

# move to the DMC directory to lead the model 
setwd("/imaging/woolgar/projects/Dorian/evaccum/evaccum-lba-dmc/DMC_190819/")
source("dmc/dmc.R")
load_model (dir_name="LBA",model_name="lba_B.R")

# go back to the root
setwd(wd)

##################################
# create our model design matrix # 
##################################
# this will make 34 models
# 1 = fixed
# 2 = varying by coherence
# 3 = varying by angle
# 4 = varying by both coherence and angle
levs <- list("1","Coh","Angle",c("Coh","Angle"))
levs2 <- list(c(0,0),c(1,0),c(0,1),c(1,1))
design_matrix <- expand.grid(c(1:3),c(1:3),c(1:3))
design_matrix <- design_matrix[-1,] %>% rbind(expand.grid(c(1,4),c(1,4),c(1,4))) %>%  as.data.frame()
colnames(design_matrix) <- c('B','mean_v','t0')
design_matrix['Count'] <- 1:nrow(design_matrix)

theFiles = Sys.glob(file.path(here::here('jobs/hierarchical/'), "hsampling_LBA*_LBA.RData"))
library(gtools)
theFilesSorted <- mixedsort(theFiles)

models <- vector(mode = "list", length = length(theFiles))
modelResults <- as.data.frame(matrix(nrow=length(theFiles),ncol=6))
colnames(modelResults) <- c("rhat","dic","min_dev","B","mean_v","t0")

for (fileNumber in 1:length(theFiles)) {

  load(theFiles[fileNumber])
  
  models[fileNumber] <- hsampling_LBA
  
  modelResults$B[fileNumber] <- design_matrix$B[fileNumber]
  modelResults$mean_v[fileNumber] <- design_matrix$mean_v[fileNumber]
  modelResults$t0[fileNumber] <- design_matrix$t0[fileNumber]
  
  # check convergence for each chain (rhat)
  tmp <- gelman.diag.dmc(hsampling_LBA)
  tmp <- unlist(lapply(tmp,function(x){x$mpsrf}))
  modelResults$rhat[fileNumber] <- mean(tmp)
  # larger than 1.1 is something to look at
  
  tmp <- h.IC.dmc(hsampling_LBA,DIC=TRUE)
  tmp <- colSums(tmp)
  modelResults$min_dev[fileNumber] <- tmp[1]
  modelResults$dic[fileNumber] <- tmp[2]

  
  
  # all chains converged well
  #   1    6    8   10    2    3    7    5    9    4 
  #1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 1.01 
  #Mean
  #[1] 1.01
  
  # visual inspection of chains
  #for (plot in 1:1) {
  #  plot.dmc(sampling_LBA,density=FALSE,smooth=FALSE,subject=i,layout=c(2,5))
    #  Sys.sleep(0.1) 
  #}
  
  # Use get.thin to look at mean and minimum effective sample sizes and an estimate of what thinning might be needed to remove all autocorrelation as actual n divided by either minimum or mean effective n.
  #get.thin(sampling_LBA)
  
  # Thinning and sample size are appropriate
  # Minimum Effective Size
  #    1    2    3    4    5    6    7    8    9   10 
  # 4076 4292 4120 4389 4274 4126 4339 4128 4404 4308 
  # Mean Effective Size
  #    1    2    3    4    5    6    7    8    9   10 
  # 4377 4399 4258 4517 4448 4389 4526 4420 4561 4600 
  # Thin
  #      1 2 3 4 5 6 7 8 9 10
  # mean 2 2 3 2 2 2 2 2 2  2
  # min  3 3 3 2 3 3 3 3 2  3
  
  


  
}

modelResults <- modelResults[order(modelResults['dic'],decreasing=FALSE),]

#    rhat      dic     min_dev  B mean_v t0
# 23 1.032174 14008.38 13253.81 3      2  3

# now let's look at the winning parameters

winningModel <- 23
load(theFiles[winningModel])
# double check it
gelman.diag.dmc(hsampling_LBA)
h.IC.dmc(hsampling_LBA,DIC=TRUE)

# this wont work because you don't have the pop.mean in this file
#h.check.recovery.dmc(hsampling_LBA,pop.mean,ptype=1,hyper=TRUE)

# summary of the hypers - h1 is mean params, h2 is the scale
paramSummary <- summary.dmc(hsampling_LBA,hyper=TRUE)
meanParams <-paramSummary$statistics
round(paramSummary$quantiles[,c(1,3,5)],3)
effectiveSize.dmc(hsampling_LBA,hyper=TRUE)
# summary of the participants
paramSummary2 <- summary.dmc(hsampling_LBA)


# ok, let's extract the parameters now
# Gets all sampling iterations of all chains for all subjects
group.inference.dist <- function (hsamples, fun) {
  # Bring longer thetas down to minimum nmc by sampling
  nmcs <- sapply(hsamples, function(x) x$nmc)
  nmc  <- min(nmcs)
  for (i in 1:length(hsamples)) if (nmcs[i] > nmc) {
    hsamples[[i]]$theta <- 
      hsamples[[i]]$theta[, , sample(1:dim(hsamples[[i]]$theta)[3], nmc)]
  }
  # inference <- lapply(hsamples, function(x) x["theta"])
  inference <- list()
  for (i in 1:length(hsamples)) {
    thetas <- hsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  apply(inf2, c(3, 1), mean)
}
# Gets average effect for each individual (i.e., one value per subject)
subject.inference.dist <- function (hsamples, fun) {
  # Bring longer thetas down to minimum nmc by sampling
  nmcs <- sapply(hsamples, function(x) x$nmc)
  nmc  <- min(nmcs)
  for (i in 1:length(hsamples)) if (nmcs[i] > nmc) {
    hsamples[[i]]$theta <- 
      hsamples[[i]]$theta[, , sample(1:dim(hsamples[[i]]$theta)[3], nmc)]
  }
  # inference <- lapply(hsamples, function(x) x["theta"])
  inference <- list()
  for (i in 1:length(hsamples)) {
    thetas <- hsamples[[i]]$theta
    inference [[i]] <- fun (thetas)
  }
  inf2 <- unlist(inference)
  dim3 <- c(dim(inference[[1]]), length(inf2)/prod(dim(inference[[1]])))
  dim(inf2) <- dim3
  effect <- apply(inf2, c(4), mean)
  names(effect) <- names(hsamples)
  data.frame(effect)
}


# boundary
fun <- function (thetas) {
  (thetas[,"B.cohE.angleE",, drop=F]) # this is just a shorthand for print
}
subject_means_b_ecer <- subject.inference.dist(hsampling_LBA, fun)
all_samples_b_ecer <- as.numeric(group.inference.dist(hsampling_LBA, fun))

fun <- function (thetas) {
  (thetas[,"B.cohE.angleH",, drop=F]) 
}
subject_means_b_echr <- subject.inference.dist(hsampling_LBA, fun)
all_samples_b_echr <- as.numeric(group.inference.dist(hsampling_LBA, fun))

fun <- function (thetas) {
  (thetas[,"B.cohH.angleE",, drop=F]) 
}
subject_means_b_hcer <- subject.inference.dist(hsampling_LBA, fun)
all_samples_b_hcer <- as.numeric(group.inference.dist(hsampling_LBA, fun))

fun <- function (thetas) {
  (thetas[,"B.cohH.angleH",, drop=F]) 
}
subject_means_b_hchr <- subject.inference.dist(hsampling_LBA, fun)
all_samples_b_hchr <- as.numeric(group.inference.dist(hsampling_LBA, fun))

# winning accumulator
fun <- function (thetas) {
  (thetas[,"mean_v.cohE.angleE.true",, drop=F]) 
}
subject_means_v1_ecer <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v1_ecer <- as.numeric(group.inference.dist(hsampling_LBA, fun))
fun <- function (thetas) {
  (thetas[,"mean_v.cohE.angleH.true",, drop=F]) 
}
subject_means_v1_echr <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v1_echr <- as.numeric(group.inference.dist(hsampling_LBA, fun))
fun <- function (thetas) {
  (thetas[,"mean_v.cohH.angleE.true",, drop=F]) 
}
subject_means_v1_hcer <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v1_hcer <- as.numeric(group.inference.dist(hsampling_LBA, fun))
fun <- function (thetas) {
  (thetas[,"mean_v.cohH.angleH.true",, drop=F]) 
}
subject_means_v1_hchr <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v1_hchr <- as.numeric(group.inference.dist(hsampling_LBA, fun))

# sd for winning accumulator (?)
fun <- function (thetas) {
  (thetas[,"sd_v.true",, drop=F]) 
}
subject_means_sdv1 <- subject.inference.dist(hsampling_LBA, fun)
all_samples_sdv1 <- as.numeric(group.inference.dist(hsampling_LBA, fun))

# losing accumulator
fun <- function (thetas) {
  (thetas[,"mean_v.cohE.angleE.false",, drop=F]) 
}
subject_means_v0_ecer <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v0_ecer <- as.numeric(group.inference.dist(hsampling_LBA, fun))
fun <- function (thetas) {
  (thetas[,"mean_v.cohE.angleH.false",, drop=F]) 
}
subject_means_v0_echr <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v0_echr <- as.numeric(group.inference.dist(hsampling_LBA, fun))
fun <- function (thetas) {
  (thetas[,"mean_v.cohH.angleE.false",, drop=F]) 
}
subject_means_v0_hcer <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v0_hcer <- as.numeric(group.inference.dist(hsampling_LBA, fun))
fun <- function (thetas) {
  (thetas[,"mean_v.cohH.angleH.false",, drop=F]) 
}
subject_means_v0_hchr <- subject.inference.dist(hsampling_LBA, fun)
all_samples_v0_hchr <- as.numeric(group.inference.dist(hsampling_LBA, fun))

# non-decision time
fun <- function (thetas) {
  (thetas[,"t0",, drop=F]) 
}
subject_means_t0 <- subject.inference.dist(hsampling_LBA, fun)
all_samples_t0 <- as.numeric(group.inference.dist(hsampling_LBA, fun))

samplingEstimates <- data.frame(mean(all_samples_b_ecer),
                                mean(all_samples_b_echr),
                                mean(all_samples_b_hcer),
                                mean(all_samples_b_hchr),
                                mean(all_samples_v1_ecer),
                                mean(all_samples_v1_echr),
                                mean(all_samples_v1_hcer),
                                mean(all_samples_v1_hchr),
                                mean(all_samples_sdv1),
                                mean(all_samples_v0_ecer),
                                mean(all_samples_v0_echr),
                                mean(all_samples_v0_hcer),
                                mean(all_samples_v0_hchr),
                                mean(all_samples_t0))
names(samplingEstimates)[1] <- "b_ecer"
names(samplingEstimates)[2] <- "b_echr"
names(samplingEstimates)[3] <- "b_hcer"
names(samplingEstimates)[4] <- "b_hchr"
names(samplingEstimates)[5] <- "v1_ecer"
names(samplingEstimates)[6] <- "v1_echr"
names(samplingEstimates)[7] <- "v1_hcer"
names(samplingEstimates)[8] <- "v1_hchr"
names(samplingEstimates)[9] <- "sdv1"
names(samplingEstimates)[10] <- "v0_ecer"
names(samplingEstimates)[11] <- "v0_echr"
names(samplingEstimates)[12] <- "v0_hcer"
names(samplingEstimates)[13] <- "v0_hchr"
names(samplingEstimates)[14] <- "t0"
write.csv(samplingEstimates,paste(wd,'meanParams.csv', sep = "/"))
          
combinedData <- data.frame(subject_means_b_ecer,
                           subject_means_b_echr,
                           subject_means_b_hcer,
                           subject_means_b_hchr,
                           subject_means_v1_ecer,
                           subject_means_v1_echr,
                           subject_means_v1_hcer,
                           subject_means_v1_hchr,
                           subject_means_sdv1,
                           subject_means_v0_ecer,
                           subject_means_v0_echr,
                           subject_means_v0_hcer,
                           subject_means_v0_hchr,
                           subject_means_t0)
names(combinedData)[names(combinedData) == "effect"] <- "b_ecer"
names(combinedData)[names(combinedData) == "effect.1"] <- "b_echr"
names(combinedData)[names(combinedData) == "effect.2"] <- "b_hcer"
names(combinedData)[names(combinedData) == "effect.3"] <- "b_hchr"
names(combinedData)[names(combinedData) == "effect.4"] <- "v1_ecer"
names(combinedData)[names(combinedData) == "effect.5"] <- "v1_echr"
names(combinedData)[names(combinedData) == "effect.6"] <- "v1_hcer"
names(combinedData)[names(combinedData) == "effect.7"] <- "v1_hchr"
names(combinedData)[names(combinedData) == "effect.8"] <- "sdv1"
names(combinedData)[names(combinedData) == "effect.9"] <- "v0_ecer"
names(combinedData)[names(combinedData) == "effect.10"] <- "v0_echr"
names(combinedData)[names(combinedData) == "effect.11"] <- "v0_hcer"
names(combinedData)[names(combinedData) == "effect.12"] <- "v0_hchr"
names(combinedData)[names(combinedData) == "effect.13"] <- "t0"



write.csv(combinedData,paste(wd,'subjectParams.csv', sep = "/"))

