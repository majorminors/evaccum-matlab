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

# Highlight and run the above functions to load them into the R environment

# Source model functions
source("dmc/dmc.R")
load_model("LBA", "lba_B.R")

# Load your samples object here 
load("samples/samples_LBA_model1.RData")

# First we need to define a function to extract only the parameters/contrast we want:
# This takes the difference in thresholds between cond1 and cond2 
# You'll need to change these to your own parameter names and the contrasts that
# you're interested in. 
fun <- function (thetas) {
  (cat(thetas[,"B.cohE.angleE",, drop=F],thetas[,"B.cohE.angleH",, drop=F],thetas[,"B.cohH.angleE",, drop=F],thetas[,"B.cohH.angleH",, drop=F])) # this is just a shorthand for print
}

# Then we can extract all sampling iterations using group.inference.dist()
# this produces a  mean that differs from the mean we get in paramSummary, so I'm not sure what it's doing differently
#all_samples <- as.numeric(group.inference.dist(hsampling_LBA, fun))
#head(all_samples)
#mean(all_samples)

# Or extract individual-subject average effect using subject.inference.dist()
# however, this gets a mean that is the same subject-wise!
subject_means <- subject.inference.dist(hsampling_LBA, fun)
head(subject_means)
mean(subject_means)




