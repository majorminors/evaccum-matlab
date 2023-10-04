set.seed(123)#for reproducibility

# dependencies
install.packages("here")
install.packages("pacman")
require(devtools)
install_version("gsl", version = "1.9-10", repos = "http://cran.us.r-project.org")
install.packages("rtdists")
# devtools::install_github("adzemski/rtnorm")
require(rtnorm)
install.packages("truncdist")
install.packages("numDeriv")
install.packages("vioplot")
install.packages("snowfall")
install.packages("rlecuyer")
install.packages("mcmc")
install.packages("statmod")
install.packages("pracma")
install.packages("coda")
install.packages("ggplot2")
install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz", repo=NULL, type="source")

# save the root directory
wd <- here::here();

# move to the DMC directory to lead the model 
setwd("/imaging/woolgar/projects/Dorian/evaccum/evaccum-lba-dmc/DMC_190819/")
source("dmc/dmc.R")
load_model (dir_name="LBA",model_name="lba_B.R")

# go back to the root
setwd(wd)

# load some required packages
pacman::p_load('tidyverse','magrittr','slurmR')

# remove old jobs
unlink(here::here('jobs', 'hierarchical_coh'),recursive = TRUE)
# set up folder for slurm jobs
library(slurmR)
opts_slurmR$set_tmp_path(here::here("jobs/"))
opts_slurmR$set_job_name("hierarchical_coh")

###############################################################
# first we take in the data and do some additional processing #
###############################################################

data <- read.csv(here::here("/imaging/woolgar/projects/Dorian/evaccum/evaccum-matlab/data/meg_pilot_4/behavioural/DMC_table_coh.csv"),na.strings=c("","NA")) %>% as.data.frame()

# get indices of  correct  trials and remove implausible RTs
data <- data %>% mutate(Correct = ifelse((S=="s1" & R == "r1") | (S=="s2" & R == "r2"),"Correct","Wrong"))  %>%
  filter(!is.na(R),.preserve = TRUE) %>%
  mutate(Plausibility = case_when(
    Correct == "Correct" & RT < 0.25 ~ "Too quick",
    Correct == "Correct" & RT > 4.5 ~ "Too slow",
    TRUE ~ "Good")) %>%
  filter(Plausibility == "Good",.preserve = TRUE)

# remove outlier RTs (mean ± 2.5sd) on  a subject-by-subject basis
data <- data %>% group_by(s) %>% mutate(Plausibility = case_when(
  
  Correct == "Correct" & RT < (mean(RT)-2*sd(RT)) ~ "Outlier",
  Correct == "Correct" & RT > (mean(RT)+2*sd(RT)) ~ "Outlier",
  TRUE ~ "Good")) %>% 
  # filter out outliers and remove  the last two columns
  ungroup() %>% filter(Plausibility == "Good",.preserve = TRUE) %>%
  select(-c("Correct","Plausibility")) 


# tell R to encode subject and factor columns as factors (this is an R thing)
data <- data %>% mutate( s = as.factor(s),
                         S = as.factor(S),
                         Coh   = as.factor(Coh),
                         R = as.factor(R))

# make a function to submit as a slurm job
# it reads in our preprocessed data (`data`)
# it also takes an index we'll use to cycle through model variants and generate savenames
LBAFit <- function(i,data) {
  
  # so we'll set the project root again, since this will be done by a worker with a confusing level of access to the R environment
  wd <- here::here();
  
  # and once again load the model we want
  setwd("/imaging/woolgar/projects/Dorian/evaccum/evaccum-lba-dmc/DMC_190819/")
  source("dmc/dmc.R")
  load_model (dir_name="LBA",model_name="lba_B.R")
  
  # before going back to the project root
  setwd(wd)
  
  # once again load some required stuff
  pacman::p_load('tidyverse','magrittr')
  
  
  ##################################
  # create our model design matrix # 
  ##################################
  # this will make 34 models
  # 1 = fixed
  # 2 = varying by coherence
  
  levs <- list("1","Coh")
  levs2 <- list(c(0,0),c(1,0))
  design_matrix <- expand.grid(c(1:2),c(1:2),c(1:2))
  design_matrix <- design_matrix[-1,] %>% rbind(expand.grid(c(1:2),c(1:2),c(1:2))) %>%  as.data.frame()
  colnames(design_matrix) <- c('B','mean_v','t0')
  design_matrix['Count'] <- 1:nrow(design_matrix)
  
  ##################################################
  # now we set up the model as in the dmc tutorial #
  ##################################################
  
  # Define factorial design 2 directions x 2 coherence levels
  # s1/s2 left and right motion direction
  # cohE/cohH easy/hard coherence
  factors <- list(S=c("s1","s2"),Coh=c("cohE","cohH"))
  # specify the responses
  responses <- c("r1","r2") # keys o and p
  
  # define a new factor M (Match) with correct and incorrect responses for each factor level
  match.map <- list(M=list(s1="r1",s2="r2")) 
  
  # define constants to allow parameter identification
  const <- c(sd_v.false=1,st0=0)
  
  
  # Specify how the factors map onto the parameters re: the design matrix (remember we use the index to cycle through these, so each worker will do one aspect of this)
  # e.g. boundary  changes with manipulations, accum-rate faster for correct
  pB <-levs[design_matrix$B[i]] %>% unlist()
  pV <-levs[design_matrix$mean_v[i]] %>% unlist()
  pV <- if(pV=='1') 'M' else c(pV,'M')
  pT0 <- levs[design_matrix$t0[i]] %>% unlist()
  
  p.map <- list(A="1",B=pB,t0=pT0,mean_v=pV,sd_v="M",st0="1")
  
  model1 <- model.dmc(p.map,match.map = match.map,constants=const,
                      responses=responses,factors=factors)
  
  ###########################################################################
  # now define our prior distributions---I think this is somewhat arbitrary #
  ###########################################################################
  
  vec_p <- function(param,factor1 = 1, value,plower,pupper)
  {
    if(!factor1){
      lbls <- param  
    }
    else if(factor1) {
      lbls <- do.call(paste, c(expand.grid(param, c("cohE","cohH")), list(sep='.')))
      
    }
    
    parv <- rep(value,length(lbls))
    names(parv) <- lbls
    plow <- rep(plower,length(parv))
    phigh <- rep(pupper,length(parv))
    
    vec <- list("parv"=parv,"lower"=plow,"upper"=phigh)
    return(vec)
    
  }
  
  vec_pV <- function(param,factor1 = 1, value_true,value_false,plower,pupper)
  {  
    true_false <- c("true","false")
    if(!factor1){
      lbls <- do.call(paste, c(expand.grid(param,true_false),list(sep ='.')))  
    }
    else if(factor1) {
      lbls <- do.call(paste, c(expand.grid(param, c("cohE","cohH"),true_false), list(sep='.')))
    }
    
    parv <- c(rep(value_true,length(lbls)/2),rep(value_false,length(lbls)/2))
    names(parv) <- lbls
    
    
    plow <- rep(plower,length(parv))
    phigh <- rep(pupper,length(parv))
    
    vec <- list("parv"=parv,"lower"=plow,"upper"=phigh)
    return(vec)
    
  }
  
  comp_vec <- function(vA,vB,vV,vsd,vt0) {
    pars<-c(vA$parv,vB$parv,vV$parv,vsd$parv,vt0$parv)
    upper<-c(vA$upper,vB$upper,vV$upper,vsd$upper,vt0$upper)
    lower<-c(vA$lower,vB$lower,vV$lower,vsd$lower,vt0$lower)
    comp <- list("pars" = pars,"upper"=upper,"lower"=lower)
    return(comp)
  }
  
  B <-levs2[design_matrix$B[i]] %>% unlist()
  V <-levs2[design_matrix$mean_v[i]] %>% unlist()
  T0 <- levs2[design_matrix$t0[i]] %>% unlist()
  
  pop_prior <- comp_vec(vec_p('A',0,0, .4,0,NA),
                        vec_p('B',B[1],0, 4,0,NA),
                        vec_pV('mean_v',V[1],0, 1,-1,NA,NA),
                        vec_p('sd_v.true',0,0, .3,0,NA),
                        vec_p('t0',T0[1],0, .4,.1,1))
  
  pop_scale <- comp_vec(vec_p('A',0,0, .1,0,NA),
                        vec_p('B',B[1],0, .2,0,NA),
                        vec_pV('mean_v',V[1],0, .2,.2,NA,NA),
                        vec_p('sd_v.true',0,0, .1,0,NA),
                        vec_p('t0',T0[1],0, .05,.1,1))
  
  
  
  
  # here, we double check that parameters match with model
  pop.mean <- pop_prior$pars
  check.p.vector(pop.mean,model1)
  pop.scale <- pop_scale$pars
  check.p.vector(pop.scale,model1)
  
  nparams = length(pop.mean)
  paramnames = names(pop.mean)
  
  # and we save it
  
  p.prior <- prior.p.dmc(dists = rep("tnorm",nparams),
                         p1=pop.mean,p2=pop.scale,
                         lower = pop_prior$lower,
                         upper = pop_prior$upper)
  
  
  # bind the data to the model object
  data.model <- data.model.dmc(data, model1)
  
  ######################################
  # now we do the actual model fitting #
  ######################################
  
  # here we are fitting the model non-hierarchically
  # so subjects are fixed effects, which I guess means we assume that are sampling is representative
  # but I'm not so sure about that, given the variability in performance in coherence thresholding
  # and the fact we don't threshold for very long
  # the point is mute however because I can't figure out how to generate a DIC for these models, just for each subject
  # although I think it's just a case of getting the DIC for each subject and summing them---we can ask Alessandro about this
  # but we do this first before fitting hierarchically (where I think I can generate a DIC using a provided function)
  
  # this specificies cores for the snowfall parallel processing package in R---I think we don't actually need
  # this now, since I'm running it on the cluster, but I won't remove it because I'm scared
  mycores <- 10
  
  # set up our model fitting object
  sampling_LBA <- h.samples.dmc(
    nmc = 400, # number of samples left AFTER thinning
    p.prior = p.prior, # list of priors
    thin = 10, #thinning to deal with autocorrelation
    data = data.model # note: Default number of chains is 3 * nparams
  )
  
  # run the model using the automatic function h.RUN.dmc
  sampling_LBA <- h.RUN.dmc(sampling_LBA, cores = mycores, verbose = TRUE)
  
  # save this sampling object---this is one model, all subjects
  savename =   sprintf('sampling_LBA%i_LBA.RData',i)
  save(sampling_LBA, file = here::here('jobs/hierarchical_coh/',savename))
  
  # so here we'll go on to fit the model hierarchically
  # this uses the results above to generate a starting point
  # here subjects are treated as a random effect
  
  groupstart <- make.hstart(sampling_LBA) # group-level parameters
  subjectstart <- make.theta1(sampling_LBA) # subject-level parameters
  # 
  # # If we perform hierarchical modelling, we assume that each parameter has a group-level
  # # distribution with some location (e.g. mean) and scale (e.g. SD). These group-level distributions
  # # then serve as priors for our parameters at the level of individual subjects.
  # # Thus, we need to define prior distributions for both the group-level mean and the group-level SD
  # # of each parameter. These are also known as hyperpriors.
  # 
  # # For the prior distributions for the location (i.e. mean) of our group-level distributions, we
  # # will simply recycle the priors we used for the fixed-effects modelling
  mu.prior <- p.prior
  # 
  # # What's new with hierarchical modelling is we also have to define prior distributions for the
  # # scale (i.e. SD) of our group-level distributions
  sigma.prior <- prior.p.dmc(
    # set the type of distribution: Gamma distribution for all parameters
    dists = rep(x = "gamma", times = nparams),
    # set the shape parameter (p1) of the distribution: Shape = 1, hence exponential distribution
    p1 = c(
      setNames(
        object = rep(x = 1, times = nparams),
        nm = paramnames
      )
    ),
    # set the scale parameter (p2) of the distribution
    p2 = c(
      setNames(
        object = rep(x = 1, times = nparams),
        nm = paramnames
      )
    )
  )
  
  # Lastly, we combine these into a hyper-prior list
  pp.prior <- list(mu.prior, sigma.prior)
  
  # set up our model fitting object
  hsampling_LBA <- h.samples.dmc(
    nmc = 400, # number of samples AFTER thinning
    pp.prior = pp.prior, # hyper-priors
    p.prior = mu.prior, # prior for each subject
    hstart.prior = groupstart,
    theta1 = subjectstart,
    thin = 10, # thinning to deal with autocorrelation
    data = data.model # note: Default number of chains is 3 * nparams
  )
  
  # Now we run the models using some automated functions from the DMC toolbox
  hsampling_LBA <- h.run.unstuck.dmc(
    samples = hsampling_LBA, cores = mycores, report = 10,
    p.migrate = 0.05, h.p.migrate = 0.05)
  hsampling_LAB <- h.run.converge.dmc(
    samples = h.samples.dmc(samples = hsampling_LBA, nmc = 120, thin = 25),
    nmc = 40, cores = mycores, report = 10, verbose = TRUE)
  
  # after the model has converged, get some "clean" samples to be used for parameter estimation
  hsampling_LBA <- h.run.dmc(
    h.samples.dmc(samples = hsampling_LBA, nmc = 500),
    report = 1, cores = mycores)
  hsampling_CamCAN <- h.run.dmc(
    h.samples.dmc(samples = hsampling_LBA, nmc = 500),
    report = 1, cores = mycores)
  
  savename =   sprintf('hsampling_LBA%i_LBA.RData',i)
  save(hsampling_LBA, file = here::here('jobs/hierarchical_coh/',savename))
  
  # so we will now have results that we should be able to apply the model selection tutorial 5.3 to:
  # the function h.IC.dmc automatically obtains the DIC per subject and returns the sum
  
  
}

# submit jobs to cluster
# this is essentially a fancy wrapper for sbatch, but it does go via parallel:mclapply first
# we give this our LBAFit function, with the index and the data
# we also need to increase the time, so we can pass a command directly to sbatch with sbatch_opt = list(key=value, time=7200) time is in mins
job <- Slurm_lapply(c(1:34), LBAFit, data = data, njobs=34, mc.cores=10, plan = "collect", sbatch_opt = list(time=7200))
#check jobs exit status 
#status(job)
