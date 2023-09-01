# honestly not a clue what this seed is for---could've been for slurm, could've been for parallel processing when I was trying to do it on my PC. Could be for something else entirely.
set.seed(123)#for reproducibility

#############################################################
# here I constructed anything that applied to all my models #
# and set up dependencies and stuff too                     #
#############################################################

# here I'm hard coding the unlinking of old, failed jobs
unlink(here::here('jobs', 'hierarchical'),recursive = TRUE)
# and then relinking the job folders for the next run
library(slurmR)
opts_slurmR$set_tmp_path(here::here("jobs/"))
opts_slurmR$set_job_name("hierarchical")
# not efficient but it worked

# now I'm constructing a function to submit my script as a slurm job
# it handed everything (e.g. `data`) to each of 34 models and then calculated the model fit for each and saved that for me to compile later
LBAFit <- function(i,data) {
  
    # I don't know if you'll need this, but I had trouble keeping track of where the worker was trying to work, so I'll leave it for you
    # so we'll set the project root again, since this will be done by a worker with a confusing level of access to the R environment
#  wd <- here::here();
  
 ######################################################################
 # so, in here you want to have your script that each worker will run #
 ######################################################################

  
}

# now I submit jobs to the cluster
# this is essentially a fancy wrapper for sbatch, but it does go via parallel:mclapply first
# we give this our LBAFit function, with the index and the data
# we also need to increase the time, so we can pass a command directly to sbatch with sbatch_opt = list(key=value, time=7200) time is in mins
job <- Slurm_lapply(c(1:34), LBAFit, data = data, njobs=34, mc.cores=10, plan = "collect", sbatch_opt = list(time=7200))
# you can check jobs' exit status like this
#status(job)
