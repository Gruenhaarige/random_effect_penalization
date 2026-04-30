require(MASS)
require(lme4)
require(glmmTMB)
require(methods)

# getMethod("getExponentialTerm", "bmerLassoDist")
glmmTMB:::.valid_prior
# version

source("generateData.R")
source("glmmLassoSim.R")
source("../fitGlmmLasso.R")
source("../fitModel.R")

job_id <- as.integer(Sys.getenv("PBS_ARRAYID"))
settings <- readRDS("glmmTMB_jobs295.rds")
settings <- subset(settings, jobid == job_id)
n_set <- nrow(settings)

for(j in seq_len(n_set)){
  set.seed(settings$seed[j])
  
  now <- Sys.time()
  
  results <- glmmLassoSim(n = settings$n[j],
                          p1 = settings$p1[j],
                          p2 = settings$p2[j],
                          intercept = 0,
                          beta_part = c(-0.5, 0.5),
                          covX_strength = settings$covX_strength[j],
                          sigma_eps = settings$sigma_eps[j],
                          slope = settings$slope[j],
                          n_grouping_var = settings$grV[j],
                          n_groups = settings$n_groups[j],
                          covU_strength = settings$covU_strength[j],
                          method = "glmmTMB",
                          max_iter = 1e5)
  
  out <- list(job_no = job_id,
              seed = settings$seed[j],
              setting = settings[j,],
              cluster_info = Sys.info(),
              runtime = difftime(Sys.time(), now),
              results = results)
  
  saveRDS(out, file = paste0("run3/results/glmmTMB_",
                             sprintf("%02d", settings$id[j]), "_",
                             sprintf("%02d", settings$seed[j]), ".rds"))
  
  if(j %% 10 == 0){
    message(j, " / ", n_set)
  }
}