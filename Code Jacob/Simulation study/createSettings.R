n_seeds <- 50L
start_seed <- 2212

# settings per job
spj <- c(lme4 = 600,
         blme = 12,
         glmmTMB_25 = 150,
         glmmTMB_95 = 20,
         glmmTMB_295 = 4) # 16GB

raw_settings <- expand.grid(n = 100L,
                            n_groups = 15L,
                            p1 = 5L,
                            p2 = c(25L, 95L, 295L),
                            covX_strength = c(0, 0.3, 0.7),
                            covU_strength = 0.3,
                            grV = c(1L, 2L),
                            slope = c(0L, 2L),
                            sigma_eps = 0.2)

n_set <- nrow(raw_settings)
raw_settings <- raw_settings[order(raw_settings$p2, seq_len(n_set)),]

raw_settings <- cbind.data.frame(id = seq_len(n_set), raw_settings)

all_jobs <- do.call(rbind.data.frame,
                    lapply(seq(start_seed, start_seed + n_seeds - 1L), 
                           function(seed){
                             cbind.data.frame(raw_settings,
                                              seed = rep(seed, n_set))
                           })
)

lme4_jobs <- blme_jobs <- glmmTMB_jobs25 <- subset(all_jobs, p2 == 25L)
glmmTMB_jobs95 <- subset(all_jobs, p2 == 95L)
glmmTMB_jobs295 <- subset(all_jobs, p2 == 295L)

n25 <- nrow(lme4_jobs) # gleich für lme4, blme und glmmTMB25
n95 <- nrow(glmmTMB_jobs95)
n295 <- nrow(glmmTMB_jobs295)

assign_groups <- function(n_elements, n_groups) {
  base_size <- n_elements %/% n_groups
  remainder <- n_elements %% n_groups
  
  group_sizes <- c(rep(base_size + 1, remainder),
                   rep(base_size, n_groups - remainder))
  
  group_numbers <- rep(seq_len(n_groups), times = group_sizes)
  
  return(group_numbers)
}

lme4_jobs$jobid <- assign_groups(n25, ceiling(n25 / spj["lme4"]))
blme_jobs$jobid <- assign_groups(n25, ceiling(n25 / spj["blme"]))
glmmTMB_jobs25$jobid <- assign_groups(n25, ceiling(n25 / spj["glmmTMB_25"]))
glmmTMB_jobs95$jobid <- assign_groups(n25, ceiling(n25 / spj["glmmTMB_95"]))
glmmTMB_jobs295$jobid <- assign_groups(n25, ceiling(n25 / spj["glmmTMB_295"]))

saveRDS(lme4_jobs, file = "lme4_jobs.rds")
saveRDS(blme_jobs, file = "blme_jobs.rds")
saveRDS(glmmTMB_jobs25, file = "glmmTMB_jobs25.rds")
saveRDS(glmmTMB_jobs95, file = "glmmTMB_jobs95.rds")
saveRDS(glmmTMB_jobs295, file = "glmmTMB_jobs295.rds")