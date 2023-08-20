# Just a quick script to extract the necessary stan code from the model we 
# specified earlier in stan_glmer

library(rstanarm)

mcmc_fit <- readRDS("posts/Variational MRP Pt7/data/fit_60k_meanfield.rds")

mcmc_fit$stanfit@stanmodel
