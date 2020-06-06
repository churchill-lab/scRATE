#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop("Six arguments must be supplied (<cntfile> <fitfile> <outfile> <depth> <cores> <seed>).n", call.=FALSE)
}

library(scRATE)

cntfile <- args[1]
fitfile <- args[2]
outfile <- args[3]
scaler <- as.integer(args[4])
nCores <- as.integer(args[5])
seed <- as.integer(args[6])

load(cntfile)  # This will load 'cntmat', 'csize', and 'covar_list'
gsymb <- rownames(cntmat)
exposure <- log(csize)
log_scaler <- log(scaler / median(csize))

models <- c('P', 'NB', 'ZIP', 'ZINB')
results <- list()
simcnts <- list()

message(sprintf("\nRunning a simulation for %s", gsymb))
tryCatch({
  model_fit <- readRDS(fitfile)
  res <- list()
  res[['orig']] <- list()
  ysim_g <- list()
  for (m in 1:4) {
    # Generate simulated counts for each model
    message(sprintf('Simulation by %s', models[m]))
    flush.console()
    if(m < 3) {
      yrep <- posterior_predict(model_fit[[m]], offset = exposure + log_scaler)
    } else {
      yrep <- posterior_predict(model_fit[[m]], newdata = data.frame(exposure = exposure + log_scaler))
    }
    ysim <- as.numeric(tail(yrep, 1))
    # Fit simulated counts with all four models
    gexpr <- data.frame(y = ysim, exposure = exposure + log_scaler)
    model_fit_ysim <- fit_count_models(gexpr, 'y ~ 1', nCores, seed, adapt_delta=0.95)
    res[[paste('sim', m, sep = '')]][['elpd_loo']] <- compare_count_models(model_fit_ysim)
    res[[paste('sim', m, sep = '')]][['mean_par']] <- get_model_params(model_fit_ysim)
    ysim_g[[paste('sim', m, sep = '')]] <- ysim
  }
  results[[gsymb]] <- res
  simcnts[[gsymb]] <- ysim_g 
}, error = function(err) {
  message(sprintf("Error while fitting %s with %s model", gsymb, models[[m]]))
})

things2save <- list(modelfit=results, simcount=simcnts)
saveRDS(things2save, file = outfile)
