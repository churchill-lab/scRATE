#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop("Four arguments must be supplied (<cntfile> <outfile> <cores> <seed>).n", call.=FALSE)
}

library(scRATE)

cntfile <- args[1]
outfile <- args[2]
nCores <- as.integer(args[3])
seed <- as.integer(args[4])

load(cntfile)  # This will load 'cntmat', 'gsurv', 'csize', and 'ctype'
gname <- rownames(cntmat)
num_genes <- length(gsurv)
exposure <- log(csize)

# models <- c('P', 'NB', 'ZIP', 'ZINB')
models <- c('NB', 'ZINB')
results <- list()

for (gg in c(1:num_genes)) {
  if(gsurv[gg]) {
    y <- round(unlist(cntmat[gg,]))
    cat(sprintf("\nFitting models for %s\n", gname[gg]))
    tryCatch({
      model_fit <- fit_count_models(y, exposure, ctype, nCores, seed, model2fit=c(2,4), adapt_delta=0.95, brms4zi=TRUE)
      elpd_loo <- compare_count_models(model_fit)
      mean_par <- get_model_params(model_fit, ctyped = !is.null(ctype))
      res_g <- list()
      res_g[['orig']] <- list()
      res_g[['orig']][['elpd_loo']] <- elpd_loo
      res_g[['orig']][['mean_par']] <- mean_par
      sim_g <- list()
      for (m in models) {
        # Generate simulated counts for each model
        cat(sprintf('Simulation by %s\n', m))
        flush.console()
        yrep <- posterior_predict(model_fit[[m]], draws=1)
        ysim <- as.numeric(tail(yrep, 1))
        # Fit simulated counts with all four models
        sim_g[[m]] <- ysim
      }
      results[[gname[gg]]] <- list(fit=res_g, sim=sim_g)
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }
}
saveRDS(results, file = outfile)
