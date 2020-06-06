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

models <- c('P', 'NB', 'ZIP', 'ZINB')
results <- list()
simcnts <- list()

for (gg in c(1:num_genes)) {
  if(gsurv[gg]) {
    y <- round(unlist(cntmat[gg,]))
    cat(sprintf("\nFitting models for %s\n", gname[gg]))
    tryCatch({
      model_fit <- fit_count_models(y, exposure, ctype, nCores, seed, adapt_delta=0.95, brms4zi=TRUE)
      elpd_loo <- compare_count_models(model_fit)
      mean_par <- get_model_params(model_fit, ctyped = !is.null(ctype))
      res <- list()
      res[['orig']] <- list()
      res[['orig']][['elpd_loo']] <- elpd_loo
      res[['orig']][['mean_par']] <- mean_par
      ysim_g <- list()
      for (m in 1:4) {
        # Generate simulated counts for each model
        cat(sprintf('Simulation by %s\n', models[m]))
        flush.console()
        yrep <- posterior_predict(model_fit[[m]], draws=1)
        ysim <- as.numeric(tail(yrep, 1))
        # Fit simulated counts with all four models
        model_fit_ysim <- fit_count_models(ysim, exposure, ctype, nCores, seed, adapt_delta=0.95, brms4zi=TRUE)
        res[[paste('sim', m, sep = '')]][['elpd_loo']] <- compare_count_models(model_fit_ysim)
        res[[paste('sim', m, sep = '')]][['mean_par']] <- get_model_params(model_fit_ysim, ctyped = !is.null(ctype))
        ysim_g[[paste('sim', m, sep = '')]] <- ysim
      }
      results[[gname[gg]]] <- res
      simcnts[[gname[gg]]] <- ysim_g 
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }
}

things2save <- list(modelfit=results, simcount=simcnts)
saveRDS(things2save, file = outfile)
