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

load(cntfile)  # This will load 'cntmat', 'gsurv', and 'csize'
gname <- rownames(cntmat)
num_genes <- length(gsurv)
exposure <- log(csize)

models <- c('P', 'NB', 'ZIP', 'ZINB')
results <- list()

for (gg in c(1:num_genes)) {
  if(gsurv[gg]) {
    y <- round(unlist(cntmat[gg,]))
    cat(sprintf("\nTesting %s\n", gname[gg]))
    tryCatch({
      model_fit <- fit_count_models(y, exposure, nCores, seed, brms4zi = TRUE)
      elpd_loo <- compare_count_models(model_fit)
      res <- list(orig=elpd_loo)
      for (m in 1:4) {
        # Generate simulated counts for each model
        cat(sprintf('Simulation by %s\n', models[m]))
        flush.console()
        yrep <- posterior_predict(model_fit[[m]], draws=1)
        ysim <- as.numeric(tail(yrep, 1))
        # Fit simulated counts with all four models
        model_fit_ysim <- fit_count_models(ysim, exposure, nCores, seed)
        res[[paste('sim', m, sep = '')]] <- compare_count_models(model_fit_ysim)
      }
      results[[gname[gg]]] <- res
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }
}

saveRDS(results, file = outfile)
