#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Five arguments must be supplied (<cntfile> <outfile> <formula_string> <cores> <seed>).n", call.=FALSE)
}

library(scRATE)

cntfile <- args[1]
outfile <- args[2]
formula_string <- args[3]
nCores <- as.integer(args[4])
seed <- as.integer(args[5])

load(cntfile)  # This will load 'cntmat', 'csize', and 'covar_list'

gname <- rownames(cntmat)
num_genes <- dim(cntmat)[1]
exposure <- log(csize)
message(sprintf("Formula: %s", formula_string))
f <- as.formula(formula_string)
covars2use <- all.vars(f)[-c(1, 2)]
if(identical(covars2use, character(0))) {
  covars2use <- NULL
}

for (gg in c(1:num_genes)) {
  y <- round(unlist(cntmat[gg,]))
  if(is.null(covars2use)) {
    gexpr <- data.frame(y, exposure)
  } else {
    gexpr <- data.frame(y, exposure, covar_list)
  }
  gsymb <- gname[gg]
  message(sprintf("\nFitting models for %s", gsymb))
  tryCatch({
    model_fit <- fit_count_models(gexpr, f, nCores, seed, adapt_delta = 0.95)
    elpd_loo <- compare_count_models(model_fit)
    mean_par <- get_model_params(model_fit, covariates=covars2use)
    res <- list()
    res[[gsymb]] <- list()
    res[[gsymb]][['orig']] <- list()
    res[[gsymb]][['orig']][['elpd_loo']] <- elpd_loo
    res[[gsymb]][['orig']][['mean_par']] <- mean_par
    saveRDS(model_fit, file=outfile)
    saveRDS(res, file=sprintf('%s.params.rds', outfile))
  }, error = function(err) {
    message(sprintf("Error while fitting %s", gsymb))
  })
}
