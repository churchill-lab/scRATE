#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One arguments must be supplied (<dataset>).n", call.=FALSE)
}
outbase <- args[1]

library(scRATE)

flist <- Sys.glob('_scrate*')

cat("Collating the entire results by gene...\n")
collate_results(loo_dir='.', loo_outfile=sprintf('%s.sim_4x4.rds', outbase))
res <- readRDS(sprintf('%s.sim_4x4.rds', outbase))
gname <- names(res)

cat("Collating model parameters...\n")
model_fit <- list()
for (gg in gname) {
  model_fit[[gg]] <- res[[gg]][['fit']]
}
saveRDS(model_fit, file=sprintf('%s.sim_4x4.model_fit.rds', outbase))

num_genes <- length(gname)
num_cells <- length(res[[1]][['sim']][[1]])
sim_model <- names(res[[1]][['sim']])

if ('P' %in% sim_model) {
  cat("Collating P simulation...\n")
  sim <- matrix(nrow=num_genes, ncol=num_cells)
  for (k in 1:num_genes) {
    gg <- gname[[k]]
    sim[k,] <- res[[gg]][['sim']][['P']]
  }
  rownames(sim) <- gname
  saveRDS(sim, file=sprintf('%s.sim_4x4.simP.rds', outbase))
}

if ('NB' %in% sim_model) {
  cat("Collating NB simulation...\n")
  sim <- matrix(nrow=num_genes, ncol=num_cells)
  for (k in 1:num_genes) {
    gg <- gname[[k]]
    sim[k,] <- res[[gg]][['sim']][['NB']]
  }
  rownames(sim) <- gname
  saveRDS(sim, file=sprintf('%s.sim_4x4.simNB.rds', outbase))
}

if ('ZIP' %in% sim_model) {
  cat("Collating ZIP simulation...\n")
  sim <- matrix(nrow=num_genes, ncol=num_cells)
  for (k in 1:num_genes) {
    gg <- gname[[k]]
    sim[k,] <- res[[gg]][['sim']][['ZIP']]
  }
  rownames(sim) <- gname
  saveRDS(sim, file=sprintf('%s.sim_4x4.simZIP.rds', outbase))
}

if ('ZINB' %in% sim_model) {
  cat("Collating ZINB simulation...\n")
  sim <- matrix(nrow=num_genes, ncol=num_cells)
  for (k in 1:num_genes) {
    gg <- gname[[k]]
    sim[k,] <- res[[gg]][['sim']][['ZINB']]
  }
  rownames(sim) <- gname
  saveRDS(sim, file=sprintf('%s.sim_4x4.simZINB.rds', outbase))
}

cat("Finished collating results.\n")
