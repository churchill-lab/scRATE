#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied (<indir> <outfile>)", call.=FALSE)
}

indir <- args[1]
outfile <- args[2]

flist <- Sys.glob(sprintf('%s/_scrate*', indir))
flist <- sort(flist, decreasing = FALSE)
results <- list()
for (k in 1:length(flist)) {
  newres <- readRDS(flist[k])
  results <- c(results, newres)
}
saveRDS(results, file = outfile)
bestmodel <- c(0, 0, 0, 0)
for (res in results) {
  if (rownames(res)[1] == 'model1') {
    bestmodel[1] = bestmodel[1] + 1
  } else if (rownames(res)[1] == 'model2') {
    bestmodel[2] = bestmodel[2] + 1
  } else if (rownames(res)[1] == 'model3') {
    bestmodel[3] = bestmodel[3] + 1
  } else if (rownames(res)[1] == 'model4') {
    bestmodel[4] = bestmodel[4] + 1
  }
}
cat(sprintf('%5d Poisson genes\n', bestmodel[1]))
cat(sprintf('%5d Negative Binomial genes\n', bestmodel[2]))
cat(sprintf('%5d Zero-Inflated Poisson genes\n', bestmodel[3]))
cat(sprintf('%5d Zero-Inflated Neg. Binomial genes\n', bestmodel[4]))
cat(sprintf('%5d Total\n', sum(bestmodel)))
