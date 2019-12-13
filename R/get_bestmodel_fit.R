#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param cntfile Expression quantity file (RData format: Use 'save' and 'load')
#' @param nCores Number of cores
#' @param seed Seed number
#' @param outfile Output file name to store ELPD_loo results (RDS format)
#' @return A list of model fits (only available when outfile is not specified)
#'
get_bestmodel_fit <- function(cntfile, nCores=NULL, seed=NULL, outfile=NULL) {

  if(is.null(nCores)) {
    nCores <- parallel::detectCores()
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  load(cntfile)  # This will load 'cntmat', 'gsurv', and 'csize'
  gname <- rownames(cntmat)
  num_genes <- length(gsurv)
  exposure <- log(csize)

  results <- list()
  for (gg in c(1:num_genes)) {
    if(gsurv[gg]) {
      y <- round(unlist(cntmat[gg,]))
      cat(sprintf("\nFitting models for %s\n", gname[gg]))
      tryCatch({
        results[[gname[gg]]] <- fit_count_models(y, exposure, nCores, seed, model2fit=model2fit[gg], brms4zi=TRUE)
      }, error = function(err) {
        cat(sprintf("Error while fitting %s\n", gname[gg]))
      })
    }
  }

  if(is.null(outfile)) {
    return(results)
  } else {
    saveRDS(results, file = outfile)
  }

}
