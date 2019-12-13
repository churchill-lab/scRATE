#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loo_dir A folder name in which leave-one-out ELPD result files reside
#' @param globstr Search string (wildcard supported) for loo result files (in RDS format)
#' @param loo_outfile Name of the file to save collated leave-one-out ELPD results (optional)
#' @return A list of collated results if 'loo_outfile' is not given
#'
collate_results <- function(loo_dir, globstr='_scrate*', loo_outfile=NULL) {
                                     
  flist <- Sys.glob(file.path(loo_dir, globstr))
  flist <- sort(flist, decreasing = FALSE)
  results <- c()
  for (i in 1:length(flist)) {
    results <- c(results, readRDS(flist[i]))
  }

  # Remove genes that failed to fit for some reason, if any
  gsurv <- names(results)
  for (g in gsurv) {
    if(length(results[[g]]) == 0) {
      cat(sprintf('Gene %s failed to fit for some reason. Check the log.\n', g))
      flush.console()
      results[[g]] <- NULL  # Remove gene name
    }
  }

  # Fix names in loo results if needed
  for (k in 1:length(results)) {
    loo_results <- results[[k]][['elpd_loo']]
    m3idx <- which(rownames(loo_results) == 'model_fit$ZIP')
    if(length(m3idx)) { rownames(results[[k]][['elpd_loo']])[m3idx] <- 'model3' }
    m4idx <- which(rownames(loo_results) == 'model_fit$ZINB')
    if(length(m4idx)) { rownames(results[[k]][['elpd_loo']])[m4idx] <- 'model4' }
  }

  if (!is.null(loo_outfile)) {
    saveRDS(results, file=loo_outfile)
  } else {
    return(results)
  }

}
