#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param cntfile Expression quantity file (RData format: Use 'save' and 'load')
#' @param nCores Number of cores
#' @param seed Seed number
#' @param adapt_delta The target average proposal acceptance probability during Stan’s adaptation period (default:0.8)
#' @param brms4zi Whether to run brms for zero-inflated models (default is to run rstan::sampling)
#' @param outfile Output file name to store ELPD_loo results (RDS format)
#' @return A list of ELPD_loo results returned by loo::compare
#'
run_model_comparison <- function(cntfile, nCores=NULL, seed=NULL, adapt_delta=0.8, brms4zi=FALSE, outfile=NULL) {

  if(is.null(nCores)) {
    nCores <- parallel::detectCores()
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  load(cntfile)  # This will load 'cntmat', 'gsurv', 'csize', and 'ctype'
  gname <- rownames(cntmat)
  num_genes <- length(gsurv)
  exposure <- log(csize)

  results <- list()
  for (gg in c(1:num_genes)) {
    if(gsurv[gg]) {
      y <- round(unlist(cntmat[gg,]))
      cat(sprintf("\nFitting models for %s\n", gname[gg]))
      tryCatch({
        model_fit <- fit_count_models(y, exposure, ctype, nCores, seed, adapt_delta = adapt_delta, brms4zi=brms4zi)
        elpd_loo <- compare_count_models(model_fit)
        mean_par <- get_model_params(model_fit, ctyped=!is.null(ctype))
        results[[gname[gg]]] <- list()
        results[[gname[gg]]][['elpd_loo']] <- elpd_loo
        results[[gname[gg]]][['mean_par']] <- mean_par
      }, error = function(err) {
        cat(sprintf("Error while fitting %s\n", gname[gg]))
      })
    }
  }

  if(!is.null(outfile)) {
    saveRDS(results, file = outfile)
  }

  return(results)

}
