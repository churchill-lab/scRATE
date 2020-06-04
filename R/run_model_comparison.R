#' Fit regression models for many genes
#'
#' @export
#' @param cntfile Expression quantity file (RData format: Use 'save' and 'load'). The count matrix is gene x cell.
#' @param formula_string A regression formula to fit the non-zero-inflated counts
#' @param nCores Number of cores
#' @param seed Seed number
#' @param adapt_delta The target average proposal acceptance probability during Stan’s adaptation period (default:0.8)
#' @param outfile Output file name to store ELPD_loo results (RDS format)
#' @return A list of ELPD_loo results and mean parameters returned by loo::compare
#'
run_model_comparison <- function(cntfile, formula_string=NULL, nCores=NULL, seed=NULL, adapt_delta=0.8, model2fit=NULL, outfile=NULL) {

  if(is.null(nCores)) {
    nCores <- parallel::detectCores()
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  # We assume 'cntfile' contains 'cntmat', 'csize', and 'covar_list' (See 'prepare_job_array' function)
  load(cntfile)

  gname <- rownames(cntmat)
  num_genes <- dim(cntmat)[1]
  exposure <- log(csize)
  covars <- names(covar_list)

  results <- list()
  for (gg in c(1:num_genes)) {
    y <- round(unlist(cntmat[gg,]))
    if(is.null(covars)) {
      gexpr <- data.frame(y, exposure)
    } else {
      gexpr <- data.frame(y, exposure, covar_list)
    }
    gsymb <- gname[gg]
    message(sprintf("\nFitting models for %s", gsymb))
    tryCatch({
      model_fit <- fit_count_models(gexpr, formula_string, nCores, seed, adapt_delta = adapt_delta, model2fit = model2fit)
      elpd_loo <- compare_count_models(model_fit)
      mean_par <- get_model_params(model_fit, covariates=covars)
      results[[gsymb]] <- list()
      results[[gsymb]][['elpd_loo']] <- elpd_loo
      results[[gsymb]][['mean_par']] <- mean_par
    }, error = function(err) {
      message(sprintf("Error while fitting %s", gsymb))
    })
  }

  if(!is.null(outfile)) {
    saveRDS(results, file = outfile)
  }

  return(results)

}
