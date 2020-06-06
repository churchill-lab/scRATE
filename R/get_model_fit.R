#' Fit regression models for many genes and save the fitting objects
#'
#' @export
#' @param cntfile Expression quantity file (RData format: Use 'save' and 'load')
#' @param formula_string A regression formula to fit the non-zero-inflated counts
#' @param model2fit A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#' @param nCores Number of cores
#' @param seed Seed number
#' @param outfile Output file name to store ELPD_loo results (RDS format)
#' @return A list of model fits (only available when outfile is not specified)
#'
get_model_fit <- function(cntfile, formula_string=NULL, model2fit=NULL, nCores=NULL, seed=NULL, outfile=NULL) {

  if(is.null(nCores)) {
    nCores <- parallel::detectCores()
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  # We assume 'cntfile' contains 'cntmat', 'csize', 'covar_list', and model2fit (See 'prepare_job_array' function)
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
      results[[gsymb]] <- fit_count_models(gexpr, formula_string, nCores, seed, adapt_delta = adapt_delta, model2fit = model2fit)
      if(!is.null(outfile)) {
        if(file.exists(outfile)) {
          reslist <- readRDS(outfile)
          reslilst[[gsymb]] <-results[[gsymb]]
          saveRDS(reslist, file = outfile)
        } else {
          saveRDS(results, file = outfile)
        }
      }
    }, error = function(err) {
      message(sprintf("Error while fitting %s", gsymb))
    })
  }

  if(is.null(outfile)) {
    return(results)
  }

}
