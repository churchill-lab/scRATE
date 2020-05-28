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
run_model_comparison <- function(cntfile, formula_string=NULL, nCores=NULL, seed=NULL, adapt_delta=0.8, outfile=NULL) {

  if(is.null(nCores)) {
    nCores <- parallel::detectCores()
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  load(cntfile)  # This will load 'cntmat', 'csize', and 'covar_list'

  gname <- rownames(cntmat)
  num_genes <- dim(cntmat)[1]
  exposure <- log(csize)
  if(is.null(formula_string)) {
    message(sprintf("Formulating the default additive model..."))
    formula_string <- 'y ~ 1 + offset(exposure)'
    covars <- names(covar_list)
    if(length(covars) > 0) {
      for (covar in covars) {
        formula_string <- paste(formula_string, sprintf(' + (1|%s)', covar))
      }
    }
  }
  message(sprintf("Formula: %s", formula_string))
  f <- as.formula(formula_string)
  covars2use <- all.vars(f)[-c(1, 2)]
  if(identical(covars2use, character(0))) {
    covars2use <- NULL
  }

  results <- list()
  for (gg in c(1:num_genes)) {
    y <- round(unlist(cntmat[gg,]))
    gexpr <- data.frame(y, exposure, covar_list)
    message(sprintf("\nFitting models for %s", gname[gg]))
    tryCatch({
      model_fit <- fit_count_models(gexpr, f, nCores, seed, adapt_delta = adapt_delta)
      elpd_loo <- compare_count_models(model_fit)
      mean_par <- get_model_params(model_fit, covariates=covars2use)
      results[[gname[gg]]] <- list()
      results[[gname[gg]]][['elpd_loo']] <- elpd_loo
      results[[gname[gg]]][['mean_par']] <- mean_par
    }, error = function(err) {
      message(sprintf("Error while fitting %s", gname[gg]))
    })
  }

  if(!is.null(outfile)) {
    saveRDS(results, file = outfile)
  }

  return(results)

}
