#' Fit regression models for many genes and save it
#'
#' @export
#' @param cntfile Expression quantity file (RData format: Use 'save' and 'load')
#' @param formula_string A regression formula to fit the non-zero-inflated counts
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

  load(cntfile)  # This will load 'cntmat', 'csize', 'covar_list', and 'model2fit'

  gname <- rownames(cntmat)
  num_genes <- dim(cntmat)[1]
  exposure <- log(csize)
  covars <- names(covar_list)
  if(is.null(formula_string)) {
    cat(sprintf("Formulating the default additive model...\n"))
    formula_string <- 'y ~ 1 + offset(exposure)'
    for (covar in covars) {
      formula_string <- paste(formula_string, sprintf(' + (1|%s)', covar))
    }
  }
  cat(sprintf("Formula: %s\n", formula_string))

  results <- list()
  for (gg in c(1:num_genes)) {
    y <- round(unlist(cntmat[gg,]))
    gexpr <- data.frame(y, exposure, covar_list)
    cat(sprintf("\nFitting models for %s\n", gname[gg]))
    tryCatch({
      results[[gname[gg]]] <- fit_count_models(gexpr, as.formula(formula_string), nCores, seed, adapt_delta=adapt_delta, model2fit=model2fit)
    }, error = function(err) {
      cat(sprintf("Error while fitting %s\n", gname[gg]))
    })
  }

  if(is.null(outfile)) {
    return(results)
  } else {
    saveRDS(results, file = outfile)
  }

}
