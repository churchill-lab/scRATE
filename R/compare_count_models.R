#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param model_fit A list of four model fits
#' @return ELPD_loo result
#
compare_count_models <- function(model_fit) {
  loo_1 <- loo(model_fit$P)
  loo_2 <- loo(model_fit$NB)
  loo_3 <- loo(model_fit$ZIP)
  loo_4 <- loo(model_fit$ZINB)
  res <- loo_compare(loo_1, loo_2, loo_3, loo_4)
  return(res)
}
