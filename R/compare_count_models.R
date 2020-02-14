#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param model_fit A list of four model fits
#' @return ELPD_loo result
#
compare_count_models <- function(model_fit) {
  models <- names(model_fit)
  if('P' %in% models) {
    loo_1 <- loo(model_fit$P)
  } else if('NB' %in% models) {
    loo_2 <- loo(model_fit$NB)
  } else if('ZIP' %in% models) {
    loo_3 <- loo(model_fit$ZIP)
  } else if('ZINB' %in% models) {
    loo_4 <- loo(model_fit$ZINB)
  }
  res <- loo_compare(loo_1, loo_2, loo_3, loo_4)
  return(res)
}
