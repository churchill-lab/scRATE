#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param model_fit A list of four model fits
#' @return ELPD_loo result
#
compare_count_models <- function(model_fit) {
  models <- names(model_fit)
  loo_list <- list()
  if('P' %in% models) {
    loo_1 <- loo(model_fit$P)
    loo_list[['P']] <- loo_1
  }
  if('NB' %in% models) {
    loo_2 <- loo(model_fit$NB)
    loo_list[['NB']] <- loo_2
  }
  if('ZIP' %in% models) {
    loo_3 <- loo(model_fit$ZIP)
    loo_list[['ZIP']] <- loo_3
  }
  if('ZINB' %in% models) {
    loo_4 <- loo(model_fit$ZINB)
    loo_list[['ZINB']] <- loo_4
  }
  res <- do.call(loo_compare, list(loo_list))
  return(res)
}
