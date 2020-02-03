#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loo_results ELPD_loo for models
#' @param margin A multiplier for confidence interval
#' @return selected_model Reports the model scRATE selects (1:P, 2:NB, 3:ZIP, 4:ZINB)
#
select_model <- function(loo_results, margin=2) {

  m1idx <- which(rownames(loo_results) == 'model_fit$P'    | rownames(loo_results) == 'model1')
  if(length(m1idx)) { rownames(loo_results)[m1idx] <- 'P' }
  m2idx <- which(rownames(loo_results) == 'model_fit$NB'   | rownames(loo_results) == 'model2')
  if(length(m2idx)) { rownames(loo_results)[m2idx] <- 'NB' }
  m3idx <- which(rownames(loo_results) == 'model_fit$ZIP'  | rownames(loo_results) == 'model3')
  if(length(m3idx)) { rownames(loo_results)[m3idx] <- 'ZIP' }
  m4idx <- which(rownames(loo_results) == 'model_fit$ZINB' | rownames(loo_results) == 'model4')
  if(length(m4idx)) { rownames(loo_results)[m4idx] <- 'ZINB' }

  if (rownames(loo_results)[1] == 'P') {
    return(1)
  } else if (rownames(loo_results)[1] == 'NB') {
    if (abs(loo_results['P',][['elpd_diff']]) < margin * loo_results['P',][['se_diff']]) {
      return(1)
    } else {
      return(2)
    }
  } else if (rownames(loo_results)[1] == 'ZIP') {
    if (abs(loo_results['P',][['elpd_diff']]) < margin * loo_results['P',][['se_diff']]) {
      return(1)
    } else if (abs(loo_results['NB',][['elpd_diff']]) < margin * loo_results['NB',][['se_diff']]) {
      return(2)
    } else {
      return(3)
    }
  } else if (rownames(loo_results)[1] == 'ZINB') {
    if (abs(loo_results['P',][['elpd_diff']]) < margin * loo_results['P',][['se_diff']]) {
      return(1)
    } else if (abs(loo_results['NB',][['elpd_diff']]) < margin * loo_results['NB',][['se_diff']]) {
      return(2)
    } else if (abs(loo_results['ZIP',][['elpd_diff']]) < margin * loo_results['ZIP',][['se_diff']]) {
      return(3)
    } else {
      return(4)
    }
  }

}
