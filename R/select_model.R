#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param loo_results ELPD_loo for models
#' @param margin A multiplier for confidence interval
#' @return selected_model Reports the model scRATE selects (1:P, 2:NB, 3:ZIP, 4:ZINB)
#
select_model <- function(loo_results, margin=2) {

  m3idx <- which(rownames(loo_results) == 'model_fit$ZIP')
  if(length(m3idx)) { rownames(loo_results)[m3idx] <- 'model3' }
  m4idx <- which(rownames(loo_results) == 'model_fit$ZINB')
  if(length(m4idx)) { rownames(loo_results)[m4idx] <- 'model4' }

  if (rownames(loo_results)[1] == 'model1') {
    return(1)
  } else if (rownames(loo_results)[1] == 'model2') {
    if (abs(loo_results['model1',][['elpd_diff']]) < margin * loo_results['model1',][['se_diff']]) {
      return(1)
    } else {
      return(2)
    }
  } else if (rownames(loo_results)[1] == 'model3') {
    if (abs(loo_results['model1',][['elpd_diff']]) < margin * loo_results['model1',][['se_diff']]) {
      return(1)
    } else if (abs(loo_results['model2',][['elpd_diff']]) < margin * loo_results['model2',][['se_diff']]) {
      return(2)
    } else {
      return(3)
    }
  } else if (rownames(loo_results)[1] == 'model4') {
    if (abs(loo_results['model1',][['elpd_diff']]) < margin * loo_results['model1',][['se_diff']]) {
      return(1)
    } else if (abs(loo_results['model2',][['elpd_diff']]) < margin * loo_results['model2',][['se_diff']]) {
      return(2)
    } else if (abs(loo_results['model3',][['elpd_diff']]) < margin * loo_results['model3',][['se_diff']]) {
      return(3)
    } else {
      return(4)
    }
  }

}
