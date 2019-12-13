#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param model_fit A list of four model fits
#' @param ctyped Whether model_fit used ctype as covariate
#' @return A list of mean parameter values for each model
#
get_model_params <- function(model_fit, ctyped=FALSE) {
  result <- list()
  if(ctyped) {
    result[['P']] <- list()
    result[['P']][['summary']] <- as.data.frame(summary(model_fit[['P']]))
    result[['P']][['fixef']] <- fixef(model_fit[['P']])
    result[['P']][['ranef']] <- ranef(model_fit[['P']])$ctype
    result[['P']][['coef']] <- coef(model_fit[['P']])$ctype
    result[['NB']] <- list()
    result[['NB']][['summary']] <- as.data.frame(summary(model_fit[['NB']]))
    result[['NB']][['fixef']] <- fixef(model_fit[['NB']])
    result[['NB']][['ranef']] <- ranef(model_fit[['NB']])$ctype
    result[['NB']][['coef']] <- coef(model_fit[['NB']])$ctype
    result[['ZIP']] <- list()
    result[['ZIP']][['fixef']] <- summary(model_fit[['ZIP']])$fixed
    result[['ZIP']][['random']] <- summary(model_fit[['ZIP']])$random$ctype
    result[['ZIP']][['ranef']] <- ranef(model_fit[['ZIP']])$ctype[ , , 1]
    result[['ZIP']][['coef']] <- coef(model_fit[['ZIP']])$ctype[ , , 1]
    result[['ZINB']] <- list()
    result[['ZINB']][['fixef']] <- summary(model_fit[['ZINB']])$fixed
    result[['ZINB']][['random']] <- summary(model_fit[['ZINB']])$random$ctype
    result[['ZINB']][['shape']] <- summary(model_fit[['ZINB']])$spec_pars
    result[['ZINB']][['ranef']] <- ranef(model_fit[['ZINB']])$ctype[ , , 1]
    result[['ZINB']][['coef']] <- coef(model_fit[['ZINB']])$ctype[ , , 1]
  } else {
    result[['P']]    <- model_fit[['P']]$coefficients
    result[['NB']]   <- colMeans(as.data.frame(model_fit[['NB']]))
    result[['ZIP']]  <- colMeans(as.matrix(model_fit[['ZIP']], pars = c("b_Intercept", "b_zi_Intercept")))
    result[['ZINB']] <- colMeans(as.matrix(model_fit[['ZINB']], pars = c("b_Intercept", "b_zi_Intercept", "shape")))
  }
  return(result)
}
