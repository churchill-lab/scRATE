#' An internal function to re-organize stan_fit parameters
#'
#' @export
#' @param model_fit A list of four model fits
#' @param covariates A list of covariates if any
#' @return A list of mean parameter values for each model
#
get_model_params <- function(model_fit, covariates=NULL) {

  result <- list()
  models <- names(model_fit)
  if(is.null(covariates)) {
    if('P' %in% models) {
      result[['P']]    <- model_fit[['P']]$coefficients
    }
    if('NB' %in% models) {
      result[['NB']]   <- colMeans(as.data.frame(model_fit[['NB']]))
    }
    if('ZIP' %in% models) {
      result[['ZIP']]  <- colMeans(as.matrix(model_fit[['ZIP']], pars = c("b_Intercept", "b_zi_Intercept")))
    }
    if('ZINB' %in% models) {
      result[['ZINB']] <- colMeans(as.matrix(model_fit[['ZINB']], pars = c("b_Intercept", "b_zi_Intercept", "shape")))
    }
  } else {
    if('P' %in% models) {
      result[['P']] <- list()
      result[['P']][['summary']] <- as.data.frame(summary(model_fit[['P']]))
      result[['P']][['fixef']] <- fixef(model_fit[['P']])
      result[['P']][['ranef']] <- list()
      result[['P']][['coef']] <- list()
      for (covar in covariates) {
        result[['P']][['ranef']][[covar]] <- ranef(model_fit[['P']])[[covar]]
        result[['P']][['coef']][[covar]] <- coef(model_fit[['P']])[[covar]]
      }
    }
    if('NB' %in% models) {
      result[['NB']] <- list()
      result[['NB']][['summary']] <- as.data.frame(summary(model_fit[['NB']]))
      result[['NB']][['fixef']] <- fixef(model_fit[['NB']])
      result[['NB']][['ranef']] <- list()
      result[['NB']][['coef']] <- list()
      for (covar in covariates) {
        result[['NB']][['ranef']][[covar]] <- ranef(model_fit[['NB']])[[covar]]
        result[['NB']][['coef']][[covar]] <- coef(model_fit[['NB']])[[covar]]
      }
    }
    if('ZIP' %in% models) {
      result[['ZIP']] <- list()
      result[['ZIP']][['fixef']] <- summary(model_fit[['ZIP']])$fixed
      result[['ZIP']][['random']] <- list()
      result[['ZIP']][['ranef']] <- list()
      result[['ZIP']][['coef']] <- list()
      for (covar in covariates) {
        result[['ZIP']][['random']][[covar]] <- summary(model_fit[['ZIP']])$random[[covar]]
        result[['ZIP']][['ranef']][[covar]] <- ranef(model_fit[['ZIP']])[[covar]][ , , 1]
        result[['ZIP']][['coef']][[covar]] <- coef(model_fit[['ZIP']])[[covar]][ , , 1]
      }
    }
    if('ZINB' %in% models) {
      result[['ZINB']] <- list()
      result[['ZINB']][['fixef']] <- summary(model_fit[['ZINB']])$fixed
      result[['ZINB']][['shape']] <- summary(model_fit[['ZINB']])$spec_pars
      result[['ZINB']][['random']] <- list()
      result[['ZINB']][['ranef']] <- list()
      result[['ZINB']][['coef']] <- list()
      for (covar in covariates) {
        result[['ZINB']][['random']][[covar]] <- summary(model_fit[['ZINB']])$random[[covar]]
        result[['ZINB']][['ranef']][[covar]] <- ranef(model_fit[['ZINB']])[[covar]][ , , 1]
        result[['ZINB']][['coef']][[covar]] <- coef(model_fit[['ZINB']])[[covar]][ , , 1]
      }
    }
  }
  return(result)
}
