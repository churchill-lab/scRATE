#' Fit a gene with regression models
#'
#' @export
#' @param gexpr A data frame that containes UMI count (y), log_depth (exposure), and covariate values
#' @param formula_string A regression formula to fit the non-zero-inflated counts (Note: Do not use 'offset()')
#' @param nCores Number of cores
#' @param seed Seed number
#' @param adapt_delta  The target average proposal acceptance probability during Stan’s adaptation period (default:0.8)
#' @param model2fit A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#' @return A list of `stanfit` returned by model(s)
#'
fit_count_models <- function(gexpr, formula_string=NULL, nCores=NULL, seed=NULL, adapt_delta=0.8, model2fit=NULL) {

  if(is.null(nCores)) {
    nCores <- min(4, parallel::detectCores())
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  # Automatically formulate a simple additive model using all the covariates in 'gexpr'
  covariates <- names(gexpr)[-c(1, 2)]
  if(is.null(formula_string)) {
    message(sprintf("Formulating the default additive model..."))
    formula_string <- 'y ~ 1'
    if(!identical(covariates, character(0))) {
      for (covar in covariates) {
        formula_string <- paste(formula_string, sprintf(' + (1|%s)', covar))
      }
    }
  }
  message(sprintf("Formula: %s", formula_string))
  f12 <- as.formula(formula_string)
  f34 <- as.formula(paste(formula_string, ' + offset(exposure)'))

  fitting <- list()

  if(is.null(model2fit) || 1 %in% model2fit) {
    message('Fitting data with Poisson model...')
    if(identical(covariates, character(0))) {
      fitting[["P"]] <-   stan_glm(f12,
                                   family = poisson,
                                   data = gexpr,
                                   offset = exposure,
                                   cores = nCores,
                                   seed = seed,
                                   refresh = 0)
    } else {
      fitting[["P"]] <- stan_glmer(f12,
                                   family = poisson,
                                   data = gexpr,
                                   offset = exposure,
                                   cores = nCores,
                                   seed = seed,
                                   refresh = 0)
    }
  }

  if(is.null(model2fit) || 2 %in% model2fit) {
    message('Fitting data with Negative Binomial model...')
    if(identical(covariates, character(0))) {
      fitting[["NB"]] <-   stan_glm(f12,
                                    family = neg_binomial_2,
                                    data = gexpr,
                                    offset = exposure,
                                    cores = nCores,
                                    seed = seed,
                                    refresh = 0)
    } else {
      fitting[["NB"]] <- stan_glmer(f12,
                                    family = neg_binomial_2,
                                    data = gexpr,
                                    offset = exposure,
                                    cores = nCores,
                                    seed = seed,
                                    refresh = 0)
    }
  }

  if(is.null(model2fit) || 3 %in% model2fit) {
    message('Fitting data with Zero-Inflated Poisson model...')
    myprior_3 <- get_prior(bf(f34, zi ~ 1),
                           family = zero_inflated_poisson(),
                           data = gexpr)
    myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
    fitting[["ZIP"]] <- brm(bf(f34, zi ~ 1),
                            family = zero_inflated_poisson(),
                            data = gexpr,
                            prior = myprior_3,
                            control = list(adapt_delta = adapt_delta),
                            cores = nCores,
                            seed = seed,
                            refresh = 500)
  }

  if(is.null(model2fit) || 4 %in% model2fit) {
    message('Fitting data with Zero-Inflated Negative Binomial model...')
    myprior_4 <- get_prior(bf(f34, zi ~ 1),
                           family = zero_inflated_negbinomial(),
                           data = gexpr)
    myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
    fitting[["ZINB"]] <- brm(bf(f34, zi ~ 1),
                             family = zero_inflated_negbinomial(),
                             data = gexpr,
                             control = list(adapt_delta = adapt_delta),
                             prior = myprior_4,
                             cores = nCores,
                             seed = seed,
                             refresh = 500)
  }

  return(fitting)

}
