#' Fit a gene with regression models
#'
#' @export
#' @param gexpr A data frame that containes UMI count (y), log_exposure (exposure), and covariate values
#' @param f A regression formula to fit the non-zero-inflated counts
#' @param nCores Number of cores
#' @param seed Seed number
#' @param adapt_delta  The target average proposal acceptance probability during Stan’s adaptation period (default:0.8)
#' @param model2fit A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#' @return A list of `stanfit` returned by model(s)
#'
fit_count_models <- function(gexpr, f, nCores=NULL, seed=NULL, adapt_delta=0.8, model2fit=NULL) {

  if(is.null(nCores)) {
    nCores <- min(4, parallel::detectCores())
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  covariates <- all.vars(f)[-c(1, 2)]
  fitting <- list()

  if(is.null(model2fit) || 1 %in% model2fit) {
    message('Fitting data with Poisson model...')
    if(identical(covariates, character(0))) {
      fitting[["P"]] <-   stan_glm(f,
                                   family = poisson,
                                   data = gexpr,
                                   cores = nCores,
                                   seed = seed,
                                   refresh = 0)
    } else {
      fitting[["P"]] <- stan_glmer(f,
                                   family = poisson,
                                   data = gexpr,
                                   cores = nCores,
                                   seed = seed,
                                   refresh = 0)
    }
  }

  if(is.null(model2fit) || 2 %in% model2fit) {
    message('Fitting data with Negative Binomial model...')
    if(identical(covariates, character(0))) {
      fitting[["NB"]] <-   stan_glm(f,
                                    family = neg_binomial_2,
                                    data = gexpr,
                                    cores = nCores,
                                    seed = seed,
                                    refresh = 0)
    } else {
      fitting[["NB"]] <- stan_glmer(f,
                                    family = neg_binomial_2,
                                    data = gexpr,
                                    cores = nCores,
                                    seed = seed,
                                    refresh = 0)
    }
  }

  if(is.null(model2fit) || 3 %in% model2fit) {
    message('Fitting data with Zero-Inflated Poisson model...')
    myprior_3 <- get_prior(bf(f, zi ~ 1),
                           family = zero_inflated_poisson(),
                           data = gexpr)
    myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
    fitting[["ZIP"]] <- brm(bf(f, zi ~ 1),
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
    myprior_4 <- get_prior(bf(f, zi ~ 1),
                           family = zero_inflated_negbinomial(),
                           data = gexpr)
    myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
    fitting[["ZINB"]] <- brm(bf(f, zi ~ 1),
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
