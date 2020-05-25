#' Fit a gene with regression models
#'
#' @export
#' @param gexpr A data frame that containes UMI count (y), log_exposure (exposure), and covariate values
#' @param f A regression formula to fit the non-zero-inflated counts
#' @param nCores Number of cores
#' @param seed Seed number
#' @param adapt_delta  The target average proposal acceptance probability during Stan’s adaptation period (default:0.8)
#' @param model2fit A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#' @param brms4zi Whether to run brms for zero-inflated models (Deprecated)
#' @return A list of `stanfit` returned by model(s)
#'
fit_count_models <- function(gexpr, f, nCores=NULL, seed=NULL, adapt_delta=0.8, model2fit=NULL, brms4zi=TRUE) {

  if(is.null(nCores)) {
    nCores <- min(4, parallel::detectCores())
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  # if(is.null(ctype)) {
  #   gexpr <- data.frame(y, exposure)
  # } else {
  #   gexpr <- data.frame(y, exposure, ctype)
  # }

  fitting <- list()

  if(is.null(model2fit) || 1 %in% model2fit) {
    message('Fitting data with Poisson model...')
    fitting[["P"]] <- stan_glmer(f,
                                 family = poisson,
                                 offset = exposure,
                                 data = gexpr,
                                 cores = nCores,
                                 seed = seed,
                                 refresh = 0)
    # if(is.null(ctype)) {
    #   fitting[["P"]] <- stan_glm(y ~ 1,
    #                              family = poisson,
    #                              offset = exposure,
    #                              data = gexpr,
    #                              cores = nCores,
    #                              seed = seed,
    #                              refresh = 0)
    # } else {
    #   fitting[["P"]] <- stan_glmer(y ~ 1 + (1 | ctype),
    #                                family = poisson,
    #                                offset = exposure,
    #                                data = gexpr,
    #                                cores = nCores,
    #                                seed = seed,
    #                                refresh = 0)
    # }
  }

  if(is.null(model2fit) || 2 %in% model2fit) {
    message('Fitting data with Negative Binomial model...')
    fitting[["NB"]] <- stan_glmer(f,
                                  family = neg_binomial_2,
                                  offset = exposure,
                                  data = gexpr,
                                  cores = nCores,
                                  seed = seed,
                                  refresh = 0)
    # if(is.null(ctype)) {
    #   fitting[["NB"]] <- stan_glm(y ~ 1,
    #                               family = neg_binomial_2,
    #                               offset = exposure,
    #                               data = gexpr,
    #                               cores = nCores,
    #                               seed = seed,
    #                               refresh = 0)
    # } else {
    #   fitting[["NB"]] <- stan_glmer(y ~ 1 + (1 | ctype),
    #                                 family = neg_binomial_2,
    #                                 offset = exposure,
    #                                 data = gexpr,
    #                                 cores = nCores,
    #                                 seed = seed,
    #                                 refresh = 0)
    # }
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
    # if(is.null(ctype)) {
    #   myprior_3 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1),
    #                          family = zero_inflated_poisson(),
    #                          data = gexpr)
    #   myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
    # } else {
    #   myprior_3 <- get_prior(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
    #                          family = zero_inflated_poisson(),
    #                          data = gexpr)
    #   myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
    # }
    # if(brms4zi) {
    #   if(is.null(ctype)) {
    #     fit_3 <- brm(bf(y ~ 1 + offset(exposure), zi ~ 1),
    #                  family = zero_inflated_poisson(),
    #                  data = gexpr,
    #                  prior = myprior_3,
    #                  control = list(adapt_delta = adapt_delta),
    #                  cores = nCores,
    #                  seed = seed,
    #                  refresh = 500)
    #   } else {
    #     fit_3 <- brm(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
    #                  family = zero_inflated_poisson(),
    #                  data = gexpr,
    #                  prior = myprior_3,
    #                  control = list(adapt_delta = adapt_delta),
    #                  cores = nCores,
    #                  seed = seed,
    #                  refresh = 500)
    #   }
    # } else {
    #   fit_3 <- rstan::sampling(stanmodels$zip,
    #                            data=list(N = length(y),
    #                                      Y = y,
    #                                      offset = exposure,
    #                                      prior_only = 0,
    #                                      df = myprior_3_values[1],
    #                                      loc = myprior_3_values[2],
    #                                      scale = myprior_3_values[3],
    #                                      control = list(adapt_delta = adapt_delta)),
    #                            cores = nCores,
    #                            seed = seed)
    # }
    # fitting[["ZIP"]] <- fit_3
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
    # if(is.null(ctype)) {
    #   myprior_4 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1),
    #                          family = zero_inflated_negbinomial(),
    #                          data = gexpr)
    #   myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
    # } else {
    #   myprior_4 <- get_prior(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
    #                          family = zero_inflated_negbinomial(),
    #                          data = gexpr)
    #   myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
    # }
    # if (brms4zi) {
    #   if(is.null(ctype)) {
    #     fit_4 <- brm(bf(y ~ 1 + offset(exposure), zi ~ 1),
    #                  family = zero_inflated_negbinomial(),
    #                  data = gexpr,
    #                  control = list(adapt_delta = adapt_delta),
    #                  prior = myprior_4,
    #                  cores = nCores,
    #                  seed = seed,
    #                  refresh = 500)
    #   } else {
    #     fit_4 <- brm(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
    #                  family = zero_inflated_negbinomial(),
    #                  data = gexpr,
    #                  control = list(adapt_delta = adapt_delta),
    #                  prior = myprior_4,
    #                  cores = nCores,
    #                  seed = seed,
    #                  refresh = 500)
    #   }
    # } else {
    #   fit_4 <- rstan::sampling(stanmodels$zinb,
    #                            data=list(N = length(y),
    #                                      Y = y,
    #                                      offset = exposure,
    #                                      prior_only = 0,
    #                                      df = myprior_4_values[1],
    #                                      loc = myprior_4_values[2],
    #                                      scale = myprior_4_values[3],
    #                                      control = list(adapt_delta = adapt_delta)),
    #                            cores = nCores,
    #                            seed = seed)
    # }
    # fitting[["ZINB"]] <- fit_4
  }

  return(fitting)

}
