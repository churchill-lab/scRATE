#' Bayesian model selection for scRNA-seq count data
#'
#' @export
#' @param y Numeric vector of UMI counts
#' @param exposure Numeric vector of cell sizes (total UMI counts per cell)
#' @param ctype Cell types
#' @param adapt_delta  The target average proposal acceptance probability during Stan’s adaptation period (default:0.8)
#' @param nCores Number of cores
#' @param seed Seed number
#' @param model2fit A specific model to fit (1:P, 2:NB, 3:ZIP, 4:ZINB, NULL:All)
#' @param brms4zi Whether to run brms for zero-inflated models (default is to run rstan::sampling)
#' @return A list of `stanfit` returned by model(s)
#'
fit_count_models <- function(y, exposure, ctype=NULL, nCores=NULL, seed=NULL, adapt_delta=0.8, model2fit=NULL, brms4zi=FALSE) {

  if(is.null(nCores)) {
    nCores <- min(4, parallel::detectCores())
  }
  if(is.null(seed)) {
    seed <- 1004
  }

  if(is.null(ctype)) {
    gexpr <- data.frame(y, exposure)
  } else {
    gexpr <- data.frame(y, exposure, ctype)
  }
  fitting <- list()

  if(is.null(model2fit) || 1 %in% model2fit) {
    message('Fitting data with Poisson model...')
    if(is.null(ctype)) {
      fitting[["P"]] <- stan_glm(y ~ 1,
                                 family = poisson,
                                 offset = exposure,
                                 data = gexpr,
                                 cores = nCores,
                                 seed = seed,
                                 refresh = 0)
    } else {
      fitting[["P"]] <- stan_glmer(y ~ 1 + (1 | ctype),
                                   family = poisson,
                                   offset = exposure,
                                   data = gexpr,
                                   cores = nCores,
                                   seed = seed,
                                   refresh = 0)
    }
  }

  if(is.null(model2fit) || 2 %in% model2fit) {
    message('Fitting data with Negative Binomial model...')
    if(is.null(ctype)) {
      fitting[["NB"]] <- stan_glm(y ~ 1,
                                  family = neg_binomial_2,
                                  offset = exposure,
                                  data = gexpr,
                                  cores = nCores,
                                  seed = seed,
                                  refresh = 0)
    } else {
      fitting[["NB"]] <- stan_glmer(y ~ 1 + (1 | ctype),
                                    family = neg_binomial_2,
                                    offset = exposure,
                                    data = gexpr,
                                    cores = nCores,
                                    seed = seed,
                                    refresh = 0)
    }
  }

  if(is.null(model2fit) || 3 %in% model2fit) {
    message('Fitting data with Zero-Inflated Poisson model...')
    if(is.null(ctype)) {
      myprior_3 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1),
                             family = zero_inflated_poisson(),
                             data = gexpr)
      myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
    } else {
      myprior_3 <- get_prior(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
                             family = zero_inflated_poisson(),
                             data = gexpr)
      myprior_3_values <- eval(parse(text=gsub("student_t", "c", myprior_3$prior[1])))
    }
    if(brms4zi) {
      if(is.null(ctype)) {
        fit_3 <- brm(bf(y ~ 1 + offset(exposure), zi ~ 1),
                     family = zero_inflated_poisson(),
                     data = gexpr,
                     prior = myprior_3,
                     control = list(adapt_delta = adapt_delta),
                     cores = nCores,
                     seed = seed,
                     refresh = 500)
      } else {
        fit_3 <- brm(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
                     family = zero_inflated_poisson(),
                     data = gexpr,
                     prior = myprior_3,
                     control = list(adapt_delta = adapt_delta),
                     cores = nCores,
                     seed = seed,
                     refresh = 500)
      }
    } else {
      fit_3 <- rstan::sampling(stanmodels$zip,
                               data=list(N = length(y),
                                         Y = y,
                                         offset = exposure,
                                         prior_only = 0,
                                         df = myprior_3_values[1],
                                         loc = myprior_3_values[2],
                                         scale = myprior_3_values[3],
                                         control = list(adapt_delta = adapt_delta)),
                               cores = nCores,
                               seed = seed)
    }
    fitting[["ZIP"]] <- fit_3
  }

  if(is.null(model2fit) || 4 %in% model2fit) {
    message('Fitting data with Zero-Inflated Negative Binomial model...')
    if(is.null(ctype)) {
      myprior_4 <- get_prior(bf(y ~ 1 + offset(exposure), zi ~ 1),
                             family = zero_inflated_negbinomial(),
                             data = gexpr)
      myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
    } else {
      myprior_4 <- get_prior(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
                             family = zero_inflated_negbinomial(),
                             data = gexpr)
      myprior_4_values <- eval(parse(text=gsub("student_t", "c", myprior_4$prior[1])))
    }
    if (brms4zi) {
      if(is.null(ctype)) {
        fit_4 <- brm(bf(y ~ 1 + offset(exposure), zi ~ 1),
                     family = zero_inflated_negbinomial(),
                     data = gexpr,
                     control = list(adapt_delta = adapt_delta),
                     prior = myprior_4,
                     cores = nCores,
                     seed = seed,
                     refresh = 500)
      } else {
        fit_4 <- brm(bf(y ~ 1 + offset(exposure) + (1 | ctype), zi ~ 1),
                     family = zero_inflated_negbinomial(),
                     data = gexpr,
                     control = list(adapt_delta = adapt_delta),
                     prior = myprior_4,
                     cores = nCores,
                     seed = seed,
                     refresh = 500)
      }
    } else {
      fit_4 <- rstan::sampling(stanmodels$zinb,
                               data=list(N = length(y),
                                         Y = y,
                                         offset = exposure,
                                         prior_only = 0,
                                         df = myprior_4_values[1],
                                         loc = myprior_4_values[2],
                                         scale = myprior_4_values[3],
                                         control = list(adapt_delta = adapt_delta)),
                               cores = nCores,
                               seed = seed)
    }
    fitting[["ZINB"]] <- fit_4
  }

  message('Finished fitting all four models (P, NB, ZIP, & ZINB).')
  return(fitting)

}
