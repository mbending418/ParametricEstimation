library(glue)

#' @title maximum.likelihood.estimator
#' @description
#' uses maximum likelihood to estimate the parameters of a distribution from data
#' 
#' @name maximum.likelihood.estimator
#' 
#' @param data the data to estimate parameters from
#' 
#' @param distribution the distribution to estimate parameters for
#' 
#' @param lower.bound the lower bound for numerical estimation for parameters
#'  with no closed form
#' @param upper.bound the upper bound for numerical estimation for parameters
#'  with no closed form
#' @param tol the desired accuracy for numerical estimation for parameters
#'  with no closed form
#' 
#' @return returns a named numeric with the parameters
#' 
#' @export
maximum.likelihood.estimator <- function(data, 
                                         distribution = c("exponential",
                                                          "weibull",
                                                          "lognormal"),
                                         lower.bound = 0,
                                         upper.bound = 100,
                                         tol = 1e-10) {
  if (distribution == "exponential") {
    mu.hat = mean(data)
    par = setNames(c(mu.hat), c("mu"))
  } else if (distribution == "weibull") {
    n = length(data)
    loglik=function(beta.hat) sum(data^beta.hat * log(data))/sum(data^beta.hat) - sum(log(data))/n - 1/beta.hat
    beta.hat = uniroot(loglik, lower=lower.bound, upper=upper.bound, tol=tol)$root
    alpha.hat = (sum(data^beta.hat)/n)^(1/beta.hat)
    par = setNames(c(alpha.hat, beta.hat), c("alpha", "beta"))
  } else if (distribution == "lognormal") {
    mu.hat = mean(log(data))
    sigma.hat = sqrt(mean((log(data)-mu.hat)^2))
    setNames(c(mu.hat, sigma.hat), c("mu", "sigma"))
  } else {
    stop(glue("distribution unknown: {distribution}"))
  }
  par
}

