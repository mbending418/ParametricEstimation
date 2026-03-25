library(glue)

#' @title mom.estimator
#' @description
#' uses method of moments to estimate the parameters of a distribution from data
#' 
#' @name mom.estimator
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
mom.estimator <- function(data, 
                          distribution = c("exponential",
                                           "shifted.exponential",
                                            "weibull"),
                          lower.bound = 1,
                          upper.bound = 100,
                          tol = 1e-2) {
  moment <- function(k) mean(data^k)
  
  if (distribution == "exponential") {
    mu.hat = moment(1)
    par = setNames(c(mu.hat), c("mu"))
  } else if (distribution == "shifted.exponential"){
    mu.hat = sqrt(moment(2) - moment(1)^2)
    gamma.hat = moment(1) - mu.hat
    par = setNames(c(mu.hat, gamma.hat), c("mu", "gamma"))
  } else if (distribution == "weibull") {
    objective <- function(beta.hat) (gamma(2/beta.hat)/gamma(1/beta.hat+1)^2) - (moment(2)/moment(1)^2)
    beta.hat = uniroot(objective, lower=lower.bound, upper=upper.bound, tol=tol)$root
    alpha.hat = moment(1)/gamma(1/beta.hat + 1)
    par = setNames(c(alpha.hat, beta.hat), c("alpha", "beta"))
  } else {
    stop(glue("distribution unknown: {distribution}"))
  }
  par
}
