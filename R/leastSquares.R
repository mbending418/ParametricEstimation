library(glue)
library(VGAM)

#' @title least.squares.estimator
#' @description
#' uses least squares to estimate the parameters of a distribution from data
#' 
#' @name least.squares.estimator
#' 
#' @param data the data to estimate parameters from
#' 
#' @param censor a censorship indicator (an array of 0's and 1's the length of data with 1 indicating not censored)
#' 
#' @param distribution the distribution to estimate parameters for
#' 
#' @param plot.fit set to TRUE to plot the least squares fit
#' 
#' @return returns a named numeric with the parameters, MTTF, slope, intercept and R.squared
#' 
#' @export
least.squares.estimator <- function(data,
                                    censor = seq(1,n),
                                    distribution = c("exponential",
                                                     "shift.exponential",
                                                     "weibull",
                                                     "shift.weibull",
                                                     "lognormal",
                                                     "birnbaum.saunders"),
                                    plot.fit = FALSE) {
  n = length(data)
  data = sort(data)
  i = seq(1, n)
  
  if (distribution=="exponential" || distribution == "shift.exponential") {
    x.var = data
    y.var = -log(1-i/(n+1))
    x.lab = "xi"
    y.lab = "-log(1-(i/(n+1))))"
  } else if (distribution == "weibull") {
    x.var = log(data)
    y.var = log(-log(1 - (i/(n+1))))
    x.lab = "log(xi)"
    y.lab = "log(-log(1-(i/(n+1)))))"
  } else if (distribution == "shift.weibull") {
    objective <- function(gam) {
      unname(least.squares.estimator(data-gam, distribution = "weibull")["R.squared"])
    }
    gam = optimize(objective, interval = (c(0, min(data))), maximum=TRUE)$maximum
    x.var = log(data-gam)
    y.var = log(-log(1 - (i/(n+1))))
    x.lab = glue("log(xi-{round(gam, digits=4)})")
    y.lab = "log(-log(1-(i/(n+1)))))"
  } else if (distribution == "lognormal") {
    x.var = log(data)
    y.var = qnorm(i/(n+1))
    x.lab = "log(xi)"
    y.lab = "inv.erf(i/(n+1))"
  } else if (distribution == "birnbaum.saunders") {
    x.var = data
    y.var = qnorm(i/(n+1)) * sqrt(data)
    x.lab = "xi"
    y.lab = "inv.erf(i/(n+1)) * sqrt(xi)"
  } else {
    stop(glue("distribution unknown: {distribution}"))
  }
  
  #removed the censored data
  x.var = x.var[censor*seq(1,n)]
  y.var = y.var[censor*seq(1,n)]
  
  model = lm(formula=y.var~x.var)
  slope=unname(model$coefficients[2])
  intercept=unname(model$coefficients[1])
  if (plot.fit) {
    plot(x.var, y.var, xlab=x.lab, ylab=y.lab)
    abline(model)
    legend("topleft", legend = glue("r={round(cor(x.var, y.var), digits=4)} |  slope={round(slope, digits=4)} | intercept={round(intercept, digits=4)}"))
  }
  if (distribution=="exponential") {
    mu.tilde = 1/slope
    res = setNames(c(mu.tilde, 1/mu.tilde, mu.tilde),
                   c("mu", "lambda", "MTTF"))
  } else if (distribution == "shift.exponential")
  {
    mu.tilde = 1/slope
    gamma.tilde = -intercept / slope
    res = setNames(c(mu.tilde, 1/mu.tilde, gamma.tilde, mu.tilde+gamma.tilde),
                   c("mu", "lambda", "gamma", "MTTF"))
  } else if (distribution == "weibull") {
    alpha.tilde = exp(-intercept/slope)
    beta.tilde = slope
    res = setNames(c(alpha.tilde, beta.tilde, alpha.tilde*gamma(1/beta.tilde + 1)),
                   c("alpha", "beta", "MTTF"))
  } else if (distribution == "shift.weibull") {
    alpha.tilde = exp(-intercept/slope)
    beta.tilde = slope
    gamma.tilde = gam
    res = setNames(c(alpha.tilde, beta.tilde, gamma.tilde, alpha.tilde*gamma(1/beta.tilde + 1) + gamma.tilde),
                   c("alpha", "beta","gamma", "MTTF"))
  } else if (distribution == "lognormal") {
    mu.tilde = -intercept/slope
    sigma.tilde = 1/slope
    res = setNames(c(mu.tilde, sigma.tilde, exp(mu.tilde+.5*sigma.tilde^2)),
                   c("mu", "sigma", "MTTF"))
  } else if (distribution == "birnbaum.saunders") {
    alpha.tilde = 1/sqrt(-slope*intercept)
    beta.tilde = -intercept/slope
    res = setNames(c(alpha.tilde, beta.tilde),
                   c("alpha", "beta"))
  } else {
    stop(glue("distribution unknown: {distribution}"))
  }
  
  res["slope"] = slope
  res["intercept"] = intercept
  res["R.squared"] = cor(x.var,y.var)^2
  res
}