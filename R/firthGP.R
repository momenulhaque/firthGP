
#' Estimate generalized Poisson regression model for count data where the data has separation, quasi-separation or near-to-quasi separation problem
#'
#' @param formula similar as glm
#' @param data a data frame or tibble
#' @param alpha level of significance for confidence interval
#'
#' @return a table of coefficients and ci
#'
#' @importFrom stats model.extract model.frame model.matrix optim qchisq qt
#'
#' @export
#'
#'
#'
#' @examples
#'
#' quasifit = firthGP(y ~ x1 + x2, data = quasiSep)
#' quasifit
#'
#'
#'
#'
#'
firthGP <- function(formula, data = parent.frame(), alpha=0.05){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m[[1]] <- as.name("model.frame")
  m <- model.frame(formula, data)
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m)
  xt <- t(x)
  k <- ncol(x); n <- nrow(x)
  init <- rep(0, k+1)
  xnames <- dimnames(x)[[2]]
  if (is.null(xnames)) {
    xnames <- paste("x", 1:k, sep = "")
    dimnames(xx) <- list(NULL, xnames)
  }
  logLikfull <- function(beta, y, x){
    # p_{i}, given X_{i} and beta
    pi <- as.vector(exp(as.matrix(x) %*% c(beta[1:k])))
    logli <- (y-1) * log(pi+beta[k+1]*y) - log(exp(pi+beta[k+1]*y)) + log(pi)
    Vi = as.vector(1 - (beta[k+1]*y*(y-1))/(pi+beta[k+1]*y)^2)
    #mi = as.vector((y*(y-1))/(pi+beta[4]*y)^2)
    #I1 = c(sum(pi*mi), sum(as.vector(xt[2])*pi*mi), sum(as.vector(xt[3])*pi*mi), sum(as.vector(y*mi)))
    I <- xt %*% diag(pi*Vi) %*% x
    #I = cbind(rbind(I2, I1[1:3]), I1)
    penal <- .5 * determinant(I, logarithm = TRUE)$modulus[1]
    return(-sum(logli) - penal)
  }

  suppressWarnings(
    fit.full <- optim(par=rep(0, k+1), fn=logLikfull, y=y, x=x, method="BFGS", hessian=T)
  )



  #fit.null <- optim(par=0, fn=logLiknull, method="BFGS", hessian=T)
  LL.0 <- -fit.full$value - qchisq(1 - alpha, 1)/2
  fit <- list(coefficients = fit.full$par, sd = sqrt(diag(solve(fit.full$hessian))), alpha = alpha,  loglik =-fit.full$value,
              terms = colnames(x), formula = formula(formula), call=match.call())

  fit$ci.lower <- fit$ci.upper <- rep(0, k+1)
  lik <- NULL
  for(i in 1:(k+1)){
    lo <- fit$coefficients[i] - fit$sd[i] * qt(1-alpha/2, n-1)
    up <- fit$coefficients[i] + fit$sd[i] * qt(1-alpha/2, n-1)

    fit$ci.lower[i] <-  lo;
    fit$ci.upper[i] <-  up;

  }

  est.swM <- data.frame(Estimate=round(fit$coefficients, 3),
                        std.err = round(fit$sd, 3),
                        lower.ci=round(fit$ci.lower, 3),
                        upper.ci=round(fit$ci.upper, 3))


  row.names(est.swM) <- c(xnames, "delta")

  fit1 <- list()
  attr(fit1, "class") <- c("poison")
  fit1$call <- call
  fit1$coefficients <- est.swM[xnames,]
  fit1$dispersion <- est.swM["delta",]
  fit1
}



