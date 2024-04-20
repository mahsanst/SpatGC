#' Gamma-Count (GC) Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the GC distribution with
#' parameters \eqn{\alpha} and \eqn{\gamma}.
#'
#' The GC distribution with parameters  \eqn{\alpha} and \eqn{\gamma}
#' has the density
#' \deqn{P(Y_T=y)=G(y\alpha,\gamma T)-G\left(\left(y+1\right)\alpha,\gamma T\right)}
#' where \deqn{G(n\alpha,\gamma T)=\frac{1}{\Gamma(n\alpha)}\int_{0}^{\gamma T} u^{n\alpha-1}\exp(-u)du}
#' for \eqn{\alpha} and \eqn{\lambda} which must be positive
#' values and \eqn{y \in \{0, 1, 2, \ldots\}}.
#'
#' @name GC
#' @param y a vector or matrix of observations for which the pdf needs to be computed.
#' @param q a vector or matrix of quantiles for which the cdf needs to be computed.
#' @param p a vector or matrix of probabilities for which the quantile needs to be computed.
#' @param n number of random values to be generated.
#' @param alpha the dispersion parameter of GC: default is 1; (shape parameter of waiting time variable);
#' @param gamma the rate parameter of GC and waiting time variable;
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(Y \ge y)}; otherwise, \eqn{P(Y > y)}.
#' @param method Character string. The method used for generating random variables from the GC distribution in `rGC`. Options are:
#'               - `"PF"`: Based on the probability function.
#'               - `"IT"`: Inverse transformation method based on the quantile function.
#'               - `"WT"`: Based on the waiting times distribution.
#'
#' @return
#'  \code{G} gives the G function,
#'  \code{pGC} gives the distribution function,
#'  \code{dGC} gives the density,
#'  \code{qGC} gives the quantile function and
#'  \code{rGC} generates random variables from the GC Distribution.
#' @references
#'      Winkelmann, R. (1995). Duration dependence and dispersion in count-data models.
#'          \emph{Journal of Business & Economic Statistics}, 13(4):467-474.
#'      Nadifar, M., Baghishani, H., and Fallah, A. (2023). A flexible
#'          generalized poisson likelihood for spatial counts constructed by renewal theory, motivated by groundwater
#'            quality assessment.
#'            \emph{Journal of Agricultural, Biological, and Environmental Statistics},
#'                   28:726-748.
#'        \emph{Neutrosophic Sets and Systems},  22, 30-38.
#' @importFrom stats pgamma rgamma runif
#' @examples
#' # In a study, the number of disease incidence, we will calculate
#' # the probability that the number of disease is zero with rate 1
#' dGC(0, alpha = 1, gamma = 1)
#'
#' @export
G <- function(alpha, gamma) {
  if (any(alpha < 0)) {
    stop(message = "alpha must be a positive value")
  }
  if (any(gamma < 0)) {
    stop(message = "gamma must be a positive value")
  }
  return(stats::pgamma(gamma, shape = alpha, rate = 1))
}

#' @name GC
#' @export
dGC <- function(y, alpha, gamma) {
  if (any(alpha < 0)) {
    stop(message = "alpha must be a positive value")
  }
  if (any(gamma < 0)) {
    stop(message = "gamma must be a positive value")
  }
  if (any(y < 0) && any(y - floor(y) == 0)) {
    stop(message = "Warning: x should be a positive integer.")
  }
  dgc <- G(y * alpha, gamma) - G((y + 1) * alpha, gamma)
  return(dgc)
}

#' @name GC
#' @examples
#' # the probability that the disease will receive at least one
#' pGC(q = 1, alpha = 1, gamma = 1, lower.tail = FALSE)
#' # the probability that the disease will receive at most three
#' pGC(q = 3, alpha = 1, gamma = 1, lower.tail = TRUE)
#' @export
pGC <- function(q, alpha, gamma, lower.tail = TRUE) {
  if (any(q < 0) && any(q - floor(q) == 0)) {
    stop(message = "Warning: q should be a  positive integer.")
  }
  if (any(alpha < 0)) {
    stop(message = "alpha must be a positive value")
  }
  if (any(gamma < 0)) {
    stop(message = "gamma must be a positive value")
  }
  T <- gamma * alpha
  pgc <- 1 - stats::pgamma(T, (q + 1) * alpha, 1)
  if (!lower.tail) {
    pgc <- 1 - pgc
  }
  return(pgc)
}

#' @name GC
#' @examples
#' # Calcaute the quantiles
#' qGC(p = c(0.25, 0.5, 0.75), alpha = 1, gamma = 1)
#' @export
qGC <- function(p, alpha = 1, gamma) {
  n <- length(p)
  if (any(p < 0) || any(p > 1)) {
    stop(message = "Warning: p must be between 0 and 1")
  }
  if (alpha < 0) {
    stop(message = "alpha must be a positive value")
  }
  if (any(gamma < 0)) {
    stop(message = "gamma must be a positive value")
  }
  if (length(p) != n) {
    ifelse(length(p) == 1, p <- rep(p, n), stop(message = "length of p is incorrect"))
  }

  if (length(gamma) != n) {
    ifelse(length(gamma) == 1, gamma <- rep(gamma, n),
      stop(message = "length of gamma is incorrect")
    )
  }
  H <- function(y, alpha, gamma, p) {
    pGC(y, alpha, gamma) - p
  }

  qgc <- numeric(n)
  for (i in 1:n) {
    I <- c(0, 40)
    if (H(I[1], alpha, gamma[i], p[i]) * H(I[2], alpha, gamma[i], p[i]) > 0 &&
      H(I[1], alpha, gamma[i], p[i]) > 0) {
      qgc[i] <- 0
    } else {
      qgc[i] <- stats::uniroot(H, I, alpha, gamma[i], p[i], extendInt = "yes")$root
    }
  }
  return(round(qgc))
}

#' @name GC
#' @examples
#' # Simulate 10 values from GC(1, 1)
#' rGC(n = 10, alpha = 1, gamma = 1, method = "PF")
#' @export
rGC <- function(n = 1, alpha = 1, gamma = gamma, method = "PF") {
  if (alpha < 0) {
    stop(message = "alpha must be a positive value")
  }
  if (any(gamma < 0)) {
    stop(message = "gamma must be a positive value")
  }
  if (method == "PF") {
    y <- numeric(n)
    M <- 2000
    prob <- numeric(M + 1)
    ifelse(length(gamma) != n, gamma <- rep(gamma, n),
      gamma <- gamma
    )
    for (i in 1:n) {
      # compute the discrete probability distribution and
      # then sample from it
      for (j in 1:M) {
        yy <- j - 1
        prob[j] <- dGC(yy, alpha, gamma[i])
        ## G(yy * alpha, gamma[i]) - G((yy + 1) * alpha, gamma[i])
      }
      y[i] <- sample(0:M, size = 1, prob = prob)
    }
  }
  if (method == "IT") {
    P <- stats::runif(n)
    qGC(P, alpha, gamma)
  }
  if (method == "WT") {
    ## To save the sum of Gamma ramdom variables.
    s <- 0
    ## To keep each Gamma random number.
    x <- numeric(0)
    repeat {
      ## Generate a Gamma.
      u <- stats::rgamma(n = 1L, shape = alpha, rate = gamma)
      ## Update the sum.
      s <- s + u
      if (s >= n) break
      ## Append the new value.
      x <- append(x, values = u)
    }
    ## Tabulate (count) occurrences in each interval.
    y <- tabulate(ceiling(cumsum(x)))
  }
  return(y)
}
