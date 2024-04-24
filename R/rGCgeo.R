#' Generate Data from GC Spatial Regression Model with Geospatial Dependency
#'
#' This function generates spatially dependent count data based on the Gamma-Count (GC) spatial regression model.
#' It uses a specified geospatial dependency model with parameters such as `sigma` for variance and `range` for spatial range.
#' The function returns a list containing the generated data and relevant information about the simulation.
#'
#' @param n Integer. The number of knots (or spatial units) for which the data should be generated.
#' @param alpha Numeric. The dispersion parameter of the Gamma-Count model.
#' @param beta0 Numeric. The intercept term for the model.
#' @param beta Numeric vector. The regression coefficients (fixed effects) for the model.
#' @param V Optional numeric. The variance of the spatial random effects for lattice data.
#' @param rho Optional numeric. The spatial correlation coefficient for the CAR model. Default is 1.
#' @param sigma Optional numeric. The variance of the spatial random effects for geospatial data with Matern covariance.
#' @param range Optional numeric. The range parameter for geospatial data with Matern covariance.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{covariate}{A matrix of covariates with the specified number of knots (`n`) and columns based on the length of `beta`.}
#'   \item{phi}{A vector of spatial random effects generated from the Matern covariance model.}
#'   \item{eta}{A vector representing the linear predictor, calculated as the dot product of the covariates and coefficients plus the spatial effects (`phi`).}
#'   \item{mu}{A vector of mean response values calculated as the product of `alpha` and the exponential of `eta`.}
#'   \item{y}{A vector of simulated count data based on the GC model and the mean response values (`mu`).}
#'   \item{ID}{A vector of knot IDs from 1 to `n`.}
#' }
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
#' \donttest{
#' # Generate data from the GC spatial regression model with geospatial dependency
#'
#' data <- rGCgeo(n = 100, alpha = 1, beta0 = 0.3, beta = c(-0.5, 0.5),
#' sigma = 1, range = 2)
#'
#' # View the generated data
#' print(data)
#' }
#'
#' @export
rGCgeo <- function(n = n, alpha, beta0, beta,
                   V = NULL, rho = 1, sigma = NULL,
                   range = NULL) {
  # Initialize data list
  data <- list()

  # Number of covariates
  P <- length(beta)

  # Generate covariates using multivariate normal distribution
  covariate <- mvtnorm::rmvnorm(n, rep(0, P), diag(P))
  colnames(covariate) <- paste0("x", 1:P)
  data$covariate <- covariate

  # Combine the intercept and coefficients
  Coef <- c(beta0, beta)

  # Create the design matrix
  X <- as.matrix(cbind(1, covariate))

  # Generate spatial random effects using geospatial model with Matern covariance
  phi <- spatGEO(n, sigma, range)
  data$phi <- phi

  # Calculate the linear predictor
  data$eta <- X %*% Coef + phi

  # Calculate the mean response
  data$mu <- alpha * exp(data$eta)

  # Generate count data based on the GC model and the mean response values
  data$y <- sapply(data$mu, FUN = rGC, n = 1, alpha = alpha, method = "WT")

  # Add IDs to the data list
  data$ID <- 1:n

  # Return the data list
  return(data)
}
