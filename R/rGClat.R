#' Generate Data from GC Spatial Regression Model with Lattice Spatial Effect
#'
#' This function generates spatially dependent count data based on the Gamma-Count (GC) spatial regression model.
#' It uses a specified spatial dependency model (either ICAR or CAR) and optional adjacency matrix or shapefile for spatial relationships.
#' The function returns a list containing the generated data and relevant information about the simulation.
#'
#' @param n Integer. The number of knots (or spatial units) for which the data should be generated. If a shapefile or adjacency matrix (`W`) is provided, this will be determined from those inputs.
#' @param alpha Numeric. The dispersion parameter of the Gamma-Count model.
#' @param beta0 Numeric. The intercept term for the model.
#' @param beta Numeric vector. The regression coefficients (fixed effects) for the model.
#' @param spatial Character. Specifies the type of spatial dependency to use. Options are "ICAR" for Intrinsic Conditional Autoregressive, or "CAR" for Conditional Autoregressive.
#' @param W Optional matrix. The adjacency matrix for lattice data. If provided, it will be used to define spatial relationships between knots.
#' @param V Optional numeric. The variance of the spatial random effects for lattice data.
#' @param rho Optional numeric. The spatial correlation coefficient for the CAR model. Default is 1.
#' @param shapefile Optional. A shapefile defining the spatial relationships between knots. If provided, it will be used to define an adjacency matrix.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{covariate}{A matrix of covariates with the specified number of knots (`n`) and columns based on the length of `beta`.}
#'   \item{phi}{A vector of spatial random effects, generated based on the specified spatial dependency model (`spatial`).}
#'   \item{eta}{A vector representing the linear predictor, calculated as the dot product of the covariates and coefficients plus the spatial effects (`phi`).}
#'   \item{y}{A vector of simulated count data based on the GC model and the linear predictor (`eta`).}
#'   \item{ID}{A vector of knot IDs from 1 to `n`.}
#' }
#' @importFrom mvtnorm rmvnorm
#' @importFrom sf st_make_valid
#'
#' @examples
#' # Generate a random adjacency matrix for a 429x429 grid
#' # (First uncomment, then run it, please)
#' # W <- rAdj(429)
#'
#' # Generate data from the GC spatial regression model with the specified parameters
#' # data <- rGClat(n = 200, alpha = 1, beta0 = 0.3, beta = c(-0.5, 0.5),
#' #  spatial = "ICAR", W = W, V = 1)
#'
#' # View the generated data
#' # print(data)
#'
#' @export
rGClat <- function(n = n, alpha, beta0, beta, spatial = "ICAR", W = NULL,
                   V = NULL, rho = 1, shapefile = NULL) {
  # Parameter validation and adjacency matrix generation
  if (!is.null(shapefile)) {
    W <- spdep::nb2mat(neighbours = spdep::poly2nb(sf::st_make_valid(shapefile)), style = "B", zero.policy = TRUE)
    n <- ncol(W)
  }
  if (!is.null(W)) {
    n <- ncol(W)
  }

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

  # Generate spatial random effects based on the specified spatial model
  if (spatial == "ICAR") {
    phi <- as.numeric(spatICAR(W, V))
  } else {
    phi <- as.numeric(spatCAR(W, V, rho))
  }

  # Store the spatial random effects in the data list
  data$phi <- phi

  # Calculate the linear predictor
  data$eta <- X %*% Coef + phi

  # Generate count data based on the GC model and the linear predictor
  data$y <- rGC(n, alpha, gamma = exp(data$eta))

  # Add IDs to the data list
  data$ID <- 1:n

  # Return the data list
  return(data)
}
