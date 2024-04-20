#' Generate Spatial Random Fields from CAR Models
#'
#' This function generates spatial random fields from Conditional Autoregressive (CAR) models
#' for lattice spatial data. Given a neighborhood matrix, a variance parameter, and a spatial
#' dependence parameter, the function produces a spatial random field.
#'
#' @param W Numeric matrix. The neighborhood matrix representing the adjacency relationships
#' between spatial units.
#'          It can be provided by using the function `rAdj`.
#' @param sig Numeric. The variance of the spatial random effects from the CAR model.
#'  Must be positive.
#' @param rho Numeric. The spatial dependence parameter for the CAR model.
#' Must be between -1 and 1, inclusive.
#'
#' @return A numeric vector representing the spatial random field from the CAR model.
#' The length of the vector
#'         is equal to the number of spatial units (rows in the neighborhood matrix).
#'
#' @details
#' The function starts by computing the diagonal matrix of the number of neighbors for each
#'  spatial unit.
#' Then, it calculates the precision matrix (Q) based on the given parameters and
#' neighborhood matrix.
#' A small constant (0.0001) is added to the diagonal to ensure the precision matrix is
#' non-singular.
#' Finally, the covariance matrix is calculated as the inverse of the precision matrix
#' multiplied by the variance parameter (`sig`).
#' The function uses multivariate normal random generation
#' (using `rmvnorm` from the `mvtnorm` package) to produce the spatial random field.
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
#' # Generate a random adjacency matrix for 5 spatial regions with a probability of 0.2
#' W <- rAdj(n = 5, p = 0.2)
#'
#' # Generate a spatial random field from the CAR model using the adjacency matrix
#' # with parameters variance = 0.1, and rho = 0.5
#' spatial_random_field <- spatCAR(W = W, sig = 0.1, rho = 0.5)
#' print(spatial_random_field)
#'
#' @export
spatCAR <- function(W, sig, rho) {
  # Parameter validation
  if (sig <= 0) {
    stop("Variance (sig) must be positive value")
  }
  if (rho < -1 || rho > 1) {
    stop("Rho must be between -1 and 1")
  }

  # Number of spatial units
  ncell <- nrow(W)

  # Calculate the number of neighbors for each spatial unit
  n_neighbors <- rowSums(W)

  # Calculate the precision matrix (Q) and add a small constant to diagonal to avoid singularity
  Q <- diag(n_neighbors) - rho * W + diag(1e-4, ncell)

  # Calculate the covariance matrix based on the precision matrix and variance
  cov_matrix <- sig * solve(Q)

  # Generate spatial random field from CAR model
  spatial_random_field <- as.vector(mvtnorm::rmvnorm(1,
    mean = rep(0, ncell),
    sigma = cov_matrix
  ))

  # Return the spatial random field
  return(spatial_random_field)
}
