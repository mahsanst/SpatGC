#' Generate Spatial Random Fields from ICAR Models
#'
#' This function generates spatial random fields from Intrinsic Conditional Autoregressive (ICAR)
#' models for lattice spatial data. Given an adjacency matrix and a variance parameter, the function
#' produces a spatial random field from the ICAR model.
#'
#' @param W Numeric matrix. The adjacency matrix representing the adjacency relationships between spatial units.
#'          It can be provided by using the function `generate_adjacency_matrix`.
#' @param sig Numeric. The variance of the spatial random effects from the ICAR model. Must be positive.
#'            Default value is 1.
#'
#' @return A numeric vector representing the spatial random field from the ICAR model.
#'         The length of the vector is equal to the number of spatial units (columns in the adjacency matrix).
#'
#' @details
#' The function starts by computing the degree of each spatial unit (number of neighbors) from the adjacency matrix.
#' Then, it constructs the precision matrix (Q) based on the adjacency matrix. The function computes the eigenvalues and eigenvectors of Q.
#' To generate the spatial random field, it uses a standard normal distribution to generate a vector of random numbers,
#' scales the random numbers by the variance and the inverse of the eigenvalues (excluding the smallest eigenvalue), and then transforms
#' the random numbers using the eigenvectors. The function returns a spatial random field vector.
#'
#' @importFrom stats rnorm
#'
#' @examples
#' # Define an adjacency matrix for 5 spatial units
#' W <- rAdj(n = 5, p = 0.2)
#'
#' # Generate spatial random field from ICAR model with variance = 1
#' spatial_random_field_icar <- spatICAR(W = W, sig = 1)
#' print(spatial_random_field_icar)
#'
#' @export
spatICAR <- function(W, sig = 1) {
  # Parameter validation
  if (sig <= 0) {
    stop("Variance (sig) must be positive")
  }

  # Number of spatial units (columns in the adjacency matrix)
  n <- ncol(W)

  # Compute the degree of each spatial unit
  delta <- rowSums(W)

  # Construct the precision matrix (Q)
  Q <- -W
  diag(Q) <- delta

  # Compute the eigenvalues and eigenvectors of the precision matrix
  eigen_decomp <- eigen(Q)
  Q_W <- eigen_decomp$vectors[, order(eigen_decomp$values)]

  # Eigenvalues of the precision matrix
  D <- sort(eigen_decomp$values)

  # Generate standard normal random values
  RN <- stats::rnorm(n - 1, mean = 0, sd = sqrt(sig * (1 / D[-1])))

  # Generate spatial random field using the ICAR model
  phi <- Q_W %*% c(0, RN)

  # Return the spatial random field as a numeric vector
  return(as.vector(phi))
}
