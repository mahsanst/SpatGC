#' Generate Spatial Random Fields from Matern Covariance Model
#'
#' This function generates spatial random fields from the Matern covariance model for geospatial data.
#' Given a number of knots, a variance parameter, and a range parameter, the function produces
#' a spatial random field from the Matern covariance model.
#'
#' @param m Integer. The number of knots (or spatial units) for which the spatial random field should be generated.
#' @param sigma Numeric. The variance of the spatial random effects with Matern covariance. Must be positive.
#' @param range Numeric. The range of the spatial random effects with Matern covariance. Must be positive.
#'
#' @return A numeric vector representing the spatial random field from the Matern covariance model.
#'         The length of the vector is equal to the number of knots specified by the parameter `m`.
#'
#' @details
#' The function starts by converting the variance (`sigma`) and range (`range`) parameters into the tau and kappa
#' parameters required for the Matern covariance model. Then, it constructs a mesh for the given number of knots,
#' defining the spatial layout of the data. The function uses INLA functions (`inla.spde2.matern` and others) to
#' create an SPDE object and set the priors for the spatial model. It then projects the mesh and constructs a predictor
#' for the model. Finally, the function calculates the precision matrix, samples from the distribution, and computes
#' the spatial random field using the A matrix and sampled values.
#'
#' @importFrom stats runif
#'
#' @examples
#'
#' # Generate a spatial random field from the Matern covariance model
#' # using 100 knots, variance = 1, and range = 1
#' # (First uncomment, then run it, please)
#' # spatial_random_field <- spatGEO(m = 100, sigma = 1, range = 1)
#' # print(spatial_random_field)
#'
#' @export
spatGEO <- function(m, sigma, range) {
  # Parameter validation
  if (sigma <= 0) {
    stop("Variance (sigma) must be positive value")
  }
  if (range <= 0) {
    stop("Range must be positive value")
  }

  # Convert variance and range into tau and kappa parameters
  kappa0 <- sqrt(8) / range
  tau0 <- 1 / (sqrt(4 * pi) * kappa0 * sigma)

  # Construct mesh with the specified number of knots
  points <- matrix(stats::runif(m * 2), m, 2)
  mesh <- INLA::inla.mesh.create.helper(
    points = points,
    cutoff = 0.05,
    offset = c(0.1, 0.4),
    max.edge = c(0.05, 0.5)
  )

  # Construct SPDE object and set priors for use with INLA
  spde <- INLA::inla.spde2.matern(
    mesh,
    B.tau = cbind(log(tau0), 1, 0),
    B.kappa = cbind(log(kappa0), 0, 1),
    theta.prior.mean = c(0, 0),
    theta.prior.prec = c(0.1, 1)
  )

  # Project the mesh and create a predictor
  proj <- INLA::inla.mesh.projector(
    mesh,
    xlim = c(0, 1),
    ylim = c(0, 1),
    dims = c(m, m)
  )

  # Create the A matrix for the mesh
  A <- INLA::inla.spde.make.A(mesh, loc = points)

  # Compute the precision matrix
  Q <- INLA::inla.spde.precision(spde, theta = c(0, 0))

  # Sample from the distribution using the precision matrix
  x <- as.vector(INLA::inla.qsample(n = 1, Q))

  # Calculate the spatial random field using the A matrix and sampled values
  phi <- as.vector(A %*% x)

  # Return the spatial random field
  return(phi)
}
