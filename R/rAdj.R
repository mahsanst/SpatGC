#' Generate a Random Adjacency Matrix
#'
#' This function generates a random adjacency matrix for a given number of regions.
#' The matrix is symmetric, and the upper triangular part (excluding the diagonal) is
#' randomly generated using the \code{\link[stats]{rbinom}} function with a specified
#' probability. The lower triangular part is filled in by reflecting the upper triangular
#' part to ensure the matrix is symmetric.
#'
#' @param n Integer. The number of regions, which determines the dimensions of the
#' adjacency matrix. Must be a positive integer.
#'
#' @param prob Numeric. The probability of an edge (connection) between two regions,
#' used as the probability for the binomial distribution. Default value is 0.2.
#' Should be between 0 and 1 (inclusive).
#'
#' @return A symmetric numeric matrix of dimensions \code{n x n} representing the
#' adjacency matrix. The values in the matrix are either 0 or 1, where 1 indicates
#' an edge (connection) between two regions.
#'
#' @examples
#' # Generate a random adjacency matrix for 5 regions with a probability of 0.3
#' random_adj_matrix <- rAdj(5, prob = 0.3)
#' print(random_adj_matrix)
#'
#' # Check if the matrix is symmetric
#' all(random_adj_matrix == t(random_adj_matrix))
#'
#' @importFrom stats rbinom
#' @export
rAdj <- function(n, prob = 0.2) {
  C <- matrix(0, n, n)
  C[upper.tri(C, diag = FALSE)] <- stats::rbinom(n * (n - 1) / 2, 1, prob)
  W <- C + t(C)
  return(W)
}
