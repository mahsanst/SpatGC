#' Fit ICAR Spatial Gamma-Count Model
#'
#' This function fits an ICAR spatial Gamma-Count (GC) model to a given dataset using the INLA package.
#' It constructs the formula based on the provided covariate data and ID variables, and fits the model using
#' the specified adjacency matrix (`W`) or a shapefile of the study region.
#'
#' @param Y Vector of response variables (counts).
#' @param ID Vector of indexes of regions (spatial units).
#' @param W Optional adjacency matrix representing spatial connections between regions.
#'          If not provided, it can be generated from a shapefile using the `shapefile` argument.
#' @param shapefile Optional shapefile representing the study region. If provided, the adjacency matrix (`W`)
#'                  will be calculated from the shapefile.
#' @param covariate Optional matrix of covariates. If not provided, the function assumes the model
#'                 is intercept-only.
#'
#' @return An object of class "inla" representing the fitted ICAR spatial GC model.
#'         The object contains model estimates, diagnostics, and other results.
#'
#' @importFrom spdep poly2nb
#' @importFrom stats as.formula
#' @importFrom sf st_make_valid
#' @examples
#' # Generate data from the GC spatial regression model with lattice spatial dependency
#' W <- rAdj(500) # Generate a random adjacency matrix
#' DDl <- rGClat(n = 500, alpha = 1, beta0 = 0.3, beta = c(-0.5, 0.5), W = W, V = 1)
#'
#' # Prepare the data
#' Y <- DDl$y
#' covariate <- DDl$covariate
#' ID <- DDl$ID
#'
#' # Fit the spatial GC model
#' ResultGC <- GClat(Y = Y, ID = ID, covariate = covariate, W = W)
#'
#' # Summary of the model fit
#' summary(ResultGC)
#'
#' @export
GClat <- function(Y, ID, W = NULL, shapefile = NULL, covariate = NULL) {
  # Construct adjacency matrix from shapefile if provided
  if (!is.null(shapefile)) {
    W <- spdep::nb2mat(
      neighbours = spdep::poly2nb(sf::st_make_valid(shapefile)),
      style = "B", zero.policy = TRUE
    )
  }

  # Prepare the formula
  vars <- colnames(covariate)
  inlaform <- paste(
    "Y ~",
    paste(
      paste(vars, collapse = "+"),
      paste("f(ID, model = 'besag', graph = W, constr = TRUE)"),
      sep = "+"
    )
  )
  formula <- stats::as.formula(inlaform)

  # Prepare the data frame
  data <- data.frame(Y, covariate, ID)

  # Fit the Gamma-Count model
  gc <- INLA::inla(
    formula = formula,
    data = data,
    control.family = list(hyper = list(
      "alpha" = list(prior = "pc.gammacount", param = c(3))
    )),
    family = "gammacount",
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE)
  )

  # Return the fitted model
  return(gc)
}
