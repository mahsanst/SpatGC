#' Fit Negative Binomial Spatial Model
#'
#' This function fits a negative binomial spatial model (including zero-inflated variations)
#' to a given dataset using the INLA package. The function constructs the formula based on
#' the provided covariate data and ID variables, and fits the model using the specified
#' adjacency matrix (`W`) or a shapefile of the study region.
#'
#' @param Y Vector of response variables (counts).
#' @param ID Vector of indexes of regions (spatial units).
#' @param W Optional adjacency matrix representing spatial connections between regions.
#'          If not provided, it can be generated from a shapefile using the `shapefile` argument.
#' @param shapefile Optional shapefile representing the study region. If provided, the adjacency matrix (`W`)
#'                  will be calculated from the shapefile.
#' @param covariate Optional matrix of covariates. If not provided, the function assumes the model
#'                 is intercept-only.
#' @param family The family of negative binomial models to use. Options are "nbinomial",
#'               "zeroinflatedbinomial0", and "zeroinflatedbinomial1".
#'
#' @return An object of class "inla" representing the fitted negative binomial spatial model.
#'         The object contains model estimates, diagnostics, and other results.
#' @importFrom spdep poly2nb
#' @importFrom stats as.formula
#' @importFrom sf st_make_valid
#' @examples
#' # Generate data from the GC spatial regression model with lattice spatial dependency
#' # W <- rAdj(500) # Generate a random adjacency matrix
#' # DDl <- rGClat(n = 200, alpha = 1, beta0 = 0.3, beta = c(-0.5, 0.5),
#' #  W = W, spatial = "lattice", V = 1)
#'
#' # Prepare the data
#' # Y <- DDl$y
#' # covariate <- DDl$covariate
#' # ID <- DDl$ID
#'
#' # Fit the spatial negative binomial model
#' # ResultNB <- NBlat(Y = Y, ID = ID, covariate = covariate, W = W, family = "nbinomial")
#'
#' # Summary of the model fit
#' # summary(ResultNB)
#'
#' @export
NBlat <- function(Y, ID, W = NULL, shapefile = NULL, covariate = NULL,
                  family = c("nbinomial", "zeroinflatedbinomial0", "zeroinflatedbinomial1")) {
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

  # Determine the family of the model
  Family <- match.arg(family)

  # Fit the negative binomial family model
  NB <- INLA::inla(
    formula = formula,
    data = data,
    family = Family,
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE)
  )

  # Return the fitted model
  return(NB)
}
