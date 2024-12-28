#' Compute the Expectation for a Discrete CDF
#' 
#' This function calculates the expectation of a distribution represented by a discrete CDF (Cumulative Distribution Function) of class "stepfun".
#' The expectation is computed along with the positive and negative contributions to the expectation.
#'
#' @param sf A stepfun object representing the discrete CDF.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{E}{The total expectation of the distribution.}
#'   \item{Eplus}{The positive part of the expectation.}
#'   \item{Eminus}{The negative part of the expectation.}
#' }
#'
#' @examples
#' # Example usage with an ECDF
#' sf <- ecdf(mtcars$disp)
#' Ex(sf)
#' Extracting only the expected value
#' Ex(sf)$E
#'
#' @export
Ex <- function( sf ) {
  x_vals <- knots(sf)
  y_vals <- sf(x_vals)
  n <- length(y_vals)
  dx <- diff(x_vals)
  
  #<->#  Máscara lógica
  mask <- x_vals[-n] < 0
  
  segment_areasMinus <- dx[mask] * y_vals[-n][mask]
  segment_areasPlus <- dx[!mask] * (1-y_vals[-n][!mask])
  
  Eplus <- sum(segment_areasPlus)
  Eminus <- sum(segment_areasMinus)
  E <- Eplus - Eminus
  
  #<->#  Retorno la esperanza 
  list(E=E, Eplus=Eplus , Eminus=Eminus )
}
 
