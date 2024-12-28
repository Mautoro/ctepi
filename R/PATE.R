
#' Partial Average Treatment Effect (PATE)
#'
#' This function estimates the Partial Average Treatment Effect (PATE), and the positive and negative parts of the PATE under ignorability assumption.
#'
#' @param Y A numeric vector of observed outcomes.
#' @param Z A numeric vector of treatment variables.
#' @param y.limit A numeric vector of length 2 specifying the limits for the outcome. Defaults to the minimum and maximum of `Y`.
#'
#' @return A list with the following elements:
#' \item{CTEemp}{The CTE function under ignorability assumption.}
#' \item{PATEemp}{The PATE function under ignorability assumption.}
#' \item{y}{A numeric vector of sequence values for the outcome range.}
#' \item{PATEempPlus}{The positive part of the PATE function under ignorability assumption.}
#' \item{PATEempMinus}{The negative part of the PATE function under ignorability assumption.}
#'
#' @examples
#' # Example usage:
#' result <- PATE(Y = c(1, 2, 3, 4, 5), Z = c(0, 1, 0, 1, 1))
#' 
#' @export
PATE <- function(Y,Z, y.limit=NULL){
  
  if ( is.null(y.limit) ) {
    y.limit <- c( min(Y,na.rm = T) , max(Y, na.rm = T) )
  } else if ( length(y.limit) !=2 ) {
    stop("y.limit must be a vector of two numbers.")
  }
  
  y <- Y |> unique() |> sort()
  func <- ctefunctions(Y,Z)
  
  PATEemp <- integrate.sf( func$CTE )
  PATEempPlus <- integrate.sf( func$CTEplus )
  PATEempMinus <- integrate.sf( func$CTEminus )

  list( CTEemp=func$CTE , PATEemp=PATEemp, y=y ,
        PATEempPlus=PATEempPlus , PATEempMinus=PATEempMinus )
}

