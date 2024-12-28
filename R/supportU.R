#' Compute Support for \eqn{U = Y_1 - Y_0}
#'
#' This function computes the support for the outcome difference \eqn{U = Y_1 - Y_0}, where \eqn{Y_1} and \eqn{Y_0} are outcomes
#' corresponding to treatment and control groups, respectively.
#'
#' @param Y A numeric vector of observed outcomes.
#' @param Z A numeric vector indicating the treatment variable (e.g., 1 for treatment and 0 for control).
#' @param y1.limit A vector of length 2 specifying the minimum and maximum values of Y1. If NULL, 
#'                  it defaults to the range of Y in the treatment group.
#' @param y0.limit A vector of length 2 specifying the minimum and maximum values of Y0. If NULL, 
#'                  it defaults to the range of Y in the control group.
#'
#' @return A list containing:
#' \item{ sopY1 }{ Sorted values of Y for the treatment group. }
#' \item{ sopY0 }{ Sorted values of Y for the control group. }
#' \item{ sopU }{ Sorted values of U = Y1 - Y0. }
#' \item{ sopXYU }{ A data frame with all combinations of Y1 and Y0, and their differences (U). }
#'
#' @export
supportU <- function(Y,Z,y1.limit=NULL, y0.limit=NULL){
  
  if ( is.null(y1.limit) ) {
    y1.limit <- c( min(Y,na.rm = T) , max(Y, na.rm = T) )
  }
  if ( is.null(y0.limit) ) {
    y0.limit <- c( min(Y,na.rm = T) , max(Y, na.rm = T) )
  }
  
  y1 <- sort( union(y1.limit , Y[Z==1]) )
  y0 <- sort( union(y0.limit , Y[Z==0]) )
  
  sopXYU <- expand.grid( y1 , y0 )
  sopU <- sopXYU[,1] - sopXYU[,2]
  sopXYU$sopU <- sopU
  
  list( sopY1  = y1 , 
        sopY0  = y0 , 
        sopU   = sort(sopU) , 
        sopXYU = sopXYU ) 
}
 
