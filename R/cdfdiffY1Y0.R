
#' Calculate the Mass and Distribution of the Difference U = Y1 - Y0
#'
#' This function computes the mass and cumulative distribution function (CDF)
#' of the difference U = Y1 - Y0 using a given copula. 
#'
#' @param Copula A copula function.
#' @param FY1 Function representing the CDF of Y1.
#' @param FY0 Function representing the CDF of Y0.
#' @param suppY1 Support of Y1.
#' @param suppY0 Support of Y0.
#'
#' @return A list containing two elements:
#'   \item{jointdist}{A data frame with the joint distribution of (Y1, Y0) under the specified copula.}
#'   \item{massU}{A data frame with the distribution of U = Y1 - Y0.}
#'   \item{cdfU}{CDF of U = Y1 - Y0.}
#'
#' @examples
#' # Define copulas: 
#' W <- function(a,b) { pmax(a+b-1,0)}
#' M <- function(a,b) { pmin(a,b)}
#' PI <- function(a,b) { a*b }
#' 
#' Y1 <- iris$Sepal.Length[iris$Species == "setosa"]
#' Y0 <- iris$Sepal.Length[iris$Species == "virginica"]
#' cdfUPI <- cdfdiffY1Y0(Copula = PI,
#'                       FY1 = ecdf(Y1),
#'                       FY0 = ecdf(Y0),
#'                       suppY1 = sort(unique(Y1)),
#'                       suppY0 = sort(unique(Y0)) )
#' 
#' plot( cdfUPI$cdfU )
#'
#'
#' @export
cdfdiffY1Y0 <- function(Copula, FY1, FY0, suppY1, suppY0){
  
  suppXYU <- expand.grid(suppY1=suppY1, suppY0=suppY0)
  suppXYU$U <- suppXYU[,1] - suppXYU[,2]
  
  # Evaluate the joint CDF values
  FCY1Y0values <- Copula( FY1(suppXYU[,1]) , FY0(suppXYU[,2]) )
  if ( abs( max(FCY1Y0values) - 1) > 2*.Machine$double.eps ) {
    warning( paste0("The maximum value of the joint cdf is different from 1 in ", abs( max(FCY1Y0values) - 1), "." ) ) 
    }
  FCY1Y0values <- cbind(suppXYU, FC=FCY1Y0values)
  
  # Construct matrix with CDF values
  numcol <- length(suppY1)
  FCmatrix <- matrix(FCY1Y0values$FC, ncol = numcol, byrow = TRUE)
  
  # Add zeros to the top and left of the matrix
  FCmatrix <- cbind(0, rbind(0, FCmatrix))
  
  # Calculate masses
  masses <- massCpp(FCmatrix)
  
  # Check if masses sum to 1
  if ( abs( sum(masses) - 1 ) > 2*.Machine$double.eps ) {
    warning( paste0("The sum of masses differs from 1 in ", abs( sum(masses) - 1), "." ) ) 
  }
  
  # Add masses to the data frame
  FCY1Y0values$massFC <- as.vector(t(masses))
  
  # Rename columns
  names(FCY1Y0values) <- c("GrillaY1", "GrillaY0", "U", "FCcdf", "FCmasa")
  
  # Aggregate masses by U
  massU <- aggregate_cpp(FCY1Y0values$FCmasa, FCY1Y0values$U, sort(unique(FCY1Y0values$U)) )
  massU <- as.data.frame(massU)
  names(massU) <- c("U", "massU")
  massU$FUcdf <- cumsum(massU$massU)
  
  # ECDF of U
  cdfU <- approxfun(massU$U, massU$FUcdf, method = "constant", 
                     yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(cdfU) <- c("ecdf", "stepfun", class(cdfU))
  assign("vals", massU$U, envir = environment(cdfU))
  assign("Fvals", massU$FUcdf, envir = environment(cdfU))
  
  # Return results
  list(jointdist = FCY1Y0values, distU = massU, cdfU=cdfU)
}

