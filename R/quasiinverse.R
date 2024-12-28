#' Calculate quasi-inverse of ecdf functions
#'
#' This function calculates the quasi-inverse of empirical cumulative distribution function (ECDF). 
#'
#' @param f A function of class \code{ecdf} or \code{ecdfPI}.
#' @param continuity Continuity side (default: right).
#' @param ymin,ymax Minimum and maximum values of the random variable. If `NULL`, the data’s minimum and maximum values are used (`NULL` by default).
#'
#'
#' @return A function of class \code{quasiinverse}.
#'
#' @examples
#' result <- ecdfPI(iris$Sepal.Length)
#' plot( quasiinverse(result, 'left') , col="blue" , 
#'       main = "Quasi-inverses of iris$Sepal.Length" ,
#'       xlab = "q" , ylab = "y")
#' lines( quasiinverse(result, 'right') , col="brown")
#' 
#' legend("topleft", 
#'        legend = c(expression(F^{-1} * (q) == "inf { y : F(y) >= q}"),
#'                   expression(F^{(-1)} * (q) == "sup { y : F(y) <= q}")),
#'        col = c("brown", "blue"), 
#'        lwd = 2, 
#'        box.lty = 0) 
#'
#' @export
quasiinverse <- function(f, continuity='right', ymin=NULL, ymax=NULL){
  vals <- get('vals', envir = environment(f$ecdf) )
  Fvals <- get('Fvals', envir = environment(f$ecdf) )
  if (continuity=="right") {
    ff <- 0
  } else if (continuity=="left") {
    ff <- 1
  } else {
    ff <- 0
    stop("Error: 'continuity' must be either 'left' or 'right'.")
    } 
  if ( is.null(ymin) ) {
    ymin <- min(vals)
  }
  if ( is.null(ymax) ) {
    ymax <- max(vals)
  }
  f <- approxfun(Fvals, vals, method = "constant", yleft = ymin, yright = ymax, 
                 f = ff, ties = "ordered")
  class(f) <- c("quasiinverse","stepfun", class(f) )
  assign("nobs", length(vals), envir = environment(f))
  assign("vals", vals, envir = environment(f))
  assign("Fvals", Fvals, envir = environment(f))
  f
}

