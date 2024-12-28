
#' Plot Method for ecdfPI Objects
#'
#' This function creates a plot for objects of class \code{ecdfPI}.
#'
#' @param x An object of class \code{ecdfPI}.
#' @param F1 Lower bound of the ECDF. If `NULL`, bounds under no assumptions are used.
#' @param F2 Upper bound of the ECDF. If `NULL`, bounds under no assumptions are used.
#' @param ... Additional arguments passed to \code{plot.stepfun}.
#' @param ylab A title for the y-axis (default is \code{"Fn(x)"}).
#' @param verticals Logical; if \code{TRUE}, vertical lines are added at data points.
#' @param col.01line Color for the horizontal lines at 0 and 1 (default is \code{"gray70"}).
#' @param pch Plotting character or symbol to use (default is \code{19}).
#' @param x.limit Numeric vector of length 2 specifying the minimum and maximum x-values for plotting bounds.
#'                Defaults to 20 times the range of the data if not specified.
#' @param col.bounds Color for the bands of identification (default is \code{"blue"}).
#' @param alpha.bounds Transparency level for the bands of partial identification (default is \code{0.1}).
#'
#' @details This function is a fork of \code{plot.ecdf} for \code{ecdfPI} objects.
#'
#' @export

plot.ecdfPI <- function (x, F1=NULL, F2=NULL, ..., ylab = "Fn(x)", verticals = FALSE, col.01line = "gray70", 
                         pch = 19, x.limit = NULL, col.bounds="blue", alpha.bounds=0.1) {
  plot.stepfun(x$ecdf, ..., ylab = ylab, verticals = verticals, 
               pch = pch)
  abline(h = c(0, 1), col = col.01line, lty = 2)
  #<-># Si no especifico cotas, uso las sin supuestos (incluídas en x)
  if ( is.null(F1) ) { F1 <- x$F1 }
  if ( is.null(F2) ) { F2 <- x$F2 }
  plotecdfbound(F1=F1, F2=F2, x.limit=x.limit, col.bounds=col.bounds, alpha.bounds=alpha.bounds)  
} 

