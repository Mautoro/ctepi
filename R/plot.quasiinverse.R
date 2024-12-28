
#' Plot Method for quasiinverse Objects
#'
#' This function creates a plot for objects of class \code{quasiinverse}. It uses the base R function \code{plot.stepfun} 
#' to plot the empirical cumulative distribution function (ECDF) and adds additional elements to the plot,
#' such as horizontal lines at 0 and 1, and bands of identification using \code{plotecdfbound()}.
#'
#' @param x An object of class \code{quasiinverse}.
#' @param ... Additional arguments passed to \code{plot.stepfun}.
#' @param ylab A title for the y-axis (default is \code{"Fn(x)"}).
#' @param verticals Logical; if \code{TRUE}, vertical lines are added at data points.
#' @param col.mmline Color for the horizontal lines at minimum and maximum (default is \code{"gray70"}).
#' @param pch Plotting character or symbol to use (default is \code{19}).
#'
#' @details This function is a fork of \code{plot.ecdf} for \code{ecdfPI} objects.
#'
#' @export

plot.quasiinverse <- function (x, F1=NULL, F2=NULL, ..., ylab = "Fn^{-1}(x)", verticals = FALSE, col.mmline = "gray70", 
                         pch = 19, x.limit = NULL, col.bounds="blue", alpha.bounds=0.1) {
  plot.stepfun(x, ..., ylab = ylab, verticals = verticals, 
               pch = pch)
  #abline(h = c(0, 1), col = col.mmline, lty = 2)
} 

