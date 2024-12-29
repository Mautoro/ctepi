
#' Plot Partial Identification Bounds
#'
#' This function plots the partial identification bounds for a given range of x values
#' based on an object of class "ecdfPI". It creates a shaded area between the lower and upper bounds.
#'
#' @param F1, F2 Lower and upper partial identification bounds. F1 and F2 must be step functios.
#' @param x.limit A vector of length 2 specifying the minimum and maximum x-values for the plot.
#'                If NULL, defaults to 20 times the range of the x-values in the object.
#' @param col.bounds Color of the shaded area representing the partial identification bounds.
#' @param alpha.bounds Transparency level of the shaded area.
#'
#' @importFrom grDevices rgb
#' @export
plotecdfbound <- function( F1, F2, x.limit = NULL, col.bounds="blue", alpha.bounds=0.1 ){
  vals1 <- get("x", envir = environment(F1) )
  vals2 <- get("x", envir = environment(F2) )
  
  if (length(x.limit) > 0) {
    if (length(x.limit) != 2) {
      stop("Error: x.limit must be a vector of length 2.")
    } else {
      #<-># Ordeno x.limit por si viene al revés (ej: c(6,4))
      x.limit <- range(x.limit)
      #<-># Restrinjo valores a los límites indicados
      vals1 <- unique(c( x.limit[1] , vals1[vals1 >= x.limit[1] & vals1 <= x.limit[2]] , x.limit[2] ))
      vals2 <- unique(c( x.limit[1] , vals2[vals2 >= x.limit[1] & vals2 <= x.limit[2]] , x.limit[2] ))
    }
  }
  if (is.null(x.limit)) {
    #<-># Para definir x.limit da lo mismo si tomo vals1 o vals2
    x.limit <- diff(range(vals1)) * c(-20,20) + range(vals1)
  }
  #<-># Grafico cotas
  xx1aux <- rep( c(x.limit[1],vals1,x.limit[2]) ,each=2)
  xx2aux <- rep( c(x.limit[1],vals2,x.limit[2]) ,each=2)
  nn1 <- length(xx1aux)
  nn2 <- length(xx2aux)
  xx1 <- xx1aux[-c(1,nn1)]
  yy1 <- F1(xx1aux[-c(nn1-1,nn1)])
  xx2 <- rev(xx2aux[-c(1,nn2)])
  yy2 <- F2( rev(xx2aux)[-c(1,2)])
  
  polygon( c(xx1, xx2) , c(yy1,yy2) ,  
           border = NA , col=AlphaCol(col.bounds,alpha.bounds))  
} 

