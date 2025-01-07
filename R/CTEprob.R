#' Calculate the distribution of the Causal Treatment Effect (CTE) for a Dichotomous Outcome
#'
#' This function calculates the CTE(0) distribution for a dichotomous outcome over all possible worlds
#' that are consistent with the observed data. 
#'
#' @param M Total number of individuals in the population.
#' @param nNA Total number of missing data in the observed outcome.
#' @param nZ1K1 Number of statistical units with Z=1 and K=1.
#' @param nZ0K1 Number of statistical units with Z=0 and K=1.
#' @param nY1Z1K1 Number of statistical units with Y=1, Z=1, and K=1.
#' @param nY1Z0K1 Number of statistical units with Y=1, Z=0, and K=1.
#' @param p Propensity values for which the CTE distribution is calculated.
#' @param dec.prec Number of decimal places for precision in the calculations (default is 13).
#'
#' @details 
#' The variable \eqn{\mathcal{M}} represents the set of labels for the individuals in the population of interest.
#' Observed causal inference data are typically presented as variables defined for each individual \eqn{m \in M}.
#' The outcome \eqn{Y} is defined over the space \eqn{\mathcal{M} \times \{0, 1\}}, where \eqn{K(m, z) = 1} if \eqn{Y(m, z)} is observed, 
#' and 0 otherwise. \eqn{Z} is the projection of \eqn{(m, z)} in \eqn{Z}, meaning \eqn{Z(m, z) = z}.
#'
#' The outcome \eqn{Y} is dichotomous, and each missing outcome value can take either 0 or 1. 
#' Consider the space of all possible worlds where for each world, there exists an outcome \eqn{Y^*},
#' where \eqn{Y^*} coincides with the observed outcome \eqn{Y_{\text{obs}}} for the observed statistical units.
#' For each \eqn{Y^*}, the CTE(0) is computed. The CTE(0) is a random variable over all possible worlds.
#' 
#' The value \eqn{p} represents the prior probability that each unobserved outcome is 1. This prior probability is known as
#' the propensity for success. Assuming a common \eqn{p} for all unobserved outcomes implies the absence of a causal effect.
#' Based on \eqn{p}, the probability of each possible world is determined. If all possible worlds are equally likely, \eqn{p = 0.5},
#' which assumes no causal effect. This assumption becomes unrealistic when the outcome's success probability is significantly 
#' different from 0.5.
#'
#' @return 
#' A list containing the CTE(0) distribution for a dichotomous outcome in the measurable space. The list includes:
#' 
#' \item{masses}{A data frame with the probability mass \eqn{P(CTE(0) = cte)}.}
#' \item{cdf}{A data frame with the cumulative probability \eqn{P(CTE(0) \leq cte)}.}
#' \item{probMZ}{Probabilities identified from the data:}
#'  \itemize{
#'    \item{probK1Z1}{\eqn{P(K = 1 | Z = 1)}}
#'    \item{probK1Z0}{\eqn{P(K = 1 | Z = 0)}}
#'    \item{probY1Z1K1}{\eqn{P(Y = 1 | Z = 1, K = 1)}}
#'    \item{probY1Z0K1}{\eqn{P(Y = 1 | Z = 0, K = 1)}}
#'    \item{p0}{The prior propensity where the probability of a positive causal effect equals the probability of a negative causal effect.}
#'  }
#'
#' @examples
#' # Example usage of the CTEprob function
#' result <- CTEprob(M = 200, nNA = 0, nZ1K1 = 65, nZ0K1 = 135, 
#'                   nY1Z1K1 = 34, nY1Z0K1 = 91, p = c(0:100)/100 )
#' result$masses
#' result$cdf
#'
#' @export
CTEprob <- function(M, nNA = 0, nZ1K1, nZ0K1, nY1Z1K1, nY1Z0K1, p, dec.prec=13) {
  
  n1 <- M - nZ1K1      # nZ0K1 + nNA
  n2 <- nZ1K1 + nNA
  
  a1 <- nY1Z1K1 / M    # P(Y=1|Z=1,K=1) P(K=1|Z=1)
  b1 <- 1.0 / M
  a2 <- nY1Z0K1 / M    # P(Y=1|Z=0,K=1) P(K=1|Z=0)
  b2 <- 1.0 / M
  
  probK1Z1   <- (M-n1) / M       # P(K=1|Z=1)
  probK1Z0   <- (M-n2) / M       # P(K=1|Z=0) 
  probY1Z1K1 <- nY1Z1K1 / nZ1K1  # P(Y=1|Z=1,K=1)
  probY1Z0K1 <- nY1Z0K1 / nZ0K1  # P(Y=1|Z=0,K=1)
  
  probMZ <- list(probK1Z1 = probK1Z1, probK1Z0 = probK1Z0, 
                 probY1Z1K1 = probY1Z1K1, probY1Z0K1 = probY1Z0K1)
  
  if (probK1Z1 == 0.5) {
    probMZ$p0 <- "p0 does not exist"
  } else {
    p0 <- (probY1Z1K1 * probK1Z1 - probY1Z0K1 * probK1Z0) / (probK1Z1 - probK1Z0)
    probMZ$p0 <- p0
  }
  
  obj <- CTEprobcpp(n1 = n1, n2 = n2, a1 = a1, b1 = b1, a2 = a2, b2 = b2, p = p)
  colnames(obj$results) <- obj$namesresults
  obj$results[,"cte"] <- round( obj$results[,"cte"] , dec.prec)
  
  results <- obj$results
  
  masses <- data.frame( cte = sort(unique(results[,"cte"])) )
  cdf <- masses
  
  for (i in 1:length(p)) {
    aux <- aggregate_cpp( values = results[,5+i],
                          groups = results[,"cte"] ,
                          levelgroup = masses[,"cte"] )
    masses <- cbind( masses, aux[,2] )
    cdf <- cbind( cdf, cumsum(aux[,2]) )
  }
  
  names(masses)[-1] <- paste0("p", p)
  names(cdf) <- names(masses)
  
  list( masses=masses, cdf=cdf , p=p , probMZ=probMZ)
}

