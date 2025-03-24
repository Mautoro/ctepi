#' Calculate the Causal Treatment Effect (CTE) distribution for real data with dichotomous outcome
#'
#' This function calculates the Causal Treatment Effect (CTE) distribution for real data with a dichotomous outcome variable.
#'
#' @param Yobs A vector of observed outcomes. A dichotomous variable. If \code{Yobs} is not dichotomous, it is 
#' transformed into a dichotomous variable using the threshold \code{y}.
#' @param Zobs A vector of observed treatments (\code{0} or \code{1}).
#' @param X A matrix of covariates for each individual. It can be used for calculating conditional \eqn{CTE(0|X=x)} distributions.
#' @param p A vector of prior propensity values for which the CTE distribution is calculated. If \code{p} is equal to "PY1K1", then the prior propensity \eqn{p} is set equal to the probability \eqn{P(Y=1|K=1)}.
#' @param y If \code{NULL}, \code{Yobs} is treated as a dichotomous variable. If provided, it is used as a threshold to create 
#' a dichotomous outcome variable: \code{Y > y}.
#' @param covariates A logical value indicating whether the \eqn{CTE(0)} distribution should be calculated for each subset defined by the covariates 
#' (\code{TRUE}) or for the whole dataset (\code{FALSE}).
#'
#' @details 
#' \code{CTEprobdata} accepts observed data \code{(Yobs, Zobs, X)}. Each record corresponds to an observation for an individual. 
#' If \eqn{\mathcal{M}} represents the set of labels for the individuals, then \code{(Yobs, Zobs, X)} are variables in \eqn{\mathcal{M}}.
#' The outcome \eqn{Y} is defined over the space \eqn{\mathcal{M} \times \{0, 1\}}. \eqn{Y} has at least \code{M} missing data points (counterfactuals).
#' Each missing outcome can take values \code{0} or \code{1}.
#' 
#' Consider the space of all possible worlds, where for each world there exists an outcome \eqn{Y^*}, where \eqn{Y^*} coincides 
#' with \code{Yobs} for the observed statistical units. For each \eqn{Y^*}, the \eqn{CTE(0)} is calculated. \eqn{CTE(0)} is a random 
#' variable over all possible worlds.
#'
#' The value \eqn{p} represents the prior probability that each unobserved outcome is \code{1}. This prior probability is 
#' called the "propensity for success." Assuming a common \eqn{p} for all unobserved outcomes implies the absence of a causal effect.
#' Based on \eqn{p}, the probability of each possible world is determined. If all possible worlds are equally likely, then \eqn{p = 0.5}, 
#' implying no causal effect. This assumption is unrealistic when the outcome has a probability of success significantly 
#' different from \code{0.5}.
#' 
#' The \eqn{CTE(0)} distribution over the space of possible worlds is calculated for different values of the prior propensity \eqn{p}. 
#' This allows for the assessment of the impact of the assumption about \eqn{p} on the CTE(0) distribution, and consequently, 
#' the potential existence of a positive or negative causal effect.
#' 
#' The function \code{plotCTEprob} can be used to visualize the probabilities of positive, negative, or no causal effect 
#' for a given prior propensity \eqn{p}.
#'
#' @return 
#' If \code{covariates} is \code{FALSE}, \code{CTEprobdata} returns a list with the CTE(0) distribution over the possible worlds
#' for each value in the \code{p} vector. The list also includes \code{probMZ}, containing probabilities identified in the 
#' \eqn{\mathcal{M} \times \{0, 1\}} space:
#' \itemize{
#'   \item{\code{P(K=1|Z=1)}}: Probability that \eqn{K=1} given \eqn{Z=1}.
#'   \item{\code{P(K=1|Z=0)}}: Probability that \eqn{K=1} given \eqn{Z=0}.
#'   \item{\code{P(Y=1|Z=1,K=1)}}: Probability that \eqn{Y=1} given \eqn{Z=1} and \eqn{K=1}.
#'   \item{\code{P(Y=1|Z=0,K=1)}}: Probability that \eqn{Y=1} given \eqn{Z=0} and \eqn{K=1}.
#'   \item{\code{p0}}: The prior propensity where the probability of a positive causal effect is equal to the probability of a negative causal effect.
#' }
#' 
#' If \code{covariates} is \code{TRUE}, \code{CTEprobdata} returns a list for each subset defined by the covariates \code{X}. 
#' The name of each list corresponds to the values of the covariates. If the covariate is a factor, the name of each list element 
#' refers to the level number rather than the label of the level. For each partition, \code{probMZ} also includes the number of labels 
#' for the partition (\code{M}) and the number of missing data (\code{NA}) in the outcome.
#'
#'
#' @export
CTEprobdata <- function( Yobs, Zobs, X=NULL, p = c(0:60)/60 , y=NULL , covariates = F, alpha=0.05, suppressMessages=FALSE) {
  if ( length(Yobs) != length(Zobs) ) stop("Yobs and Zobs must have the same length")
  
  if ( !is.null(X) ) {
    if ( is.null(dim(X)) ) stop("X must be a matrix (or a dataframe).")
    if ( nrow(X) != length(Yobs) ) stop("X must have the same number of rows that Yobs.")
    Xunique <- X |> unique()
  } else {
    if( covariates ) stop("X is NULL and covariates is TRUE.")
  }
  if ( any(is.character(p)) & ( length(p) > 1 | ! any(p == "PY1K1") ) ) {
    stop("Insert a valid value for p.")
  }
  
  if ( is.null(y) ) {
  } else {
    if( !suppressMessages ) { message(paste0("The outcome variable is 1 if Yobs > y, 0 if Yobs ≤ y, and NA if Yobs is NA. y=",y,".")) }
    Yobs <- 1*(Yobs > y)
  }
  
  nNAZ1 <- sum(is.na(Yobs[Zobs==1]))
  nNAZ0 <- sum(is.na(Yobs[Zobs==0]))
  
  if ( !covariates ) { 
    
    if ( length(p) > 1 ) {
      pp <- p
    } else {
      if ( p == "PY1K1" ) {
        pp <- mean(Yobs==1, na.rm = T)
      } else {
        pp <- p
      }
    }
    
    probs <- CTEprob( M = length(Yobs) , nNAZ1 = nNAZ1, nNAZ0 = nNAZ0, 
                      nZ1K1 = sum(Zobs[ !is.na(Yobs) ]),
                      nZ0K1 = sum(Zobs[ !is.na(Yobs) ] == 0),
                      nY1Z1K1 = sum(Yobs==1 & Zobs == 1 , na.rm = T), 
                      nY1Z0K1 = sum(Yobs==1 & Zobs == 0 , na.rm = T), 
                      p = pp , alpha=alpha )
    probs$data <- list( Yobs=Yobs, Zobs=Zobs, X=X, y=y)
    # gamma factor: proportion between CTEign and E(Pi)
    if ( p == "PY1K1" ) {
      n1 <- length(Yobs) - sum(Zobs[ !is.na(Yobs) ])
      n2 <- sum(Zobs[ !is.na(Yobs) ]) + (nNAZ1 + nNAZ0)
      probK1Z1 <- (length(Yobs) - n1) / length(Yobs)
      probK1Z0 <- (length(Yobs) - n2) / length(Yobs)
      probs$gamma <- 0.5 * ( 1/probK1Z1 + 1/probK1Z0 )
    }
  } else {
    
    probs <- list()
    for ( i in 1:nrow(Xunique) ) {
      filtro <- apply(X, 1, function(row) all(row == Xunique[i, ])) 
      Yaux <- Yobs[filtro]
      Zaux <- Zobs[filtro]
      Xaux <- X[filtro,]
      
      if ( p == "PY1K1" ) {
        pp <- mean(Yaux==1, na.rm = T)
      } else {
        pp <- p
      }
      
      probs[[i]] <- CTEprob( M = sum(filtro) , nNAZ1 = sum(is.na(Yaux[Zaux==1])),
                             nNAZ0 = sum(is.na(Yaux[Zaux==0])),
                             nZ1K1 = sum(Zaux[ !is.na(Yaux) ]),
                             nZ0K1 = sum(Zaux[ !is.na(Yaux) ] == 0),
                             nY1Z1K1 = sum(Yaux==1 & Zaux == 1 , na.rm = T), 
                             nY1Z0K1 = sum(Yaux==1 & Zaux == 0 , na.rm = T), 
                             p = pp  , alpha=alpha )
      names(probs)[i] <- paste0(paste0(colnames(X),".", Xunique[i, ]), collapse = ",")
      probs[[i]]$probMZ$M   <- sum(filtro)
      probs[[i]]$probMZ$nNAZ1 <- sum(is.na(Yaux[Zaux==1]))
      probs[[i]]$probMZ$nNAZ0 <- sum(is.na(Yaux[Zaux==0]))
      probs[[i]]$data <- list( Yobs=Yaux, Zobs=Zaux, X=Xaux, y=y)
      # gamma factor: proportion between CTEign and E(Pi)
      if ( p == "PY1K1" ) {
        n1 <- probs[[i]]$probMZ$M - sum(Zaux[ !is.na(Yaux) ])
        n2 <- sum(Zaux[ !is.na(Yaux) ]) + (probs[[i]]$probMZ$nNAZ1 + probs[[i]]$probMZ$nNAZ0)
        probK1Z1 <- (probs[[i]]$probMZ$M - n1) / probs[[i]]$probMZ$M
        probK1Z0 <- (probs[[i]]$probMZ$M - n2) / probs[[i]]$probMZ$M
        probs[[i]]$gamma <- 0.5 * ( 1/probK1Z1 + 1/probK1Z0 )
      }
    }
  }
  
  probs
}

