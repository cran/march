# Comparison indicators for markovian models.
# 
# Authors: Ogier Maitre, Andre Berchtold
###############################################################################

#' Compute Bayesian Information Criterion (BIC).
#' 
#' The BIC (Bayesian Information Criterion) is computed for
#'  a given \code{\link{march.Model-class}} according to the data
#'  used during construction.
#' 
#' @param model The model for which the BIC has to be computed.
#'
#' @return The number of parameters of the given model and its BIC.
#'
#' @author Ogier Maitre 
#' @example tests/examples/march.BIC.example.R
#' @export 
march.BIC <- function(model){
  nbParams <- march.nbParams(model)
  
	BIC <- -2*model@ll + nbParams*log(model@dsL)
	new("march.BIC",nbParams=nbParams,BIC=BIC)
	
}


#'	Compute Akaike Information Criterion (AIC).
#'  
#'  The AIC (Akaike Information Criterion) is computed for
#'  a given \code{\link{march.Model-class}} according to the data
#'  used during construction.
#'
#' @param model The model for which the AIC has to be computed.
#'
#' @return The number of parameters of the given model and its AIC.
#'
#' @author Ogier Maitre
#' @example tests/examples/march.AIC.example.R
#' @export
march.AIC <- function(model){
	nbParams <- march.nbParams(model)
	
	AIC <- 2*nbParams-2*model@ll
	new("march.AIC",nbParams=nbParams,AIC=AIC)
}