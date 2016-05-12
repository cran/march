# Performance indicators for markovian models.
# 
# Author: Ogier Maitre
###############################################################################

#' Compute Bayesian Information Criterion (BIC).
#' 
#' The BIC (Bayesian Information Criterion) is computed for
#'  a given \code{\link{march.Model-class}} according to the data
#'  used during construction.
#' 
#' @param model The model for which the BIC has to be computed.
#'
#' @author Ogier Maitre 
#' @example tests/examples/march.BIC.example.R
#' @export 
march.BIC <- function(model){
	-2*model@ll + march.nbParams(model)*log(model@dsL)
}


#'	Compute Akaike Information Criterion (AIC).
#'  
#'  The AIC (Akaike Information Criterion) is computed for
#'  a given \code{\link{march.Model-class}} according to the data
#'  used during construction.
#'
#' @param model The model for which the AIC has to be computed.
#'
#' @return The AIC of the given model.
#'
#' @author Ogier Maitre
#' @example tests/examples/march.AIC.example.R
#' @export
march.AIC <- function(model){
	nbParams <- march.nbParams(model)
	
	2*nbParams-2*model@ll
}