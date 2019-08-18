#
# Author: Ogier Maitre
###############################################################################


#' Construct an independence model (zero-order Markov chain).
#' 
#' Construct a \code{\link{march.Indep-class}} model from a given \code{\link{march.Dataset-class}},
#' the first \emph{maxOrder} elements of each sequence being truncated in order to return a model
#' which can be compared with other Markovian models of visible order maxOrder. 
#' 
#' @param y the \code{\link{march.Dataset-class}} from which construct the model.
#' @param maxOrder the maximum visible order among the set of Markovian models to compare.
#' 
#' @return The \code{\link{march.Indep-class}} constructed using dataset y and maxOrder.
#' 
#' 
#' @seealso \code{\link{march.Indep-class}}, \code{\link{march.Model-class}}, \code{\link{march.Dataset-class}}.
#' @example tests/examples/march.indep.construct.example.R
#' @author Ogier Maitre
#' @export
march.indep.construct <- function(y,maxOrder=0){
  maxOrder <- march.h.paramAsInteger(maxOrder)
	indC <- matrix(0,1,y@K);
	
	# numbering all the case
	for( i in 1:y@N ){
		if( y@T[i]>=0 ){ 
			for( j in march.h.seq(maxOrder+1,y@T[i])){
        if( is.na(y@y[[i]][j])==FALSE ){ indC[y@y[[i]][j]]=indC[y@y[[i]][j]]+y@weights[i]; }
			}
		}
	}
	
	total <- sum(indC);
	indP <- indC/total;
	
	ll <- 0;
	for( i in 1:y@K ){
		if( indC[i]>0 ){
			ll=ll+indC[i]*log(indP[i]);
		}
	}
	
	new("march.Indep",indP=indP,indC=indC,ll=ll,y=y,dsL=sum(indC))
}