#
# Author: ogiermaitre
###############################################################################

#' Construct an homogeneous Markov Chain.
#'
#' A \code{\link[=march.Mc-class]{march.Mc}} object of order \emph{order} is constructed from
#' the dataset \emph{y}. The first maxOrder-order elements of each sequence of the dataset are truncated in order to return a model
#' which can be compared with other Markovian models of visible order maxOrder.
#'
#' @param y the \code{\link[=march.Dataset-class]{march.Dataset}} from which the homogeneous Markov chain will be constructed.
#' @param order the order of the constructed Markov Chain.
#' @param maxOrder the maximum visible order among the set of Markovian models to compare.
#'
#' @return the \code{\link[=march.Mc-class]{march.Mc}} of order \emph{order} constructed \emph{w.r.t} the dataset \emph{y} and maxOrder.
#'
#' @author Ogier Maitre
#' @example tests/examples/march.mc.construct.example.R
#' @seealso \code{\link[=march.Mc-class]{march.Mc-class}}, \code{\link{march.Model-class}}, \code{\link[=march.Dataset-class]{march.Dataset-class}}.
#' @export
march.mc.construct <- function(y,order,maxOrder=order) {

  order <- march.h.paramAsInteger(order)
  maxOrder <- march.h.paramAsInteger(maxOrder)

  if( order>maxOrder ){
    stop("maxOrder should be greater or equal than order")
  }

	rt <- array(0,c(y@K^order,y@K))
	for(i in 1:y@N){
		if(y@T[i]>order){
			for(t in march.h.seq(maxOrder-order+1,y@T[i]-order)){
				past <- y@y[[i]][t:(t+order-1)]
				present <- y@y[[i]][t+order]
				row <- past[1]
				if(order>1){
					for(g in 2:order){
						row <- row+y@K^(g-1)*(past[g]-1)
					}
				}
				rt[row,present] <- rt[row,present] + y@weights[i]
			}
		}
	}

	tot <- rt%*%array(1,c(y@K,1))
	r <- array(0,c(y@K^order,y@K))
	for(i in 1:y@K^order){
		if(tot[i]!=0){
			r[i,] <- rt[i,]/tot[i]
		}
	}

	ll <- 0
	for(i in 1:(y@K^order)){
		for(j in 1:y@K){
			if(rt[i,j]!=0){
				ll <- ll + rt[i,j]*log(r[i,j])
			}
		}
	}

	new("march.Mc",RC=r,order=order,ll=ll,RT=rt,y=y,dsL=sum(rt),nbZeros=length(which(rt==0)))
}


#march.mc.bailey <- function( indep, y, alpha ){
#
#	n <- sum(y@T)
#	n_i <- rowSums(mc@RT)
#
#	A <- march.ci.h.A(n,n_i)
#	B <- march.ci.h.B(n,n_i)
#	C <- march.ci.h.C(alpha,y@K,n)
#
#	pm <- array(NA,c(length(n_i)))
#	pp <- array(NA,c(length(n_i)))
#
#	for( i in 1:length(n_i)){
#		if( n_i[i]>0 ){
#			pm[i] <- (sqrt(A[i])-sqrt(C*(C+1-A[i])))^2/(C+1)^2
#			pp[i] <- (sqrt(B[i])+sqrt(C*(C+1-B[i])))^2/(C+1)^2
#		}
#	}
#
#	list(pm=pm,pp=pp)
#}
