# TODO: Add comment
# 
# Author: SSP24298
###############################################################################

# Initialize a random genome
# This function is called to populate the initial parent population only.
# 
# Parameter: 
#   p : an initialization parameter (see march.AllClasses.R) havine type march.dcmm.ea.InitParameters.
# Return:
#   The newly created individual.
#
march.dcmm.ea.initialization <- function(p){
	if( p@AConst ){
		print("Error, this part is not yet implemented")
	}
	else{
		d <- march.dcmm.constructEmptyDcmm(p@M,p@y,p@orderVC,p@orderHC)
		
		# randomly initialize the A matrix
		RA <- matrix(runif(d@M^d@orderHC*d@M),nrow=d@M^d@orderHC,ncol=d@M)
		RA <- RA/rowSums(RA)
		d@A <- march.dcmm.h.expandRA(d,RA)
			
		# randomly initialize the B 3D matrix
		RB <- array(runif(d@M*d@y@K^d@orderVC*d@y@K),c(d@M,p@y@K^d@orderVC,p@y@K))
   
    for( i in 1:d@M ){
      if( p@orderVC==0 ){ sRB <- sum(RB[i,1,]) }
      else{ sRB <- rowSums(RB[i,,]) }
      
      RB[i,,] <- RB[i,,]/sRB
    }
			
		d@RB <- RB # which is stored as a reduced form into the model
		
		#randomly initialize the Pi array
		Pi <- array(0,c(d@orderHC,d@M^(d@orderHC-1),d@M))
		
		for( i in march.h.seq(1,d@orderHC) ){
			Pi[i,1:(d@M^(i-1)),] <- 0.1+matrix(runif(d@M^(i-1)*d@M),nrow=d@M^(i-1))
			if( i==1 ){
				Pi[i,1:(d@M^(i-1)),] <- Pi[i,1:(d@M^(i-1)),]/sum(Pi[i,1:(d@M^(i-1)),])
			}
			else{
				Pi[i,1:(d@M^(i-1)),] <- Pi[i,1:(d@M^(i-1)),]/rowSums(Pi[i,1:(d@M^(i-1)),])
			}
		}
		d@Pi <- Pi
	}
	d
}

# Compute fitness value, based on Log likelihood (minus the invert)
# TODO : handle multi-sequence case
#
# Parameters:
#	  d: dcmm to evaluate
#	  p: the evaluate parameters, see march.dcmm.ea.EvalParameters in march.AllClasses.R file.
# Return:
#	  the opposite of the invert of the log likelihood
#
march.dcmm.ea.evaluation <- function(d,p){
	d@ll <- march.dcmm.h.computeLL(d=d,y=p@ds)
  d@dsL <- sum(p@ds@T)
  
	-1/d@ll
}

# This function implements crossover for dcmm objects.
# It uses sbx applied to every row of A,C and Pi.
#
# Paramters:
#	  d1: first parent
#	  d2: second parent
# return: (c1,c2)
#	  c1: first child, produced by crossover
# 	c2: second child, produced by crossover.
#
march.dcmm.ea.crossover <- function(d1, d2){
	
	c1 <- march.dcmm.constructEmptyDcmm(d1@M,d1@y,d1@orderVC,d1@orderHC)			
	c2 <- march.dcmm.constructEmptyDcmm(d1@M,d1@y,d1@orderVC,d1@orderHC)
	
	# compute reduce transition matrix for parent 1 and 2
	d1RA <- march.dcmm.h.compactA(d1)
	d2RA <- march.dcmm.h.compactA(d2)
	
	c1RA <- array(0,dim(d1RA))
	c2RA <- array(0,dim(d2RA))
	
	# crossover parent 1 and 2 reduce transition matrix
	for( i in 1:d1@M^d1@orderHC ){
		r <- march.ea.h.sbx(d1RA[i,],d2RA[i,],0.1)
		c1RA[i,] <- r$c1
		c2RA[i,] <- r$c2
	}
	
	# expand into standard forms
	c1@A <- march.dcmm.h.expandRA(c1,c1RA)
	c2@A <- march.dcmm.h.expandRA(c2,c2RA)
		
	# crossover emission matrix
	for( i in 1:d1@M ){
		for( j in 1:d1@y@K^d1@orderVC ){
			r <- march.ea.h.sbx(d1@RB[i,j,],d2@RB[i,j,],2)
			c1@RB[i,j,] <- r$c1
			c2@RB[i,j,] <- r$c2
		}
	}
	
	# crossover the initial probability
	for( i in 1:d1@orderHC ){
		for( j in 1:(d1@M^(i-1))){
			r <- march.ea.h.sbx(d1@Pi[i,j,],d2@Pi[i,j,],2)
			c1@Pi[i,j,] <- r$c1
			c2@Pi[i,j,] <- r$c2
		}
	}
	
	# return both children
	list(c1=c1,c2=c2)
}

# Mutation for dcmm. It is a gaussian mutation applied to each gene, 
# with a probability pMut.
#
# Parameters:
#	  d: dcmm to mutate
#	  p: corresponding mutation parameter object, see march.dcmm.ea.MutationParameters in march.AllClasses.R
#
# Return: the mutated individual
#
march.dcmm.ea.mutation <- function(d,p){
	
	RA <- march.dcmm.h.compactA(d)
	for( i in 1:d@M^d@orderHC ){
		RA[i,] <- RA[i,]+(runif(d@M)<p@pMut)*rnorm(d@M)
		RA[i,] <- march.ea.checkBounds(v=RA[i,],lb=0,ub=1)
		RA[i,] <- march.dcmm.ea.h.scale(RA[i,],d@M)
	}
	d@A <- march.dcmm.h.expandRA(d,RA)
	
	for( i in d@M ){
		for( j in d@y@K^d@orderVC ){
			d@RB[i,j,] <- d@RB[i,j,]+(runif(d@y@K)<p@pMut)*rnorm(d@y@K)
			d@RB[i,j,] <- march.ea.checkBounds(v=d@RB[i,j,],lb=0,ub=1)
			d@RB[i,j,] <- march.dcmm.ea.h.scale(d@RB[i,j,],d@y@K)
		}
	}
	
	for( i in 1:d@orderHC ){
		for( j in 1:(d@M^(i-1))){
			d@Pi[i,j,] <- d@Pi[i,j,]+(runif(d@M)<p@pMut)*rnorm(d@M)
			d@Pi[i,j,] <- march.ea.checkBounds(v=d@Pi[i,j,],lb=0,ub=1)
			d@Pi[i,j,] <- march.dcmm.ea.h.scale(d@Pi[i,j,],d@M)
		}
	}
	d
}


march.dcmm.ea.h.scale <- function(v,j){
	mA <- sum(v)
	
	if( mA==0 ){ v <- matrix(1,ncol=1,nrow=j)*(1/j) }
	else{ v <- v/mA }
	v
}


# Compute the log-likelyhood for a multi-sequence discrete-valued series.
march.dcmm.h.computeLL <- function(d,y){
  ll <- 0
  
  for( i in 1:y@N){ 
    if( y@T[i]>(d@orderVC+d@orderVC) ){
      seq <- march.dataset.h.extractSequence(y,i)
      ll <- ll+march.dcmm.forward(d,seq)$LL 
    }
  }
  
  ll
}

# Optimizing function. It first constraints each row of each element to be sumed to 1,
# then it applies Baum-Welch algorithm.
#
# Parameters:
#	  d: dcmm to optimize
#	  p: corresponding optimizing parameter object, see march.dcmm.ea.OptimizingParameters in march.AllClasses.R
#
# Return: the optimized dcmm
#
march.dcmm.ea.optimizing <- function(d,p){
	# maintain constraints on each attribute for child1 and child2
	for( i in 1:d@M^d@orderHC ){ 
    if( sum(d@A[i,])!=0 ){
      d@A[i,] <- d@A[i,]/sum(d@A[i,])       
    }
	}
	for( i in d@M ){ 
    for( j in d@y@K^d@orderVC ){ 
      if( sum(d@RB[i,j,])!=0 )
        d@RB[i,j,] <- d@RB[i,j,]/sum(d@RB[i,j,]) 
	  }
	}
  
	for( i in 1:d@orderHC ){ 
    for( j in 1:(d@M^(i-1))){ 
      if( sum(d@Pi[i,j,])!=0 ){
        d@Pi[i,j,] <- d@Pi[i,j,]/sum(d@Pi[i,j,]) 
      }
	  } 
	}	
  
  opt <- d  # the optimized individual.
  
  referenceLL <- march.dcmm.h.computeLL(opt,p@ds)
  for( i in march.h.seq(1,p@iterBw)){
    opt <- march.dcmm.bw(opt,p@ds)
    
    currentLL <- march.dcmm.h.computeLL(opt,p@ds)
    if( abs(currentLL-referenceLL)<=p@stopBw ){ break }
  }
  opt
}