#' Construct a double chain Markov model (DCMM).
#'
#' Construct a \code{\link{march.Dcmm-class}} object, with visible order \emph{orderVC}, hidden order \emph{orderHC} and \emph{M} hidden states, according to a \code{\link{march.Dataset-class}}.
#' The first \emph{maxOrder}-\emph{orderVC} elements of each sequence are truncated in order to return a model
#' which can be compared with other Markovian model of visible order maxOrder. The construction is performed either by an evolutionary algorithm (EA) or by improving an existing DCMM.
#' The EA performs \emph{gen} generations on a population of \emph{popSize} individuals. The EA behaves as a Lamarckian evolutionary algorithm, using a Baum-Welch algorithm as
#' optimization step, running until log-likelihood improvement is less than \emph{stopBw} or for \emph{iterBw} iterations. Finally only the best individual from the population is returned as solution.
#' If a seedModel is provided, the only step executed is the optimization step, parameters related to the EA do not apply in this case.
#'
#' @param y the dataset from which the Dcmm will be constructed \code{\link{march.Dataset-class}}.
#' @param orderHC the order of the hidden chain of the constructed Dcmm.
#' @param orderVC the order of the visible chain of the constructed Dcmm (0 for a HMM).
#' @param M the number of hidden states of the Dcmm.
#' @param gen the number of generations performed by the EA.
#' @param popSize the number of individuals stored into the population.
#' @param maxOrder the maximum visible order among the set of Markovian models to compare.
#' @param seedModel a model to optimize using Baum-Welch algorithm.
#' @param iterBw the number of iterations performed by the Baum-Welch algorithm.
#' @param stopBw the minimum increase in quality (log-likelihood) authorized in the Baum-Welch algorithm.
#' @param Amodel the modeling of the hidden transition matrix (mtd, mtdg or complete)
#' @param Cmodel the modeling of the visible transition matrix (mtd, mtdg or complete)
#' @param AMCovar vector of the size Ncov indicating which covariables are used into the hidden process (0: no, 1:yes)
#' @param CMCovar vector of the size Ncov indicating which covariables are used into the visible process (0: no, 1:yes)
#' @param AConst logical, indicating whether or not the hidden transition matrix has the identity (diagonal) constraint
#' @param pMut mutation probability for the evolutionary algorithm
#' @param pCross crossover probability for the evolutionary algorithm
#'
#' @return the best \code{\link{march.Dcmm-class}} constructed by the EA or the result of the Baum-Welch algorithm on \emph{seedModel}.
#'
#' @author Emery Kevin
#' @example tests/examples/march.dcmm.construct.example.R
#' @seealso \code{\link{march.Dcmm-class}}, \code{\link{march.Model-class}}, \code{\link{march.Dataset-class}}.
#'
#' @export
march.dcmm.construct <- function(y,orderHC,orderVC,M=2,gen=5,popSize=4,maxOrder=orderVC,seedModel=NULL,iterBw=2,stopBw=0.1,Amodel="mtd",Cmodel="mtd",AMCovar=0,CMCovar=0, AConst=FALSE, pMut=0.05, pCross=0.5){
	
  if(is.logical(AConst)==FALSE){
    stop("AConst should be a logical")
  }
  
  if( is.null(seedModel) ){
    if(Amodel!="complete" & Amodel!="mtd" & Amodel!="mtdg"){
		  stop("Amodel should be equal to complete, mtd or mtdg")
	  }
	
	  if(Cmodel!="complete" & Cmodel!="mtd" & Cmodel!="mtdg"){
		  stop("Cmodel should be equal to complete mtd or mtdg")
	  }
	
	  if(length(AMCovar)!=y@Ncov & length(AMCovar)>1){
		  stop("AMCovar should have his length equal to the number of covariates in y")
	  }
	
	  if(length(CMCovar)!=y@Ncov & length(CMCovar)>1){
		  stop("CMCovar should have his length equal to the number of covariates in y")
	  }	
	
    if(AConst==TRUE & (Amodel=="mtd" | Amodel=="mtdg"|sum(AMCovar)>0)){
      stop("When the hidden transition matrix is constraint to the identity matrix, no mtd modeling can be used, nor covariates and the hidden order must be one")
    }
	
	  #Checking
	  if(M==1){
		  AMCovar <- rep(0,max(1,y@Ncov))
		  Amodel <- "complete"
		  orderHC <- 1
	  }
	
	  if(orderVC==0){
		  Cmodel <- "complete"
	  }
	
	  if(orderHC==0){
	  	stop("orderHC should be greater than 0")
		  #Amodel <- "complete"
	  }
	
	  if(orderHC==1 & Amodel=="mtdg"){
		  Amodel <- "mtd"
	  }
	
	  if(orderVC==1 & Cmodel=="mtdg"){
		  Cmodel <- "mtd"
	  }
		

		orderHC <- march.h.paramAsInteger(orderHC)
		orderVC <- march.h.paramAsInteger(orderVC)
		M <- march.h.paramAsInteger(M)
		maxOrder <- march.h.paramAsInteger(maxOrder)
    
    if( orderVC>maxOrder ){
      	stop("maxOrder should be greater or equal than orderVC")
    }
  }
  
  
  iterBw <- march.h.paramAsInteger(iterBw)
  
  
  gen <- march.h.paramAsInteger(gen)
  popSize <- march.h.paramAsInteger(popSize)
  
  if( is.null(seedModel) ){
    y <- march.dataset.h.filtrateShortSeq(y,maxOrder+1)
    y <- march.dataset.h.cut(y,maxOrder-orderVC)
  }
  
  if( is.null(seedModel)==FALSE ){
  	y <- march.dataset.h.filtrateShortSeq(y,maxOrder+1)
    y <- march.dataset.h.cut(y,maxOrder-seedModel@orderVC)
  	op <- new("march.dcmm.cov.ea.OptimizingParameters",fct=march.dcmm.cov.ea.optimizing,ds=y,stopBw=stopBw,iterBw=iterBw)
  	m <- march.dcmm.cov.ea.optimizing(seedModel,op)
  	
  	m@dsL <- sum(y@T-orderVC)
  	m@y <- y
  	m@ll <- march.dcmm.h.cov.computeLL(m,y)
  	m
  }else{
    ep <- new("march.dcmm.cov.ea.EvalParameters", ds=y, fct=march.dcmm.cov.ea.evaluation)
    ip <- new("march.dcmm.cov.ea.InitParameters",AConst=AConst,y=y,orderVC=orderVC,orderHC=orderHC,M=M,Amodel=Amodel,Cmodel=Cmodel,AMCovar=AMCovar,CMCovar=CMCovar,fct=march.dcmm.cov.ea.initialization)
    mp <- new("march.dcmm.cov.ea.MutationParameters",pMut=pMut,AConst=AConst,fct=march.dcmm.cov.ea.mutation)
    cp <- new("march.ea.cov.CrossoverParameters",fct=march.dcmm.cov.ea.crossover)
    
    op <- new("march.dcmm.cov.ea.OptimizingParameters",fct=march.dcmm.cov.ea.optimizing,ds=y,stopBw=stopBw,iterBw=iterBw)
    p <- new("march.ea.cov.Parameters",optimizing=TRUE,
             initParameters=ip,evalParameters=ep,mutationParameters = mp, optimizingParameters=op,crossoverParameters=cp,
             populationSize=popSize,crossoverProb=pCross,generation=gen)
    res <- march.loop.cov(p)
    
    res$best@dsL <- sum(y@T-orderVC)
    res$best@y <- y
    
    res$best@ll <- march.dcmm.h.cov.computeLL(res$best,y)
    res$best
    
  }
  
}


march.dcmm.cov.constructEmptyDcmm <- function(M,y,orderVC,orderHC,AMCovar,CMCovar,Amodel,Cmodel){
  
	KCovar <- y@Kcov
  	K <- y@K
  	NbAMCovar <- sum(AMCovar)
  	NbCMCovar <- sum(CMCovar)
  	placeACovar <- which(AMCovar==1)
   	placeCCovar <- which(CMCovar==1)
  
  	AtmCovar <- 1
  	if(NbAMCovar>0){
  		for (i in 1:NbAMCovar){
      		AtmCovar <- AtmCovar*KCovar[placeACovar[i]]
      	}
  	}
  
	  CtmCovar <- 1
  	if(NbCMCovar>0){
  		for (i in 1:NbCMCovar){
      		CtmCovar <- CtmCovar*KCovar[placeCCovar[i]]
    	}
  	}
 	
 	if(orderHC==0){
 		  Pi <- array(0,c(AtmCovar,M,1))
 	}else{
  		Pi <- array(0,c(M^(orderHC-1)*AtmCovar,M,orderHC))
 	}
  	
  	if(Amodel=="complete"){
    	AQ <- array(0,c(1,M^max(orderHC,1),M))
  	}else if(Amodel=="mtd"){
    	AQ <- array(0,c(1,M,M))
  	}else{
    	AQ <- array(0,c(orderHC,M,M))
  	}
  
 	if(Cmodel=="complete"){
    	CQ <- array(0,c(1,K^orderVC,K,M))
  	}else if(Cmodel=="mtd"){
    	CQ <- array(0,c(1,K,K,M))
  	}else{
    	CQ <- array(0,c(orderVC,K,K,M))
  	}
  
  	ATCovar <- list()
  	CTCovar <- list()
  
  	if(NbAMCovar>0){
  		for(i in 1:NbAMCovar){
  			ATCovar[[i]] <- array(0,c(KCovar[placeACovar[i]],M))
  		}
  	}
  
  	if(NbCMCovar>0){
  		for(i in 1:NbCMCovar){
  			CTCovar[[i]] <- array(0,c(KCovar[placeCCovar[i]],K,M))
  		}
  	}
  
  	if(Amodel=="complete"){
    	APhi <- array(1,c(1,1+NbAMCovar))
  	}else{
    	APhi <- array(1,c(1,orderHC+NbAMCovar))
  	}
  	if(Cmodel=="complete"){
      	CPhi <- array(1,c(1,1+NbCMCovar,M))
  	}else{
    	CPhi <- array(1,c(1,orderVC+NbCMCovar,M))
  	}
  
  	if(Amodel=="complete"){
  		AProbT <- array(0,c(M^(max(orderHC,1)+1)*AtmCovar,1+NbAMCovar))
  	}else{
  		AProbT <- array(0,c(M^(max(orderHC,1)+1)*AtmCovar,orderHC+NbAMCovar))
  	}
  
  	if(Cmodel=="complete"){
  		CProbT <- array(0,c(K^(orderVC+1)*CtmCovar,1+NbCMCovar,M))
  	}else{
  		CProbT <- array(0,c(K^(orderVC+1)*CtmCovar,orderVC+NbCMCovar,M))
  	}
  
  	A <- array(1,c(M^orderHC*AtmCovar,M^orderHC))
  	RB <- array(1,c(K^orderVC*CtmCovar,K,M))
  
  
	new("march.Dcmm",Pi=Pi,A=A,M=M,orderVC=orderVC,orderHC=orderHC,RB=RB,y=y,APhi=APhi,CPhi=CPhi,ATCovar=ATCovar,CTCovar=CTCovar,AMCovar=AMCovar,CMCovar=CMCovar,
  AQ=AQ, CQ=CQ, Amodel=Amodel,Cmodel=Cmodel,AProbT=AProbT,CProbT=CProbT)
}


march.dcmm.cov.h.compactA <- function(d){
	
	NbAMCovar <- sum(d@AMCovar)  
  	placeACovar <- which(d@AMCovar==1)
   
   	AtmCovar <- 1
  	if(NbAMCovar>0){
  		for (i in 1:NbAMCovar){
      		AtmCovar <- AtmCovar*d@y@Kcov[placeACovar[i]]
      	}
  	}

	
  	RA <- array(0,c(d@M^d@orderHC*AtmCovar,d@M))
  
	if(d@orderHC==0){
		RA <- d@A
	}else{
  		for( r2 in 1:(d@M*AtmCovar )){
    		for( r1 in 1:d@M^(d@orderHC-1) ){
      			for( c in 1:d@M ){
        			RA[(r1-1)*(d@M*AtmCovar)+r2,c] <- d@A[(r1-1)*(d@M*AtmCovar)+r2,r1+(c-1)*d@M^(d@orderHC-1)]
        		}
      		}
    	}
  	}
  
  	RA
}

march.dcmm.cov.h.expandRA <- function(d,RA){
    NbAMCovar <- sum(d@AMCovar)  
   	placeACovar <- which(d@AMCovar==1)
   
   	AtmCovar <- 1
   	if(NbAMCovar>0){
   		for (i in 1:NbAMCovar){
      		AtmCovar <- AtmCovar*d@y@Kcov[placeACovar[i]]
      	}
  	}

  	A <- array(0,c(d@M^d@orderHC*AtmCovar,d@M^d@orderHC))
  	offset <- 0
   	for(r1 in 1:(d@M^(d@orderHC-1))){
   		for(r2 in 1:d@M){
   			for(off in 1:AtmCovar){
   				offset <- offset+1
   				for(c in 1:d@M){
   					A[offset,r1+(c-1)*d@M^(d@orderHC-1)] <- RA[offset,c]
   				}
   			}
   		}
   	}
   	A

}

#' Viterbi algorithm for a DCMM model.
#'
#' @param d The \code{\link[=march.Dcmm-class]{march.Dcmm-class}} on which to compute the most likely sequences of hidden states.
#'
#' @return A list of vectors containing the most likely sequences of hidden states, considering the given model for each sequence of the given dataset.
#' @author Kevin Emery
#' @example tests/examples/march.dcmm.viterbi.example.R
#'
#'@export
march.dcmm.viterbi <- function(d){
  
  placeACovar <- which(d@AMCovar==1)
  placeCCovar <- which(d@CMCovar==1)
  NbAMCovar <- sum(d@AMCovar)
  NbCMCovar <- sum(d@CMCovar)
  
  AtmCovar <- 1
  if(NbAMCovar>0){
    for (i in 1:NbAMCovar){
      AtmCovar <- AtmCovar*d@y@Kcov[placeACovar[i]]
    }
  }
  
  
  CtmCovar <- 1
  if(NbCMCovar>0){
    for (i in 1:NbCMCovar){
      CtmCovar <- CtmCovar*d@y@Kcov[placeCCovar[i]]
    }
  }
  
  l <- list()
  
  for(n in 1:d@y@N){
    
    #Initialization
    s <- march.dataset.h.extractSequence(d@y,n)
    delta <- matrix(0,d@M^d@orderHC,s@N)
    avdel <- rep(0,s@N)
    X <- rep(0,s@N)
    
    t <- d@orderVC+1
    
    pos <- s@y[(t-d@orderVC):(t-1)]
    rowC <- 1
    for(r in march.h.seq(1,d@orderVC)){
      rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
    }
    
    rowCovar <- 1
    if(AtmCovar>1){
      r <- 0
      tKC <- 1
      for(i in d@y@Ncov:1){
        if(d@AMCovar[i]>0){
          r <- r+1
          rowCovar <- rowCovar+tKC*(d@y@cov[n,t,i]-1)
          tKC <- tKC*d@y@Kcov[i]
        }
      }
    }
    
    rowCovarC <- 1
    if(CtmCovar>1){
      r <- 0
      tKC <- 1
      for(i in d@y@Ncov:1){
        if(d@CMCovar[i]>0){
          r <- r+1
          rowCovarC <- rowCovarC +tKC*(d@y@cov[n,t,i]-1)
          tKC <- tKC*d@y@Kcov[i]
        }
      }
    }
    
    for(i in 1:d@M){
      delta[i,t] <- d@Pi[rowCovar,i,1]*d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t],i]
    }
    
    avdel[t] <- sum(delta[1:d@M,t])/d@M
    
    if(avdel[t]==0){
      delta[1:d@M,t] <- rep(1,d@M)/d@M
    }else{
      delta[1:d@M,t] <- delta[1:d@M,t]/avdel[t]
    }
    
    #Induction
    
    #t=orderVC+2,...,orderVC+orderHC
    
    if(d@orderHC>1){
      for(t in (d@orderVC+2):(d@orderVC+d@orderHC)){
        
        pos <- s@y[(t-d@orderVC):(t-1)]
        rowC <- 1
        for(r in march.h.seq(1,d@orderVC)){
          rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
        }
        
        rowCovar <- 1
        if(AtmCovar>1){
          r <- 0
          tKC <- 1
          for(i in d@y@Ncov:1){
            if(d@AMCovar[i]>0){
              r <- r+1
              rowCovar <- rowCovar+tKC*(d@y@cov[n,t,i]-1)
              tKC <- tKC*d@y@Kcov[i]
            }
          }
        }
        
        rowCovarC <- 1
        if(CtmCovar>1){
          r <- 0
          tKC <- 1
          for(i in d@y@Ncov:1){
            if(d@CMCovar[i]>0){
              r <- r+1
              rowCovarC <- rowCovarC +tKC*(d@y@cov[n,t,i]-1)
              tKC <- tKC*d@y@Kcov[i]
            }
          }
        }
        
        for(j in 1:d@M^(t-d@orderVC)){
          j0 <- floor((j-1)/(d@M^(t-d@orderVC-1)))+1
          
          calc <- rep(0,d@M^(t-d@orderVC-1))
          for(k in 1:d@M^(t-d@orderVC-1)){
            calc[k] <- delta[k,t-1]*d@Pi[(k-1)*AtmCovar+rowCovar,j0,t-d@orderVC]
          }
          calc2 <- max(calc)
          
          delta[j,t] <- d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t],j0]*calc2
        }
        avdel[t] <- sum(delta[1:d@M^(t-d@orderVC)],t)/(d@M^(t-d@orderVC))
        
        if(avdel[t]==0){
          delta[1:d@M^(t-d@orderVC),t] <- rep(1,d@M^(t-d@orderVC))/(d@M^(t-d@orderVC))
        }else{
          delta[1:d@M^(t-d@orderVC),t] <- delta[1:d@M^(t-d@orderVC),t]/avdel[t]
        }
      }
    }
    
    #t=orderVC+orderHC+1,...,s@N
    
    for(t in march.h.seq(d@orderVC+d@orderHC+1,s@N)){
      
      pos <- s@y[(t-d@orderVC):(t-1)]
      rowC <- 1
      for(r in march.h.seq(1,d@orderVC)){
        rowC <- rowC+d@y@K^(r-1)*(pos[r]-1)
      }
      
      rowCovar <- 1
      if(AtmCovar>1){
        r <- 0
        tKC <- 1
        for(i in d@y@Ncov:1){
          if(d@AMCovar[i]>0){
            r <- r+1
            rowCovar <- rowCovar+tKC*(d@y@cov[n,t,i]-1)
            tKC <- tKC*d@y@Kcov[i]
          }
        }
      }
      
      rowCovarC <- 1
      if(CtmCovar>1){
        r <- 0
        tKC <- 1
        for(i in d@y@Ncov:1){
          if(d@CMCovar[i]>0){
            r <- r+1
            rowCovarC <- rowCovarC +tKC*(d@y@cov[n,t,i]-1)
            tKC <- tKC*d@y@Kcov[i]
          }
        }
      }
      
      for(j in 1:d@M^d@orderHC){
        j0 <- floor((j-1)/(d@M^(d@orderHC-1)))+1
        
        calc <- rep(0,d@M^d@orderHC)
        for(k in 1:d@M^d@orderHC){
          calc[k] <- delta[k,t-1]*d@A[(k-1)*AtmCovar+rowCovar,j]
        }
        calc2 <- max(calc)
        
        delta[j,t] <- d@RB[(rowC-1)*CtmCovar+rowCovarC,s@y[t],j0]*calc2
      }
      avdel[t] <- sum(delta[,t])/(d@M^d@orderHC)
      
      if(avdel[t]==0){
        delta[,t] <- rep(1,d@M^d@orderHC)/(d@M^d@orderHC)
      }else{
        delta[,t] <- delta[,t]/avdel[t]
      }
    }
    
    #Sequence of hidden states
    
    calc4 <- max(delta[,s@N])
    calc5 <- which(calc4==delta[,s@N])
    calc6 <- calc5[1]
    calc7 <- floor((calc6-1)/(d@M^(d@orderHC-1)))+1
    X[s@N] <- calc7
    
    for(t in march.h.seq(s@N-1,d@orderHC+d@orderVC,-1)){
      
      rowCovar <- 1
      if(AtmCovar>1){
        r <- 0
        tKC <- 1
        for(i in d@y@Ncov:1){
          if(d@AMCovar[i]>0){
            r <- r+1
            rowCovar <- rowCovar+tKC*(d@y@cov[n,t,i]-1)
            tKC <- tKC*d@y@Kcov[i]
          }
        }
      }
      
      calc <- rep(0,d@M^d@orderHC)
      for(k in 1:d@M^d@orderHC){
        calc[k] <- delta[k,t]*d@A[(k-1)*AtmCovar+rowCovar,calc6]
      }
      
      calc4 <- max(calc)
      calc5 <- which(calc==calc4)
      calc6 <- calc5[1]
      calc7 <- floor((calc6-1)/(d@M^(d@orderHC-1)))+1
      X[t] <- calc7
    }
    
    for(t in march.h.seq(d@orderVC+d@orderHC-1,d@orderVC+1,-1)){
      
      rowCovar <- 1
      if(AtmCovar>1){
        r <- 0
        tKC <- 1
        for(i in d@y@Ncov:1){
          if(d@AMCovar[i]>0){
            r <- r+1
            rowCovar <- rowCovar+tKC*(d@y@cov[n,t,i]-1)
            tKC <- tKC*d@y@Kcov[i]
          }
        }
      }
      
      calc <- rep(0,d@M^(t-d@orderVC))
      for(k in 1:d@M^(t-d@orderVC)){
        calc[k] <- delta[k,t]*d@Pi[(k-1)*AtmCovar+rowCovar,calc7,t-d@orderVC+1]
      }
      
      calc4 <- max(calc)
      calc5 <- which(calc==calc4)
      calc6 <- calc5[1]
      calc7 <- floor((calc6-1)/(d@M^(t-d@orderVC-1)))+1
      X[t] <- calc7
    }
    l[[n]] <- X
  }
  l
}