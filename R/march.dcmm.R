#
# Author: Ogier Maitre
###############################################################################

march.dcmm.constructEmptyDcmm <- function(M,y,orderVC,orderHC){
  RB <- array(1/y@K,c(M,y@K^orderVC,y@K)) # emit
  Pi <- array(1/M,c(orderHC,M^(orderHC-1),M))

  RA <- array(1/M,c(M^orderHC,M))
  A <- march.dcmm.h.expandRAInternal(RA,M,orderHC)

  new("march.Dcmm",Pi=Pi,A=A,M=M,orderVC=orderVC,orderHC=orderHC,RB=RB,y=y)
}


march.dcmm.h.encodeOutput <- function(s,t,d){

  if( d@orderVC>0 ){
    code <- array(0,c(d@orderVC-1))
    i <- march.h.seq(0,d@orderVC-1)
    code = d@y@K^(i)*(s@y[i+(t-d@orderVC)]-1)
    sum(code)+1
  }
  else{
    1
  }
}

march.dcmm.h.scaleAlpha <- function(j,alphat){
  mA <- sum(alphat[1:j])/j
  if( mA==0 ){
    alphat[1:j] <- matrix(1,ncol=1,nrow=j)*(1/j)
    l <- 1
  }
  else{
    alphat[1:j] <- alphat[1:j]/mA
    l <- mA
  }
  list(scaled=alphat,scale=l)
}


march.dcmm.h.expandRA <- function(d,RA){
  march.dcmm.h.expandRAInternal(RA,d@M,d@orderHC)
}

march.dcmm.h.expandRAInternal <- function(RA,M,orderHC){
  A <- array(0,c(M^orderHC,M^orderHC))

  for( r2 in 1:M ){
    for( r1 in march.h.seq(1,M^(orderHC-1)) ){
      for( c in 1:M ){
        A[(r1-1)*M+r2,r1+(c-1)*M^(orderHC-1)] <- RA[(r1-1)*M+r2,c]
      }
    }
  }

  A
}

march.dcmm.h.compactA <- function(d){
  RA <- array(0,c(d@M^d@orderHC,d@M))

  for( r2 in 1:d@M ){
    for( r1 in 1:d@M^(d@orderHC-1) ){
      for( c in 1:d@M ){
        RA[(r1-1)*d@M+r2,c] <- d@A[(r1-1)*d@M+r2,r1+(c-1)*d@M^(d@orderHC-1)]
      }
    }
  }

  RA
}


#' Viterbi algorithm for a DCMM model.
#'
#' @param d The \code{\link[=march.Dcmm-class]{march.Dcmm-class}} on which to compute the most likely sequences of hidden states.
#' @param y The \code{\link[=march.Dataset-class]{march.Dataset-class}} to consider.
#'
#' @return A list of vectors containing the most likely sequences of hidden states, considering the given model for each sequence of the given dataset.
#' @author Ogier Maitre
#' @example tests/examples/march.dcmm.viterbi.example.R
#'
#'@export
march.dcmm.viterbi <- function(d,y){

  l <- list()

  for( n in 1:y@N ){
    s <- march.dataset.h.extractSequence(y,n)
    delta <- matrix(0,nrow=s@N,ncol=d@M^d@orderHC)

    # intial step t=orderVC+1
    t = d@orderVC+1
    for( j0 in 1:d@M ){
      delta[t,j0] <- d@Pi[1,1,j0]*d@RB[j0,march.dcmm.h.encodeOutput(s,t,d),s@y[t]]
    }

    delta[t,] <- march.dcmm.h.scaleAlpha(d@M,delta[t,])$scaled

    if( d@orderHC>1 ){
      # initial step to<=orderVC+orderHC
      for( t in march.h.seq(d@orderVC+2,d@orderVC+d@orderHC)){
        rowC <- march.dcmm.h.encodeOutput(s,t,d)
        for( j in 1:(d@M^(t-d@orderVC)) ){
          j0 = floor((j-1)/(d@M^(t-d@orderVC-1)))+1;
          delta[t,j] <- max(delta[t-1,]*d@Pi[t-d@orderVC,,j0])*d@RB[j0,rowC,s@y[t]]
        }
        delta[t,] <- march.dcmm.h.scaleAlpha(d@M^(t-d@orderVC),delta[t,])$scaled
      }
    }

    # Induction to t=length(y)
    for( t in (d@orderVC+d@orderHC+1):s@N){ #
      rowC <- march.dcmm.h.encodeOutput(s,t,d)
      for( j in 1:(d@M^d@orderHC)){
        j0 = floor((j-1)/(d@M^(d@orderHC-1)))+1;
        delta[t,j] <- max(delta[t-1,]*d@A[,j])*d@RB[j0, rowC,s@y[t]]
      }
      delta[t,] <- march.dcmm.h.scaleAlpha(d@M^d@orderHC,delta[t,])$scaled
    }


    X <- vector(mode="integer",s@N)
    e_t = which.max(delta[t,])
    x_t = floor((e_t-1)/(d@M^(d@orderHC-1)))+1

    X[s@N] <- x_t

    # find the best sequence of hidden step
    for( t in seq(s@N-1,d@orderVC+d@orderHC) ){
      e_t <- which.max(delta[t,]*d@A[,e_t])
      x_t <- floor((e_t-1)/(d@M^(d@orderHC-1)))+1
      X[t] <- x_t
    }

    if (d@orderHC>1){
      for( t in seq(d@orderHC+d@orderVC-1,d@orderVC+1)){
        e_t <- which.max(delta[t,]*d@Pi[t-d@orderVC+1,x_t,])
        X[t] <- floor((e_t-1)/(d@M^(d@orderHC-1)))+1
      }
    }

    l[[n]] <- X
  }
  l
}

#
# parameters :
#	d : a Dcmm object (see Dcmm class)
#	y : a vector containing variable sequence.
#
# return :
#
march.dcmm.forward <- function(d,s){
  alpha <- matrix(0,nrow=s@N,ncol=d@M^d@orderHC)
  l <- vector("numeric",s@N)

  # intial step t=orderVC+1
  t <- d@orderVC+1;
  alpha[t,1:d@M] <- d@RB[,march.dcmm.h.encodeOutput(s,t,d),s@y[t]]*d@Pi[1,1,]

  scaled <- march.dcmm.h.scaleAlpha(d@M,alpha[t,])
  alpha[t,]<- scaled$scaled
  l[t] <- log(scaled$scale)

  if( d@orderHC>1 ){
    # initial step to<=orderVC+orderHC
    for( t in march.h.seq((d@orderVC+2),(d@orderVC+d@orderHC)) ){
      rowC <- march.dcmm.h.encodeOutput(s,t,d)
      for( i in march.h.seq(1,d@M^(t-d@orderVC-1)) ){
        j<-1:d@M
        j0 <- d@M^(t-d@orderVC-1)*(j-1)+i
        alpha[t,j0] <- d@RB[j,rowC,s@y[t]]*d@Pi[t-d@orderVC,i,j]*alpha[t-1,i]
      }
      scaled <- march.dcmm.h.scaleAlpha(d@M^(t-d@orderVC),alpha[t,])
      alpha[t,] <- scaled$scaled
      l[t] <- log(scaled$scale)
    }
  }

  # Induction to t=length(y)
  for( t in march.h.seq(d@orderVC+d@orderHC+1,s@N) ){
    rowC <- march.dcmm.h.encodeOutput(s,t,d)
    for( j in 1:(d@M^d@orderHC)){
      j0 = floor((j-1)/(d@M^(d@orderHC-1)))+1;
      alpha[t,j] <- t(alpha[t-1,])%*%d@A[,j]*d@RB[j0, rowC,s@y[t]]
    }
    scaled <- march.dcmm.h.scaleAlpha(d@M^d@orderHC,alpha[t,])
    alpha[t,]<- scaled$scaled
    l[t] <- log(scaled$scale)
  }

  LL <- log(sum(matrix(1,ncol=1,nrow=d@M^d@orderHC)*alpha[s@N,]))+sum(l*matrix(1,nrow=length(s),ncol=1))
  list(alpha=alpha,l=l,LL=LL)
}



march.dcmm.backward <- function(d,s){
  beta <- matrix(0,nrow=s@N,ncol=d@M^d@orderHC)
  l <- vector("numeric",s@N)

  beta[s@N,] <- matrix(1,ncol=d@M^d@orderHC,nrow=1)


  for( t in march.h.seq(s@N-1,d@orderVC+d@orderHC,-1)){
    rowC <- march.dcmm.h.encodeOutput(s,t+1,d)
    for( i in 1:(d@M^(d@orderHC))){
      for( j in 1:d@M ){
        col <- d@M^(d@orderHC-1)*(j-1)+floor((i-1)/d@M)+1
        beta[t,i] <- beta[t,i]+d@A[i,col]*d@RB[j, rowC ,s@y[t+1]]*beta[t+1,col]
      }
    }
    scaled <- march.dcmm.h.scaleAlpha(d@M^d@orderHC,beta[t,])
    beta[t,] <- scaled$scaled
    l[t] <- log(scaled$scale)
  }

  if( d@orderHC>1 ){
    for( t in march.h.seq(d@orderVC+d@orderHC-1,d@orderVC+1,-1)){
      rowC <- march.dcmm.h.encodeOutput(s,t+1,d)
      for( i in 1:d@M^(t-d@orderVC)){
        for( j in 1:d@M ){
          col <- d@M^(t-d@orderVC)*(j-1)+i
          beta[t,i] <- beta[t,i]+d@Pi[t-d@orderVC+1,i,j]*d@RB[j, rowC ,s@y[t+1]]*beta[t+1,col]
        }
      }
      scaled <- march.dcmm.h.scaleAlpha(d@M^(t-d@orderVC),beta[t,])
      beta[t,] <- scaled$scaled
      l[t] <- log(scaled$scale)
    }
  }

  list(beta=beta,l=l)
}


march.dcmm.epsilon <- function(d,s,alpha,beta,SAlog,SBlog,LLAlpha){
  epsilon <- array(0,c(s@N-1,d@M^d@orderHC,d@M))
  zepsilon <- array(1,c(s@N-1,d@M^d@orderHC,d@M))

  for( t in march.h.seq(d@orderVC+1,d@orderVC+d@orderHC-1)){
    rowC <-  march.dcmm.h.encodeOutput(s,t+1,d)
    for( i in 1:d@M^(t-d@orderVC)){
      j <- 1:d@M
      col <- d@M^(t-d@orderVC)*(j-1)+i
      zepsilon[t,i,j] <- alpha[t,i] & beta[t+1,col] & d@Pi[t-d@orderVC+1,i,j] & d@RB[j, rowC,s@y[t+1]]
    }
  }

  for( t in march.h.seq(d@orderHC+d@orderVC,s@N-1)){
    rowC <- march.dcmm.h.encodeOutput(s,t+1,d)
    for( i in 1:d@M^d@orderHC ){
      j<-1:d@M
      col <- d@M^(d@orderHC-1)*(j-1)+floor((i-1)/d@M)+1
      zepsilon[t,i,j] <- alpha[t,i] & beta[t+1,col] & d@A[i,col] & d@RB[j, rowC,s@y[t+1]]
    }
  }

  LSAfact <- 0
  for( t in march.h.seq(d@orderVC+1,d@orderVC+d@orderHC-1)){
    rowC <- march.dcmm.h.encodeOutput(s,t+1,d)
    LSAfact = LSAfact+SAlog[t]
    for( i in 1:d@M^(t-d@orderVC)){
      for( j in 1:d@M ){
        if( zepsilon[t,i,j] ){
          epsilon[t,i,j] <- LSAfact+log(alpha[t,i])+log(d@Pi[t-d@orderVC+1,i,j])+log(d@RB[j,rowC,s@y[t+1]])
        }
      }
    }
  }

  for( t in march.h.seq(d@orderHC+d@orderVC,s@N-1)){
    rowC <- march.dcmm.h.encodeOutput(s,t+1,d)
    LSAfact = LSAfact+SAlog[t]

    for( i in 1:(d@M^d@orderHC) ){
      for( j in 1:d@M ){
        col <- d@M^(d@orderHC-1)*(j-1)+floor((i-1)/d@M)+1
        if( zepsilon[t,i,j] ){
          epsilon[t,i,j] <- LSAfact+log(alpha[t,i])+log(d@A[i,col])+log(d@RB[j,rowC,s@y[t+1]])
        }
      }
    }
  }

  LSBfact <- 0
  for( t in march.h.seq(s@N-1,d@orderHC+d@orderVC,-1)){
    LSBfact <- LSBfact+SBlog[t+1]
    for( i in 1:d@M^d@orderHC ){
      for( j in 1:d@M ){
        if( zepsilon[t,i,j] ){
          col <- d@M^(d@orderHC-1)*(j-1)+floor((i-1)/d@M)+1
          epsilon[t,i,j] <- exp(epsilon[t,i,j] + LSBfact+log(beta[t+1,col])-LLAlpha)
        }
      }
    }
  }

  if( d@orderHC> 1){
    for( t in march.h.seq(d@orderHC+d@orderVC-1,d@orderVC+1,-1)){
      LSBfact <- LSBfact+SBlog[t+1]
      for( i in 1:d@M^(t-d@orderVC) ){
        for( j in 1:d@M ){
          if( zepsilon[t,i,j] ){
            col <- d@M^(t-d@orderVC)*(j-1)+i
            epsilon[t,i,j] <- exp(epsilon[t,i,j]+LSBfact+log(beta[t+1,col])-LLAlpha)
          }
        }
      }
    }
  }

  epsilon
}

march.dcmm.gamma <- function(d,s,alpha,beta,SAlog,SBlog,LL,epsilon){
  gamma <- array(0,c(s@N,d@M^d@orderHC))

  for( t in march.h.seq(d@orderVC+1,d@orderVC+d@orderHC-1)){
    for( i in 1:d@M^(t-d@orderVC)){
      gamma[t,i] <- sum(epsilon[t,i,])
    }
  }

  for( t in march.h.seq(d@orderVC+d@orderHC,s@N-1)){
    for( i in 1:d@M^d@orderHC){
      gamma[t,i] <- sum(epsilon[t,i,])
    }
  }

  t <- s@N
  for( i in 1:d@M^d@orderHC){
    gamma[t,i] <- exp(log(alpha[t,i])+sum(SAlog)+SBlog[t]-LL)
  }
  gamma
}

#
# Baum-Welch algorithm applied to Dcmm
#
march.dcmm.bw <- function(d,y){

  RA <- array(0,c(d@M^d@orderHC,d@M))
  RB <- array(0,c(d@M,d@y@K^d@orderVC,d@y@K))
  Pi <- array(0,c(d@orderHC,d@M^(d@orderHC-1),d@M))

  for( n in 1:y@N ){
    # extract the current sequence
    s <- march.dataset.h.extractSequence(y,n)

    # First compute all the needed data, depending on d and y
    a <- march.dcmm.forward(d,s)
    b <- march.dcmm.backward(d,s)
    epsilon <- march.dcmm.epsilon(d,s,a$alpha,b$beta,a$l,b$l,a$LL)
    gamma <- march.dcmm.gamma(d,s,a$alpha,b$beta,a$l,b$l,a$LL,epsilon)

    # Reestimation of the RA
    #print(RA)
    for( t in march.h.seq(d@orderVC+d@orderHC,s@N-1)){
      RA <- RA+epsilon[t,,]
    }

    tot <- rowSums(RA)
    #print(d)
    for( i in 1:d@M^d@orderHC ){
      if( tot[i]>0 ){
        RA[i,] <- RA[i,]/tot[i]
      }
    }


    # Reestimation of RB
    for( t in march.h.seq((d@orderVC+1),(d@orderHC+d@orderVC-1))){
      rowC <- march.dcmm.h.encodeOutput(s,t,d)
      for( m in 1:d@M){
        RB[m,rowC,s@y[t]] <- RB[m,rowC,s@y[t]]+sum(gamma[t,((m-1)*d@M^(t-d@orderVC)+1):(m*d@M^(t-d@orderVC))])
      }
    }

    for( t in (d@orderHC+d@orderVC):s@N){
      rowC <- march.dcmm.h.encodeOutput(s,t,d)
      for( m in 1:d@M){
        RB[m,rowC,s@y[t]] <- RB[m,rowC,s@y[t]]+sum(gamma[t,((m-1)*d@M^(d@orderHC-1)+1):(m*d@M^(d@orderHC-1))])
      }
    }

    for( m in 1:d@M ){
      tot <- RB[m,,]%*%march.h.ones(d@y@K,1)

      for( j in 1:(d@y@K^d@orderVC )){
        if( tot[j] ){
          RB[m,j,] <- RB[m,j,]/tot[j]
        }
      }
    }

    # Reestimation of Pi
    Pi[1,1,] <- Pi[1,1,]+gamma[d@orderVC+1,1:d@M]

    for( t in march.h.seq(2,d@orderHC) ){
      for( i in 1:(d@M^(t-1))){
        if( gamma[d@orderVC+t-1,i]>0 )
          Pi[t,i,] <- Pi[t,i,]+epsilon[d@orderVC+t-1,i,]/gamma[d@orderVC+t-1,i]
      }
    }
  }
  Pi <- Pi/y@N

  # Expansion of RA into A
  A <- march.dcmm.h.expandRA(d,RA)


  # finally return a new Dcmm with the reestimated probability distributions
  new("march.Dcmm",Pi=Pi,A=A,M=d@M,y=d@y,orderVC=d@orderVC,orderHC=d@orderHC,RB=RB)
}

#' Construct a double chain Markov model (DCMM).
#'
#' Construct a \code{\link[=march.Dcmm-class]{march.Dcmm}} object, with visible order \emph{orderVC}, hidden order \emph{orderHC} and \emph{M} hidden states, according to a \code{\link[=march.Dataset-class]{march.Dataset}}.
#' The first \emph{maxOrder}-\emph{orderVC} elements of each sequence are truncated in order to return a model
#' which can be compared with other Markovian model of visible order maxOrder. The construction is performed either by an evolutionary algorithm (EA) or by improving an existing DCMM.
#' The EA performs \emph{gen} generations on a population of \emph{popSize} individuals. The EA behaves as a Lamarckian evolutionary algorithm, using a Baum-Welch algorithm as
#' optimization step, running until log-likelihood improvement is less than \emph{stopBw} or for \emph{iterBw} iterations. Finally only the best individual from the population is returned as solution.
#' If a seedModel is provided, the only step executed is the optimization step, parameters related to the EA does not apply in this case.
#'
#' @param y the dataset from which the Dcmm will be constructed \code{\link[=march.Dataset-class]{march.Dataset}}.
#' @param orderHC the order of the hidden chain of the constructed Dcmm.
#' @param orderVC the order of the visible chain of the constructed Dcmm (0 for a HMM).
#' @param M the number of hidden state of the Dcmm.
#' @param gen the number of generation performed by the EA.
#' @param popSize the number of individual stored into the population.
#' @param maxOrder the maximum visible order among the set of Markovian models to compare.
#' @param seedModel a model to optimize using Baum-Welch algorithm.
#' @param iterBw the number of iteration performed by the Baum-Welch algorithm.
#' @param stopBw the minimum increase in quality (log-likelihood) authorized in the Baum-Welch algorithm.
#'
#' @return the best \code{\link[=march.Dcmm-class]{march.Dcmm}} constructed by the EA or the result of the Baum-Welch algorithm on \emph{seedModel}.
#'
#' @author Ogier Maitre
#' @example tests/examples/march.dcmm.construct.example.R
#' @seealso \code{\link[=march.Dcmm-class]{march.Dcmm-class}}, \code{\link{march.Model-class}}, \code{\link[=march.Dataset-class]{march.Dataset-class}}.
#'
#' @export
march.dcmm.construct <- function(y,orderHC,orderVC,M,gen=5,popSize=4,maxOrder=orderVC,seedModel=NULL,iterBw=2,stopBw=0.1){

  if( is.null(seedModel) ){
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

  #
  if( is.null(seedModel) ){
    y <- march.dataset.h.filtrateShortSeq(y,maxOrder+1)
    y <- march.dataset.h.cut(y,maxOrder-orderVC)
  }

  if( is.null(seedModel)==FALSE ){

    # AB
    y <- march.dataset.h.filtrateShortSeq(y,maxOrder+1)
    y <- march.dataset.h.cut(y,maxOrder-orderVC)
    # \AB

    op <- new("march.dcmm.ea.OptimizingParameters",fct=march.dcmm.ea.optimizing,ds=y,stopBw=stopBw,iterBw=iterBw)
    m <- march.dcmm.ea.optimizing(seedModel,op)

    # AB
    m@dsL <- sum(y@T-orderVC)
    #m@dsL <- sum(y@T)
    # \AB

    m@y <- y
    #m@ll <- march.dcmm.forward(m,march.dataset.h.extractSequence(y,1))$LL
    m@ll <- march.dcmm.h.computeLL(m,y)
    m
  }
  else{
    ep <- new("march.dcmm.ea.EvalParameters",ds=y,fct=march.dcmm.ea.evaluation)
    ip <- new("march.dcmm.ea.InitParameters",AConst=FALSE,y=y,orderVC=orderVC,orderHC=orderHC,M=M,fct=march.dcmm.ea.initialization)
    mp <- new("march.dcmm.ea.MutationParameters",pMut=as.numeric(0.05),fct=march.dcmm.ea.mutation)
    cp <- new("march.ea.CrossoverParameters",fct=march.dcmm.ea.crossover)

    op <- new("march.dcmm.ea.OptimizingParameters",fct=march.dcmm.ea.optimizing,ds=y,stopBw=stopBw,iterBw=iterBw)
    p <- new("march.ea.Parameters",optimizing=TRUE,
             initParameters=ip,evalParameters=ep,mutationParameters = mp, optimizingParameters=op,crossoverParameters=cp,
             populationSize=popSize,crossoverProb=0.5,generation=gen)
    res<-march.loop(p)

    # AB
    res$best@dsL <- sum(y@T-orderVC)
    #res$best@dsL <- sum(y@T)
    # \AB
    res$best@y <- y
    res$best@ll <-march.dcmm.h.computeLL(res$best,y)
    res$best
  }
}
