

march.mtd.h.constructEmptyMtd <- function(order,k){
  Q <- array(0,c(k^order,k^order))
  phi <- array(0,c(order))

  new("march.Mtd",Q=Q,phi=phi,K=k,order=order)
}

# Construct nt vector, holding the number of data items per sequence (lenght of data series)
BuildArrayNumberOfDataItems <- function(x){
  n_rows_data <- dim(x)[1] # number of rows (number of data sequences)
  n_cols_data <- dim(x)[2] # number of rows (number of data sequences)
  nt <- array(0,c(1,n_rows_data)) # holds the number of data items per sequence
  for (i in 1:n_rows_data){
    nt[i] <- match(0,x[i,]) # first match of 0
    if(is.na(nt[i])) nt[i] <- n_cols_data # otherwise, just the number of columns of the matrix (since if there is no 0 all positions are valid data)
  }
  return(nt)
}

#% Computation of the contingency table (Cg in page 385 of Berchtold, 2001)
BuildContingencyTable <- function(y,order){
  n_rows_data <- y@N # number of rows (number of data sequences)
  c <- array(data=0,dim=c(order,y@K,y@K)) # crosstable (rt, Cg in page 385 of Berchtold, 2001)
  for (g in 1:order){
    for (i in 1:n_rows_data){
      for (t in 1:(y@T[i]-g)){ # mc_lag is g in Berchtold, 2001
        past <- y@y[[i]][t]
        present <- y@y[[i]][t+g]
        if(length((which(c(past,present)<1) | (which(c(past,present)>y@K))))==0){
          c[g,past,present] <- c[g,past,present] + y@weights[i]
        }
      }
    }
  }
  return(c)
}

NormalizeTable <- function(x){
  nx <- array(NA,dim=dim(x))
  m <- dim(x)[1]
  sums <- rowSums(x)
  for (i in 1:m){
    if (sums[i]!=0) nx[i,] <- x[i,]/sums[i]
  }
  return(nx)
}

#% Computation of the measure u for a crosstable
# See p. 387 of Berchtold, 2001
CalculateTheilU <- function(m,l,c){

  u <- array(data=NA,dim=c(1,l))

  for (g in 1:l){
    cg <- c[g,,]
    tc <- sum(cg) # sum of elements (TCg)
    sr <- rowSums(cg) # vector of sums of rows [Cg(.,j)]
    sc <- colSums(cg) # vector of sums of columns [Cg(i,.)]

    # the following lines implement equation 14
    num <- 0
    for (i in 1:m){
      for (j in 1:m){
        if (cg[i,j]!=0){ # if cg[i,j], there is nothing to be added
          num <- num + cg[i,j] * (log2(sc[i]) + log2(sr[j]) - log2(cg[i,j]) - log2(tc))
        }
      }
    }

    den <- 0

    for (j in 1:m){
      if(sc[j]!=0 & tc!=0){
        den <- den + sr[j] * (log2(sr[j]) - log2(tc))
      }
    }

    if(den != 0) u[g] = num/den
    else u[g] = NaN

  }
  return(u)
}

InitializeParameters <- function(u,init_method,c,is_mtdg,m,order){
  # Initialization of the lag parameters (Eq. 15 of Berchtold, 2001)
  phi <- u/(sum(u)) # TODO : Supprimer la division pour coller a march
  if (is_mtdg){
    q <- array(NA,c(order,m,m))
    for (g in 1:order){
      q[g,,] <- NormalizeTable(c[g,,])
    }
    return(list(phi=phi,q=q))
  } else {
    # Initialization of the initial matrix Q (see page 387 of Berchtold, 2001)
    if (init_method == "weighted"){
      # Initialization of Q as a weighted sum of matrices Q1,...,Ql
      q_tilde <- array(0,dim=c(1,m,m))
      for (g in 1:order){
        q_tilde[1,,] <- q_tilde[1,,] + phi[g] * c[g,,]
      }
      q <- NormalizeTable(q_tilde)
    } else if (init_method == "best"){
      # Initialization of Q as Q=Qk where k = argmax(ug)
      k <- which.max(u)
      q <- array(0,c(1,m,m))
      q[1,,] <- NormalizeTable(c[k,,])
    } else if (init_method == "random"){
      # Initialization of Q as a random matrix
      q <- 0.1 + array(runif(m*m),c(1,m,m))
      q[1,,] <- q/rowSums(q[1,,])
    } else{
      stop("Init parameter should be either, \"best\", \"random\" or \"weighted\"",call.=FALSE)
    }
    return(list(phi=phi,q=q))
  }
}

# Construct the array i0_il with all possible combinations of states 1...m in a time window of size l+1
BuildArrayCombinations <- function(m,l){
  i0_il <- array(0,c(m^(l+1),l+1))
  values <- 1:m
  for(i in 1:(l+1)){
    i0_il[,i] <- t(kronecker(values,rep(1,m^(l+1)/m^i)))
    values <- kronecker(rep(1,m),values)
  }
  return(i0_il)
}

# Construct array n_i0_il where n_i0_il[i0,...il] is the number of sequences of the form X(t-l)=il,...,X(t)=i0
BuildArrayNumberOfSequences <- function(y,order){
  n_i0_il <- array(0,dim=rep(y@K,order+1))
  for(i in 1:y@N){
    for(t in march.h.seq(1,y@T[i]-order)){
      ind <- y@y[[i]][t:(t+order)]
      n_i0_il[rbind(ind)] <- n_i0_il[rbind(ind)] + 1
    }

  }
  # Transform to a one-dimensional array
  n_i0_il <- c(n_i0_il)
  return(n_i0_il)
}

# Transform transition matrix q in such a way that the log-likelihood (Eq. 6) and the partial derivatives (equations in page 382) are easy to calculate
BuildArrayQ <- function(m,l,i0_il,n_i0_il,q){
  q_i0_il <- array(0,c(m^(l+1),l))
  for (j in 1:length(n_i0_il)){
    if( dim(q)[1]>1){
      for (g in 1:l){
        q_i0_il[j,g] <- q[g,i0_il[j,g+1],i0_il[j,1]]
      }
    }
    else {
      for (k in 1:l){
        q_i0_il[j,k] <- q[1,i0_il[j,k+1],i0_il[j,1]]
      }
    }
  }
  return(q_i0_il)
}

CalculateLogLikelihood <- function(n_i0_il,q_i0_il,phi){
  ll <- 0
  for (i in 1:length(n_i0_il)){
    if(n_i0_il[i] > 0){
      ll <- ll + n_i0_il[i] * log(q_i0_il[i,] %*% t(phi) )
    }
  }
  ll
}

# Calculate the partial derivatives of log(L) with respect to phi_k (page 382)
PartialDerivativesPhi <- function(n_i0_il,q_i0_il,l,phi){
  pd_phi <- rep(0,l)
  for (i in 1:length(n_i0_il)){
    for (g in 1:l){
      if (n_i0_il[i] > 0 && (q_i0_il[i,] %*% t(phi)) != 0){
        pd_phi[g] <- pd_phi[g] + n_i0_il[i] * q_i0_il[i,g] / (q_i0_il[i,] %*% t(phi))
      }
    }
  }
  return(pd_phi)
}

# Calculate the partial derivatives of log(L) with respect to q_ik_i0 (page 382)
PartialDerivativesQ <- function(n_i0_il,i0_il,q_i0_il,m,phi,order){
  pd_q <- array(0,dim=c(m,m))
  for (i in 1:length(n_i0_il)){
    for (j in 1:order){
      if (n_i0_il[i] > 0 && (q_i0_il[i,] %*% t(phi)) != 0){
        pd_q[i0_il[i,j+1],i0_il[i,1]] <- pd_q[i0_il[i,j+1],i0_il[i,1]] + n_i0_il[i] * phi[j] / ( q_i0_il[i,] %*% t(phi) )
      }
    }
  }
  return(pd_q)
}

OptimizePhi <- function(phi,pd_phi,delta,is_constrained,delta_stop,ll,n_i0_il,q_i0_il,k){

  delta_it <- delta[1]

  i_inc <- which.max(pd_phi) # index of the phi parameter to increase (the one corresponding to the largest derivative)
  i_dec <- which.min(pd_phi) # index of the phi parameter to decrease (the one corresponding to the smallest derivative)
  par_inc <- phi[i_inc]
  par_dec <- phi[i_dec]

  if (is_constrained){

    if (par_inc == 1){ # "If phi_plus = 1, this parameter cannot be increased and the algorithm cannot improve the log-likelihood through a reestimation of the vector phi" (page 382)
      # "If the distribution was not reevaluated because the parameter to increase was already set to 1, delta is not changed." (page 383)
      return(list(phi=phi,ll=ll,delta=delta))
    }
    #if ( par_inc + delta_it > 1){ # "If phi_plus + delta > 1, the quantity delta is too large and it must be set to delta = 1 - phi_plus." (page 382)
    #    delta_it <- 1 - par_inc
    #}
    if ( par_dec == 0 ){ # "If phi_minus = 0, this parameter cannot be decreased. Then, we decrease the parameter corresponding to the smallest strictly positive derivative"
      pd_phi_sorted <- sort(pd_phi,index.return=TRUE)
      i_dec <- pd_phi_sorted$ix[min(which(phi[pd_phi_sorted$ix]>0))]
      par_dec <- phi[i_dec]
    }
    #if ( par_dec - delta_it < 0){ # "If phi_minus - delta < 0, the quantity delta is too large and it must be set to delta = 1 - phi_minus."
    #    delta_it <- 1 - par_dec
    #}
    delta_it <- min(c(delta_it,1-par_inc,par_dec))
  }

  new_phi <- phi

  while(TRUE){
    new_phi[i_inc] <- par_inc + delta_it
    new_phi[i_dec] <- par_dec - delta_it
    if (!is_constrained){ # constraints Eq. 4 Berchtold, 2001 p. 380
      t <- sum(new_phi[new_phi >= 0])
      for (i in 1:k){
        q_min <- min(q[i,])
        q_max <- max(q[i,])
        if (t*q_min + (1-t)*q_max < 0){
          return(list(phi=new_phi,ll=new_ll,delta=delta))
        }
      }
    }
    new_ll <- CalculateLogLikelihood(n_i0_il,q_i0_il,new_phi) # "Once the new vector phi is known, we can compute the new log-likelihood of the model." (page 382)
    if (new_ll > ll){ # "If it is larger than the previous value, the new vector phi is accepted and the procedure stops." (page 382)
      if (delta_it == delta[1]){ # "If the distribution was reevaluated with the original value of delta [...], delta is set to 2*delta" (page 383)
        delta[1] <- 2*delta[1]
      } # "if the distribution was reevaluated with a value of delta smaller than its value at the beginning of Step 2, delta keeps its present value"
      return(list(phi=new_phi,ll=new_ll,delta=delta))
    } else { # In the other case, delta is divided by 2 and the procedure iterates
      if (delta_it <= delta_stop) { # "When delta becomes smaller than a fixed threshold, we stop the procedure, even if phi was not reevaluated"
        delta[1] <- 2*delta[1] # "If the algorithm was completed without reestimation, delta is set to twice the value reached at the end of Step 2"
        return(list(phi=phi,ll=ll,delta=delta))
      }
      delta_it <- delta_it/2
    }
  }

}

OptimizeQ <- function(q,j,pd_q,delta,delta_stop,ll,n_i0_il,q_i0_il,phi,i0_il,k,l,g){
  delta_it <- delta[j+1]

  i_inc <- which.max(pd_q[j,]) # index of the phi parameter to increase (the one corresponding to the largest derivative)
  i_dec <- which.min(pd_q[j,]) # index of the phi parameter to decrease (the one corresponding to the smallest derivative)
  par_inc <- q[g,j,i_inc]
  par_dec <- q[g,j,i_dec]

  if (par_inc == 1){ # "If phi_plus = 1, this parameter cannot be increased and the algorithm cannot improve the log-likelihood through a reestimation of the vector phi" (page 382)
    # "If the distribution was not reevaluated because the parameter to increase was already set to 1, delta is not changed." (page 383)
    return(list(q=q,q_i0_il=q_i0_il,ll=ll,delta=delta))
  }
  #if ( par_inc + delta_it > 1){ # "If phi_plus + delta > 1, the quantity delta is too large and it must be set to delta = 1 - phi_plus." (page 382)
  #   delta_it <- 1 - par_inc
  #}
  if ( par_dec == 0 ){ # "If phi_minus = 0, this parameter cannot be decreased. Then, we decrease the parameter corresponding to the smallest strictly positive derivative"
    pd_q_sorted <- sort(pd_q[j,],index.return=TRUE)
    i_dec <- pd_q_sorted$ix[min(which(q[g,j,pd_q_sorted$ix]>0))]
    par_dec <- q[g,j,i_dec]
  }
  #if ( par_dec - delta_it < 0){ # "If phi_minus - delta < 0, the quantity delta is too large and it must be set to delta = 1 - phi_minus."
  #    delta_it <- 1 - par_dec
  #}

  delta_it <- min(c(delta_it,1-par_inc,par_dec))

  new_q_row <- q[g,j,]

  while(TRUE){
    new_q_row[i_inc] <- par_inc + delta_it
    new_q_row[i_dec] <- par_dec - delta_it
    new_q <- q
    new_q[g,j,] <- new_q_row
    new_q_i0_il <- BuildArrayQ(k,l,i0_il,n_i0_il,new_q)
    new_ll <- CalculateLogLikelihood(n_i0_il,new_q_i0_il,phi) # "Once the new vector phi is known, we can compute the new log-likelihood of the model." (page 382)
    if (new_ll > ll){ # "If it is larger than the previous value, the new vector phi is accepted and the procedure stops." (page 382)
      if (delta_it == delta[j+1]){ # "If the distribution was reevaluated with the original value of delta [...], delta is set to 2*delta" (page 383)
        delta[j+1] <- 2*delta[j+1]
      } # "if the distribution was reevaluated with a value of delta smaller than its value at the beginning of Step 2, delta keeps its present value"
      return(list(q=new_q,q_i0_il=new_q_i0_il,ll=new_ll,delta=delta))
    } else { # In the other case, delta is divided by 2 and the procedure iterates
      if (delta_it <= delta_stop) { # "When delta becomes smaller than a fixed threshold, we stop the procedure, even if phi was not reevaluated"
        delta[j+1] <- 2*delta[j+1] # "If the algorithm was completed without reestimation, delta is set to twice the value reached at the end of Step 2"
        return(list(q=q,q_i0_il=q_i0_il,ll=ll,delta=delta))
      }
      delta_it <- delta_it/2
    }
  }

}

#' Construct a Mixture Transition Distribution (MTD) model.
#'
#' A Mixture Transition Distribution model (\code{\link[=march.Mtd-class]{march.Mtd}}) object of order \emph{order} is constructed
#' according to a given \code{\link[=march.Dataset-class]{march.Dataset}} \emph{y}. The first \emph{maxOrder}-\emph{order}
#' elements of each sequence are truncated in order to return a model
#' which can be compared with other Markovian models of visible order maxOrder.
#'
#' @param y the dataset (\code{\link[=march.Dataset-class]{march.Dataset}}) from which to construct the model.
#' @param order the order of the constructed model.
#' @param maxOrder the maximum visible order among the set of Markovian models to compare.
#' @param mtdg flag indicating whether the constructed model should be a MTDg using a different transition matrix for each lag (value: \emph{TRUE} or \emph{FALSE}).
#' @param init the init method, to choose among \emph{best}, \emph{random} and \emph{weighted}.
#' @param deltaStop the delta below which the optimization phases of phi and Q stop.
#' @param llStop the ll increase below which the EM algorithm stop.
#' @param maxIter the maximal number of iterations of the optimisation algorithm (zero for no maximal number).
#'
#' @author Ogier Maitre
#' @example tests/examples/march.mtd.construct.example.R
#' @seealso \code{\link[=march.Mtd-class]{march.Mtd-class}}, \code{\link{march.Model-class}}, \code{\link[=march.Dataset-class]{march.Dataset-class}}.
#' @export
march.mtd.construct <- function(y,order,maxOrder=order,mtdg=FALSE,init="best", deltaStop=0.0001, llStop=0.01, maxIter=0){

  order <- march.h.paramAsInteger(order)
  maxOrder <- march.h.paramAsInteger(maxOrder)

  if( order>maxOrder ){
    stop("maxOrder should be greater or equal than order")
  }

  ySave <- y

  # here we create a new data set by truncating the first (maxOrder-order+1) element
  # this data set will be used to construct the mtd model, instead of the original one.
#   for( i in 1:y@N ){
#     y@y[[i]] <- y@y[[i]][(maxOrder-order+1):length(y@y[[i]])]
#     y@T[i] <- length(y@y[[i]])
#   }

  y <- march.dataset.h.filtrateShortSeq(y,maxOrder+1)
  y <- march.dataset.h.cut(y,maxOrder-order)

  is_constrained <- TRUE # if FALSE the model is unconstrained and the constraints given by Eq. 4 are activated instead of those given by Eq. 3 (Berchtold, 2001, p. 380)
  is_mtdg <- mtdg
  init_method <- init

  # 1.1. Choose initial values for all parameters
  # nt <- BuildArrayNumberOfDataItems(y) This is no longer needed, as y now contains the T field
  c <- BuildContingencyTable(y,order)
  u <- CalculateTheilU(y@K,order,c)
  init_params <- InitializeParameters(u=u,init_method=init_method,c=c,is_mtdg=is_mtdg,m=y@K,order=order)
  phi <- init_params$phi
  q <- init_params$q


  # 1.2. Choose a value for delta and a criterion to stop the algorithm
  delta <- array(0.1,dim=c(1,y@K+1)) # a different delta is used for each of the m+1 sets of parameters (vector phi + each row of the transition matrix Q)
  delta_stop <- deltaStop
  ll_stop <- llStop

  # 2. Iterations
  # 2.1. Reestimate the vector phi by modifying two of its elements
  # 2.2. Reestimate the transition matrix Q by modifying two elements of each row

  # 3. End criterion
  # 3.1. If the increase of the LL since the last iteration is greater than the stop criterion, go back to step 2
  # 3.2. Otherwise, end the procedure

  i0_il <- BuildArrayCombinations(y@K,order)
  n_i0_il <- BuildArrayNumberOfSequences(y,order)
  q_i0_il <- BuildArrayQ(m=y@K,l=order,i0_il=i0_il,n_i0_il = n_i0_il,q=q)
  new_ll <- CalculateLogLikelihood(n_i0_il=n_i0_il,q_i0_il=q_i0_il,phi=phi)

  iter <- 0
  while (TRUE){
    if( maxIter>0 && iter>=maxIter ){ break }
    else{ iter <- iter+1 }


    ll <- new_ll
    # "Reestimate the vector phi by modifying two of its elements." (p. 383)
    pd_phi <- PartialDerivativesPhi(n_i0_il=n_i0_il,q_i0_il=q_i0_il,l=order,phi = phi)
    res_opt_phi <- OptimizePhi(phi=phi,pd_phi=pd_phi,delta=delta,is_constrained=is_constrained,delta_stop=delta_stop,ll=ll,n_i0_il=n_i0_il,q_i0_il=q_i0_il)
    phi <- res_opt_phi$phi
    new_ll <- res_opt_phi$ll
    delta <- res_opt_phi$delta


    # "Reestimate the transition matrix Q by modifying two elements of each row." (p. 383)
    for (j in 1:y@K){
      qTmp <- array(NA,c(dim(q)[1],y@K,y@K))
      llTmp <- array(NA,c(dim(q)[1]))
      pd_q <- PartialDerivativesQ(n_i0_il,i0_il,q_i0_il,y@K,phi=phi,order=order)

      # w.r.t. the partial derivative, find which matrix q_g should be modified
      for( g in 1:dim(q)[1] ){
        res_opt_q <- OptimizeQ(q=q,j=j,pd_q=pd_q,delta=delta,delta_stop=delta_stop,ll=new_ll,n_i0_il=n_i0_il,q_i0_il=q_i0_il,phi=phi,i0_il = i0_il,k=y@K,l=order,g=g)
        qTmp[g,,] <- res_opt_q$q[g,,]
        q_i0_il <- res_opt_q$q_i0_il
        llTmp[g] <- res_opt_q$ll
        delta <- res_opt_q$delta
      }
      # use the one maximizing the ll
      llMaxId <- which.max(llTmp)
      q[llMaxId,,] <- qTmp[llMaxId,,]
    }
    new_ll <- llTmp[llMaxId]
    if (new_ll - ll < ll_stop){ break }
  }

  nbZeros <- length(which(q==0))+length(which(phi==0))

  # construct and return the final object.
  new("march.Mtd",order=order,Q=q,phi=phi,ll=ll,y=ySave,dsL=sum(y@T-order),nbZeros=nbZeros)
}
