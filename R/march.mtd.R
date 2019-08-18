
# march.mtd.h.constructEmptyMtd <- function(order,k){
#   Q <- array(0,c(k^order,k^order))
#   phi <- array(0,c(order))
# 
#   new("march.Mtd",Q=Q,phi=phi,K=k,order=order)
# }

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

#Build the contingency tables. The first order table concerns the contingency between the lag g and the present. The next tables
#concern the contingency table between the covariate and the 
 BuildContingencyTable <- function(y,order){
   n_rows_data <- y@N # number of rows (number of data sequences)
   l<-list() # crosstable (rt, Cg in page 385 of Berchtold, 2001)
   for (g in 1:order){
     c=array(0,c(y@K,y@K))
     for (i in 1:n_rows_data){
       for (t in 1:(y@T[i]-g)){ # mc_lag is g in Berchtold, 2001
         past <- y@y[[i]][t]
         present <- y@y[[i]][t+g]
         if(length((which(c(past,present)<1) | (which(c(past,present)>y@K))))==0){
           c[past,present] <- c[past,present] + y@weights[i]
         }
       }
     }
     l[[g]]<-c
   }

   if(y@Ncov>0){
     for(j in 1:y@Ncov){
       kcov=y@Kcov[j]
       CT=matrix(0,kcov,y@K)

       for(n in 1:y@N){
         for (i in 1:y@T[n]){
          row=y@cov[n,i,j]
           col=y@y[[n]][i]
           CT[row,col]=CT[row,col]+y@weights[n]
        }
       }
       l[[order+j]]<-CT
     }
   }
   return(l)
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
CalculateTheilU <- function(y,order,c){

  u <- array(data=NA,dim=c(1,order+y@Ncov))

  for (g in 1:order){
    cg <- c[[g]]
    tc <- sum(cg) # sum of elements (TCg)
    sr <- rowSums(cg) # vector of sums of rows [Cg(.,j)]
    sc <- colSums(cg) # vector of sums of columns [Cg(i,.)]

    # the following lines implement equation 14
    num <- 0
    for (i in 1:y@K){
      for (j in 1:y@K){
        if (cg[i,j]!=0){ # if cg[i,j], there is nothing to be added
          num <- num + cg[i,j] * (log2(sc[i]) + log2(sr[j]) - log2(cg[i,j]) - log2(tc))
        }
      }
    }

    den <- 0

    for (j in 1:y@K){
      if(sc[j]!=0 & tc!=0){
        den <- den + sc[j] * (log2(sc[j]) - log2(tc))
      }
    }

    if(den != 0) u[g] = num/den
    else u[g] = NaN

  }
  if(y@Ncov>0){
  for (g in 1:y@Ncov){
    cg <- c[[order+g]]
    tc <- sum(cg) # sum of elements (TCg)
    sr <- rowSums(cg) # vector of sums of rows [Cg(.,j)]
    sc <- colSums(cg) # vector of sums of columns [Cg(i,.)]
    
    # the following lines implement equation 14
    num <- 0
    for (i in 1:y@Kcov[g]){
      for (j in 1:y@K){
        if (cg[i,j]!=0){ # if cg[i,j], there is nothing to be added
          num <- num + cg[i,j] * (log2(sc[j]) + log2(sr[i]) - log2(cg[i,j]) - log2(tc))
        }
      }
    }
    
    den <- 0
    
    for (j in 1:y@K){
      if(sc[j]!=0 & tc!=0){
        den <- den + sc[j] * (log2(sc[j]) - log2(tc))
      }
    }
    
    if(den != 0) u[order+g] = num/den else u[order+g] = NaN
    
  }
}
  return(u)
}

InitializeParameters <- function(u,init_method,c,is_mtdg,m,order,kcov,ncov){
  # Initialization of the lag parameters (Eq. 15 of Berchtold, 2001)
  phi <- u/(sum(u)) # TODO : Supprimer la division pour coller a march
  if (is_mtdg){
    q <- array(NA,c(order,m,m))
    for (g in 1:order){
      q[g,,] <- NormalizeTable(c[[g]])
    }
  } else {
    # Initialization of the initial matrix Q (see page 387 of Berchtold, 2001)
    if (init_method == "weighted"){
      # Initialization of Q as a weighted sum of matrices Q1,...,Ql
      q <- array(0,c(1,m,m))
      q_tilde <- array(0,dim=c(1,m,m))
      for (g in 1:order){
        q_tilde[1,,] <- q_tilde[1,,] + u[g] * c[[g]]
      }
      q[1,,] <- NormalizeTable(q_tilde[1,,])
    } else if (init_method == "best"){
      # Initialization of Q as Q=Qk where k = argmax(ug)
      u_tilde<-u[1:order]
      k <- which.max(u_tilde)
      q <- array(0,c(1,m,m))
      q[1,,] <- NormalizeTable(c[[k]])
    } else if (init_method == "random"){
      # Initialization of Q as a random matrix
      q <- 0.1 + array(runif(m*m),c(1,m,m))
      q[1,,] <- q/rowSums(q[1,,])
    } else{
      stop("Init parameter should be either, \"best\", \"random\" or \"weighted\"",call.=FALSE)
    }
  }
  S=list()
  #Initialisation of the matrices of transition between covariates and the dependant variable by normalizing the crosstables.
  if(ncov>0){
    for(i in 1:ncov){
      S[[i]]=matrix(0,kcov[i],m)
      S[[i]]<-NormalizeTable(c[[order+i]])
    }
  }
    
    return(list(phi=phi,q=q,S=S))
  
}

# Construct the array i0_il with all possible combinations of states 1...m together with the covariates in a time window of size l+1
BuildArrayCombinations <- function(m,l,kcov,ncov){
  tCovar <- 1
  
  if(prod(kcov)>0){
  tCovar<-prod(kcov)
  }
  
  i0_il <- array(0,c(m^(l+1)*tCovar,l+1+ncov))
  values <- 1:m
  for(i in 1:(l+1)){
    i0_il[,i] <- t(kronecker(values,rep(1,m^(l+1)*tCovar/m^i)))
    values <- kronecker(rep(1,m),values)
  }
  
  if(ncov>0){
    totm <- tCovar
    totp <- m^(l+1)
    for(i in 1:ncov){
      totm <- totm/kcov[i]
      values <- 1:kcov[i]
      i0_il[,l+1+i] <- kronecker(rep(1,totp),kronecker(values,rep(1,totm)))
      totp <- totp*kcov[i]
    }
  }
  return(i0_il)
}

BuildArrayNumberOfSequences <- function(y,order){
  
  if(y@Ncov>0){
  tCovar<-prod(y@Kcov)
  }else{
    tCovar=1
  }
  
  n_i0_il <- array(0,dim=c(y@K^(order+1)*tCovar,1))
  for(n in 1:y@N){
    for(t in march.h.seq(1,y@T[n]-order)){
      pos <- y@y[[n]][t:(t+order)]
      row1=pos[1]
      for(g in 2:(order+1)){
        row1=row1+(y@K)^(g-1)*(pos[g]-1)
      }
      if(y@Ncov>0){
        posCov <- y@cov[n,t+order,]
        Covar<-tCovar
        row2=posCov[y@Ncov]
        totKC=1
        if(y@Ncov>1){
          for(i in (y@Ncov-1):1 ){
            totKC=totKC*y@Kcov[i+1]
            row2=row2+totKC*(y@cov[n,t+order,i]-1)
          }
        }
      }
      if(y@Ncov>0){
        row1=(row1-1)*tCovar+row2
      }
      n_i0_il[row1]<-n_i0_il[row1]+1
    }

  }
  # Transform to a one-dimensional array
  #n_i0_il <- c(n_i0_il)
  return(n_i0_il)
}

BuildArrayQ <- function(m,l,i0_il,n_i0_il,q,kcov,ncov,S){
  
  if(ncov>0){
    tCovar<-prod(kcov)
  }else{
    tCovar=1
  }
  
  q_i0_il <- array(0,c(m^(l+1)*tCovar,l+ncov))
  for (i in 1:length(n_i0_il)){
    if( dim(q)[1]>1){
      for (j in 1:l){
        q_i0_il[i,j] <- q[j,i0_il[i,j+1],i0_il[i,1]]
      }
    }else {
      for (j in 1:l){
        q_i0_il[i,j] <- q[1,i0_il[i,j+1],i0_il[i,1]]
      }
    }
    if(ncov>0){
      for(j in 1:ncov){
        q_i0_il[i,l+j] <- S[[j]][i0_il[i,l+1+j],i0_il[i,1]]
      }
    }
  }
  return(q_i0_il)
}

CalculateLogLikelihood <- function(n_i0_il,q_i0_il,phi){
  ll <- 0
  for (i in 1:length(n_i0_il)){
    if(n_i0_il[i] > 0){
      ll <- ll + n_i0_il[i] * log(q_i0_il[i,]%*%t(phi) )
    }
  }
  ll
}

# Calculate the partial derivatives of log(L) with respect to phi_k (page 382)
PartialDerivativesPhi <- function(n_i0_il,q_i0_il,l,phi,ncov){
  pd_phi <- rep(0,l+ncov)
  for (k in 1:length(n_i0_il)){
   for (m in 1:(l+ncov)){
      if (n_i0_il[k] > 0 && (q_i0_il[k,] %*% t(phi)) != 0){
        pd_phi[m] <- pd_phi[m] + n_i0_il[k] * q_i0_il[k,m] / (q_i0_il[k,] %*% t(phi))
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

PartialDerivativesS<- function(CCov,ColVT,k,kcov,n_i0_il,q_i0_il,i0_il,phi){
  pd_s<-array(0,dim=c(kcov,k))
  nc=length(n_i0_il)
  
  for (k in 1:nc){
    if (n_i0_il[k] > 0 && (q_i0_il[k,] %*% t(phi)) != 0){
      pd_s[i0_il[k,ColVT],i0_il[k,1]]<-pd_s[i0_il[k,ColVT],i0_il[k,1]]+n_i0_il[k]*phi[CCov]/(q_i0_il[k,] %*% t(phi))
    }
  }
  return(pd_s)
  
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
  }

  

  while(TRUE){
    if(is_constrained){
    delta_it <- min(c(delta_it,1-par_inc,par_dec))
    }
    new_phi <- phi
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

OptimizeQ <- function(q,j,pd_q,delta,delta_stop,ll,n_i0_il,q_i0_il,phi,i0_il,k,l,g,kcov,ncov,S){
  delta_it <- delta[j+1]

  i_inc <- which.max(pd_q[j,]) # index of the Q parameter to increase (the one corresponding to the largest derivative)
  i_dec <- which.min(pd_q[j,]) # index of the Q parameter to decrease (the one corresponding to the smallest derivative)
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


  while(TRUE){
    delta_it <- min(c(delta_it,1-par_inc,par_dec))
    
    new_q_row <- q[g,j,]
    new_q_row[i_inc] <- par_inc + delta_it
    new_q_row[i_dec] <- par_dec - delta_it
    new_q <- q
    new_q[g,j,] <- new_q_row
    new_q_i0_il <- BuildArrayQ(k,l,i0_il,n_i0_il,new_q, kcov,ncov,S)
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

OptimizeS <-function(order,k,kcov,ncov,S,Tr,phi,pcol,ll,pd_s,delta,delta_stop,n_i0_il,i0_il,q_i0_il,q,n){
  
  delta_it<-delta
  i_inc<-which.max(pd_s)
  i_dec<-which.max(pd_s)
  par_inc<-S[[n]][Tr,i_inc]
  par_dec<-S[[n]][Tr,i_dec]
  
  if(par_inc==1){
    return(list(S=S,ll=ll,delta=delta,q_i0_il=q_i0_il))
  }
  if(par_dec==0){
    pd_s_sorted<-sort(pd_s,index.return=TRUE)
    i_dec<-pd_s_sorted$ix[min(which(S[Tr,pd_s_sorted$ix]>0))]
    par_dec<-S[[n]][Tr,i_dec]
  }
  
  delta_it <- min(c(delta_it,1-par_inc,par_dec))
  new_S_row<-S[[n]][Tr,]
 
  while(TRUE){
    new_S_row[i_inc]<-par_inc+delta_it
    new_S_row[i_dec]<-par_dec-delta_it
    new_S<-S
    new_S[[n]][Tr,]<-new_S_row
    new_q_i0_il<-BuildArrayQ(k,order,i0_il,n_i0_il,q,kcov,ncov,new_S)
    new_ll <- CalculateLogLikelihood(n_i0_il,new_q_i0_il,phi)
    if(new_ll>ll){
      if(delta_it==delta){
        delta<-2*delta
      }
      return(list(S=new_S,ll=new_ll,delta=delta,q_i0_il=new_q_i0_il))
    }else{
      if(delta_it<=delta_stop){
        delta<-2*delta
        return(list(S=S,ll=ll,delta=delta,q_i0_il=q_i0_il))
      }
      delta_it<-delta_it/2
    }
  }
}
#' Construct a Mixture Transition Distribution (MTD) model.
#'
#' A Mixture Transition Distribution model (\code{\link{march.Mtd-class}}) object of order \emph{order} is constructed
#' according to a given \code{\link{march.Dataset-class}} \emph{y}. The first \emph{maxOrder}-\emph{order}
#' elements of each sequence are truncated in order to return a model
#' which can be compared with other Markovian models of visible order maxOrder.
#'
#' @param y the dataset (\code{\link{march.Dataset-class}}) from which to construct the model.
#' @param order the order of the constructed model.
#' @param maxOrder the maximum visible order among the set of Markovian models to compare.
#' @param mtdg flag indicating whether the constructed model should be a MTDg using a different transition matrix for each lag (value: \emph{TRUE} or \emph{FALSE}).
#' @param MCovar vector of the size Ncov indicating which covariables are used (0: no, 1:yes)
#' @param init the init method, to choose among \emph{best}, \emph{random} and \emph{weighted}.
#' @param deltaStop the delta below which the optimization phases of phi and Q stop.
#' @param llStop the ll increase below which the EM algorithm stop.
#' @param maxIter the maximal number of iterations of the optimisation algorithm (zero for no maximal number).
#'
#' @author Ogier Maitre, Kevin Emery
#' @example tests/examples/march.mtd.construct.example.R
#' @seealso \code{\link{march.Mtd-class}}, \code{\link{march.Model-class}}, \code{\link{march.Dataset-class}}.
#' @export
march.mtd.construct <- function(y,order,maxOrder=order,mtdg=FALSE,MCovar=0,init="best", deltaStop=0.0001, llStop=0.01, maxIter=0){
  order <- march.h.paramAsInteger(order)
  if(order<1){
    stop('Order should be greater or equal than 1')
  }
  maxOrder <- march.h.paramAsInteger(maxOrder)

  if( order>maxOrder ){
    stop("maxOrder should be greater or equal than order")
  }

  ySave <- y
  y <- march.dataset.h.filtrateShortSeq(y,maxOrder+1)
  y <- march.dataset.h.cut(y,maxOrder-order)
  
  if(sum(MCovar)>0){
    placeCovar <- which(MCovar==1)
    y@cov <- y@cov[,,placeCovar,drop=FALSE]
    y@Ncov <- as.integer(sum(MCovar))
    y@Kcov <- y@Kcov[placeCovar]
  }else{
    y@cov <- array(0,c(1,1))
    y@Ncov <- as.integer(0)
    y@Kcov <- 0
  }
  
  is_constrained <- TRUE # if FALSE the model is unconstrained and the constraints given by Eq. 4 are activated instead of those given by Eq. 3 (Berchtold, 2001, p. 380)
  is_mtdg <- mtdg
  init_method <- init
  # 1.1. Choose initial values for all parameters
  # nt <- BuildArrayNumberOfDataItems(y) This is no longer needed, as y now contains the T field
  c <- BuildContingencyTable(y,order)
  u <- CalculateTheilU(y,order,c)
  init_params <- InitializeParameters(u=u,init_method=init_method,c=c,is_mtdg=is_mtdg,m=y@K,order=order,kcov=y@Kcov,ncov=y@Ncov)
  phi <- init_params$phi
  q <- init_params$q
  S <- init_params$S


  # 1.2. Choose a value for delta and a criterion to stop the algorithm
  if(y@Ncov>0){
    maxkcov <- max(y@Kcov)
  }else{
    maxkcov=0
  }
  delta <- array(0.1,dim=c(1,y@K+1+maxkcov*y@Ncov)) # a different delta is used for each of the m+1 sets of parameters (vector phi + each row of the transition matrix Q)
  delta_stop <- deltaStop
  ll_stop <- llStop
  # 2. Iterations
  # 2.1. Reestimate the vector phi by modifying two of its elements
  # 2.2. Reestimate the transition matrix Q by modifying two elements of each row

  # 3. End criterion
  # 3.1. If the increase of the LL since the last iteration is greater than the stop criterion, go back to step 2
  # 3.2. Otherwise, end the procedure
  i0_il <- BuildArrayCombinations(y@K,order,y@Kcov,y@Ncov)
  n_i0_il <- BuildArrayNumberOfSequences(y,order)
  q_i0_il <- BuildArrayQ(y@K,l=order,i0_il=i0_il,n_i0_il = n_i0_il,q=q,kcov=y@Kcov,ncov=y@Ncov,S=S)
  new_ll <- CalculateLogLikelihood(n_i0_il=n_i0_il,q_i0_il=q_i0_il,phi=phi)
 

  iter <- 0
  while (TRUE){
    if( maxIter>0 && iter>=maxIter ){ break }
    else{ iter <- iter+1 }


    ll <- new_ll
    # "Reestimate the vector phi by modifying two of its elements." (p. 383)
    pd_phi <- PartialDerivativesPhi(n_i0_il=n_i0_il,q_i0_il=q_i0_il,l=order,phi = phi,ncov=y@Ncov)
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
        res_opt_q <- OptimizeQ(q=q,j=j,pd_q=pd_q,delta=delta,delta_stop=delta_stop,ll=new_ll,n_i0_il=n_i0_il,q_i0_il=q_i0_il,phi=phi,i0_il = i0_il,k=y@K,l=order,g=g, kcov=y@Kcov, ncov=y@Ncov, S=S)
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
    
    if(y@Ncov>0){
      tot=0
      for (i in 1:y@Ncov){
        tot=tot+1
        pd_s<-PartialDerivativesS(CCov=order+i,ColVT=order+i+1,k=y@K,kcov=y@Kcov[i],n_i0_il=n_i0_il,q_i0_il=q_i0_il,i0_il=i0_il,phi)
        for (g in 1:y@Kcov[i]){
          res_opt_S<-OptimizeS(order,y@K,y@Kcov,y@Ncov,S,g,phi,order+i,new_ll,pd_s[g,],delta[1+y@K+(tot-1)*maxkcov+g],delta_stop,n_i0_il,i0_il,q_i0_il,q,i)
          new_ll<-res_opt_S$ll
          S<-res_opt_S$S
          delta[1+y@K+(tot-1)*maxkcov+g]<-res_opt_S$delta
        }
      }
    }
    
    if (new_ll - ll < ll_stop){ break }
  }
  
  #Computation of the number of zeros
  nbZeros <- 0
  nbZeros <- length(which(q==0))+length(which(phi==0))
  if(sum(MCovar)>0){
    for(i in 1:sum(MCovar)){
      nbZeros <- nbZeros+length(which(S[[i]]==0))
    }
  }
  ll <- as.numeric(ll)
  
  #Computation of the full transition matrix
  #Remind that y is the dataset whithout the unused covariates. We put at the beginning
  #y in ySave
  tmCovar <- 1
  if(y@Ncov>0){
    for (i in 1:y@Ncov){
      tmCovar <- tmCovar*y@Kcov[i]
    }
  }
  NSS <- matrix(0,y@K^(order+1)*tmCovar,1)
  ValT <- BuildArrayCombinations(y@K,order,y@Kcov,y@Ncov)
  ProbT <- BuildArrayQ(y@K,order,ValT,NSS,q,y@Kcov,y@Ncov,S)
  l <- GMTD_tm_cov(order,y@K,phi,matrix(ProbT,y@K^(order+1)*tmCovar,order+y@Ncov))
  
  RA <- l$CRHOQ  
  
  #Return the final model
  new("march.Mtd",RA=RA,order=order,Q=q,phi=phi,S=S,MCovar=MCovar,ll=ll,y=ySave,dsL=sum(y@T-order),nbZeros=nbZeros)
}
