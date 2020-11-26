# This file is part of March.
# It contains functions for the computation of confidence intervals.
# 
# March is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyrigth 2014-2020, Ogier Maitre, Kevin Emery, André Berchtold 
# andre.berchtold@unil.ch 




##############################################################################
##############################################################################
# Tool functions used fo CIs

march.bailey.ci <- function(ni,N,alpha,k){
  
  if (ni>0){
    pMinus <- max((ni-(1/8))/(N+(1/8)),0)
    pPlus <- min((ni+(7/8))/(N+(1/8)),1)
    
    B <- qchisq(1-(alpha/k),1)
    C <- B/(4*N)
    
    t1 <- (sqrt(pMinus)-sqrt(C*(C+1-pMinus)))^2
    t2 <- (sqrt(pPlus)+sqrt(C*(C+1-pPlus)))^2
    t3 <- (C+1)^2
    
    LowerBd <- t1/t3
    UpperBd <- t2/t3
    
    list(LowerBd,UpperBd) 
  } else {
    list(0,0)
  }
  
}


march.ci.h.A <- function(n,ni){
  (ni-1/8)/(n+1/8)
}

march.ci.h.B <- function(n,ni){
  (ni+7/8)/(n+1/8)
}

march.ci.h.C <- function(alpha,K,n){
  qchisq(1-alpha/K,df=1)/(4*n)
}

# Get the d^2*n value as defined in Thompson 1987 p43 paper 
# ("Sample Size for Estimating Multinomial Proportions").
# The value is available for some value of alpha (as defined in the vector a).
# If the alpha value is not defined, the default value (for alpha=0.05) is
# returned.
#
# Parameters : 
#   alpha : the significance level (alpha)
#
# Returns :
#  the d^2*n associated with the alpha value, if defined, or the one
#	associated with alpha=0.05.
#
march.ci.h.d2n <- function(alpha){
  a <- c( 0.5,0.4,0.3,0.2,0.1,0.05,0.025,0.02,0.01,0.005,0.001,0.0005,0.0001 )
  d2n <- c( 0.44129,0.50729,0.60123,0.74739,1.00635,1.27359,1.55963,1.65872,1.96986,2.28514,3.02892,3.33530,4.11209)
  
  id <- which(a==alpha)
  
  if( !length(id) ){
    warning("alpha value has not been found, using 0.05 instead.",call.=FALSE)
    id <- which(a==0.05)
  }
  
  d2n[id]
}

# the expectation of Z_t(g)
# march.mtd.h.z <- function(mtd,y,t,g){
#   s <- 0
#   for( k in 1:mtd@order ){
#     s <- s+mtd@phi[k]*mtd@Q[1,y@y[t-k],y@y[t]]
#   }
#   if(y@Ncov>0){
#     for(j in 1:y@Ncov)
#       s <- s+mtd@phi[order+j]*mtd@S[[j]][cov[n,t,j],y@y[t]]
#   }
#   mtd@phi[g]*mtd@Q[1,y@y[t-g],y@y[t]]/s
# }
# 
# 
# # the weight coefficient for the g-th lag 
# march.mtd.h.l <- function(mtd,y,g){
#   s <- 0
#   for( i in 1:y@N ){
#     ys <- march.dataset.h.extractSequence(y,i)
#     for( t in march.h.seq(mtd@order+1,ys@N)){
#       s <- s+ march.mtd.h.z(mtd,ys,t,g)
#     }
#   }
#   s
# }

# Matrix containing the estimations of the number of data used to compute each element of Q
# march.mtd.h.n <- function(mtd,y){
#   nki_0 <- array(0,c(mtd@y@K,mtd@y@K))
#   for( k in 1:mtd@y@K ){
#     for( i0 in 1:mtd@y@K){
#       s <- 0
#       for( i in 1:y@N ){
#         ys <- march.dataset.h.extractSequence(y,i)
#         for( t in march.h.seq(mtd@order+1,ys@N) ){
#           if( ys@y[t]==i0 ){
#             for( g in 1:mtd@order ){
#               if( ys@y[t-g]==k)
#                 s <- s+march.mtd.h.z(mtd,ys,t,g,n)
#             }
#           }      
#         }
#       }
#       nki_0[k,i0] <- s
#     }
#   }
#   nki_0
# }

#Matrices containing the estimations of the number of data used to compute each element of Q
march.mtd.h.n <- function(mtd,y,is_mtdg){
  
  #Initialization
  if(is_mtdg==FALSE){
    nki_0 <- array(0,c(mtd@y@K,mtd@y@K))
  }else{
    nki_0 <- array(0,c(mtd@order,mtd@y@K,mtd@y@K))
  }
  numcov <- 0
  if(sum(mtd@MCovar)>0){
    placeCovar <- which(mtd@MCovar==1)
    numcov <- array(0,c(sum(mtd@MCovar),max(y@Kcov[placeCovar]),mtd@y@K))
  }
  
  #Computation of the matrices
  for(n in 1:y@N){
    ys <- march.dataset.h.extractSequence(y,n)
    for( t in march.h.seq(mtd@order+1,ys@N)){
      #Computation of the denominator (see p.9 Confidence Intervals for Markovian Models, Berchtold)
      tot <- march.mtd.h.z.tot(mtd,ys,t,n,is_mtdg)
      for(ord in 1:mtd@order){
        if(is_mtdg==FALSE){
          nki_0[ys@y[t-ord],ys@y[t]] <- nki_0[ys@y[t-ord],ys@y[t]]+mtd@phi[ord]*mtd@Q[1,ys@y[t-ord],ys@y[t]]/tot
        }else{
          nki_0[ord,ys@y[t-ord],ys@y[t]] <- nki_0[ord,ys@y[t-ord],ys@y[t]]+mtd@phi[ord]*mtd@Q[ord,ys@y[t-ord],ys@y[t]]/tot
        }
      }
      if(sum(mtd@MCovar)>0){
        for(i in 1:sum(mtd@MCovar)){
          numcov[i,y@cov[n,t,placeCovar[i]],ys@y[t]] <- numcov[i,y@cov[n,t,placeCovar[i]],ys@y[t]]+mtd@phi[mtd@order+i]*mtd@S[[i]][y@cov[n,t,placeCovar[i]],ys@y[t]]/tot
        }
      }
    }
  }
  list(nki_0=nki_0,numcov=numcov)
}

march.mtd.h.z.tot <- function(mtd,ys,t,n,is_mtdg){
  
  s <- 0
  placeCovar <- which(mtd@MCovar==1)
  
  if(is_mtdg==FALSE){
    for( k in 1:mtd@order ){
      s <- s+mtd@phi[k]*mtd@Q[1,ys@y[t-k],ys@y[t]]
    }
  }else{
    for( k in 1:mtd@order ){
      s <- s+mtd@phi[k]*mtd@Q[k,ys@y[t-k],ys@y[t]]
    }
  }
  
  if(sum(mtd@MCovar)>0){
    for(j in 1:sum(mtd@MCovar))
      s <- s+mtd@phi[mtd@order+j]*mtd@S[[j]][mtd@y@cov[n,t,placeCovar[j]],ys@y[t]]
  }
  s
}


alphat <- function(d,s){
  a <- array(0,c(s@N,2))
  
  a[1,] <- d@RB[,1,s@y[1]]*d@Pi[1,1,]
  for( t in 2:s@N ){
    for( g in 1:d@M ){
      for( i in 1:d@M ){
        a[t,g] <- a[t,g]+a[t-1,i]*d@A[i,g]*d@RB[g,1,s@y[t]]
      }
    }
  }
  a
}

betat <- function(d,s){
  b <- array(0,c(s@N,2))
  
  b[s@N,] <- c(1,1)
  for( t in (s@N-1):1 ){
    for( g in 1:d@M ){
      for( i in 1:d@M ){
        b[t,i] <- b[t,i]+d@A[i,g]*d@RB[g, 1 ,s@y[t+1]]*b[t+1,g]
      }
    }
  }
  b
}

march.dcmm.exp_z <- function(d,alpha,beta,t,g,n){
  L <- sum(alpha[n,])
  
  alpha[t,g]*beta[t,g]/L
}

#
# march.dcmm.test <- function(d){
#   PoA <-  array(0,c(d@M^d@orderHC))
#   PoCt <- array(0,c(d@M^d@orderHC))
#   for( i in 1:d@y@N ){
#     # number of point for A
#     s <- march.dataset.h.extractSequence(d@y,i)
#     
#     a <- march.dcmm.forward(d,s)
#     b <- march.dcmm.backward(d,s)
#     epsilon <- march.dcmm.epsilon(d,s,a$alpha,b$beta,a$l,b$l,a$LL)
#     gamma <- march.dcmm.gamma(d,s,a$alpha,b$beta,a$l,b$l,a$LL,epsilon)
#     
#     PoA <- PoA+colSums(gamma[d@orderHC:(d@y@T[i]),])
#     
#     if( d@M>2 ){
#       PoCt <- PoCt+colSums(gamma[1:(d@M-1),])+colSums(gamma[d@M:d@y@T[i],])  
#     }
#     else{
#       PoCt <- PoCt+gamma[1,]+colSums(gamma[d@M:d@y@T[i],])  
#     }
#   }
#   
#   PoC <- array(0,c(1,d@M))
#   for( i in 0:(d@M-1) ){
#     PoC[i+1] <- sum(PoCt[(i*d@orderHC+1):((i+1)*d@orderHC)])
#   } 
#   list(PoA,PoC)
# }