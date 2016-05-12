# This file is part of March
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
# Copyrigth 2014, Ogier Maitre, Andre Berchtold
# ogier.maitre@epfl.ch andre.berchtold@unil.ch
#


# the imported methods which needed for march are added to the 
# NAMMESPACE file with the following Code.
#' @importFrom stats C qchisq rnorm runif
#' @importFrom utils read.table write.table
#' @importFrom methods new


###############################################################################
# Show methods here, are the implementation of generic methods to print
# models to the standard output. These methods are implemented for all
# models (indep, mc, mtd, dcmm).
###############################################################################

# show method for model (abstract) object
march.model.show <- function(object){
  cat("ll : ",object@ll,"\n")
  cat("dsL: ",object@dsL,"\n")
}

# show method for indep object
march.indep.show <- function(object){
  cat(march.name(object))
  cat("\n")

  s <- array(0,c(1,length(object@indP)))
  colnames(s) <- 1:object@y@K

  s[1,] <- object@indP
  rownames(s) <- c("indP")
  print(s)
  cat("\n")

  s[1,] <- object@indC
  rownames(s) <- c("indC")
  print(s)
  cat("\n")

  march.model.show(object)
}

# show method for mc object
march.mc.show <- function(object){
  cat(march.name(object))
  cat("\nRC : \n")

  s <- march.h.mc.printableMatrix(object@RC,object@order,object@y@K)
  print(s)
  cat("\n")

  cat("RT : \n")
  s <- march.h.mc.printableMatrix(object@RT,object@order,object@y@K)
  print(s)
  cat("\n")
  cat("order : ",object@order,"\n")

  march.model.show(object)
}

# show method for mtd object
march.mtd.show <- function(object){
  cat(march.name(object))
  cat("\n")
  if( dim(object@Q)[1]==1 ){

    cat("Q : \n")
    for( i in 1:object@y@K ){
      cat(sprintf("%.4f", object@Q[1,i,]))
      cat("\n")
    }
  }
  else{
    for( g in 1:dim(object@Q)[1]){
      cat(sprintf("Q%d : \n",g))
      for( i in 1:object@y@K ){
        cat(sprintf("%.4f", object@Q[g,i,]))
        cat("\n")
      }
    }
  }

  cat("\n")
  cat("phi : ")
  cat(sprintf("%.4f",object@phi))
  cat("\n\n")
  cat("order : ",object@order,"\n")
  march.model.show(object)
}

# show method for dcmm object
march.dcmm.show <- function(object){
  cat(march.name(object))
  cat("\nRA : \n")
  s <- march.h.mc.printableMatrix(march.dcmm.h.compactA(object),object@orderHC,object@M)
  print(s)
  cat("\n")

  cat("RB : \n")
  for( i in 1:object@M ){
    cat(i,":\n");
    if( object@y@K^object@orderVC==1 ){
      s <- t(data.frame(object@RB[i,,]))
      s <- march.h.mc.printableMatrix(s,1,object@y@K)
    }
    else{
      s <- march.h.mc.printableMatrix(object@RB[i,,],object@orderVC,object@y@K)
    }
    print(s)
  }

  cat("\nPi : \n")
  for( i in 1:object@orderHC ){
    for( j in 1:(object@M^(i-1))){
      cat(sprintf("%.4f",object@Pi[i,j,]))
      cat("\n")
    }
    cat("\n")
  }

  march.model.show(object)
}

# this part describes how a call to show method (print) should be
# redirected to the rigth method, depending on the object considered.
setMethod(f="show",signature="march.Indep",definition=march.indep.show)
setMethod(f="show",signature="march.Mc",definition=march.mc.show)
setMethod(f="show",signature="march.Mtd",definition=march.mtd.show)
setMethod(f="show",signature="march.Dcmm",definition=march.dcmm.show)

###############################################################################
# nbParams methods here, are the implementation of generic methods to obtain
# the number of parameters that have been used during the model construction.
# These methods are implemented for all models (indep, mc, mtd, dcmm).
###############################################################################

# nbParams method for model (abstract) object
march.model.nbParams <- function(object){
  0
}

# nbParams method for mc object
march.mc.nbParams <- function(object){
  object@y@K^object@order*(object@y@K-1)-object@nbZeros+length(which(rowSums(object@RT)==0))
}

# nbParams method for indep object
march.indep.nbParams <- function(object){
  object@y@K-1
}

# nbParams method for mtd object
march.mtd.nbParams <- function(object){
  dim(object@Q)[1]*object@y@K*(object@y@K-1) + object@order-1-object@nbZeros+length(which(rowSums(object@Q)==0));
}

# nbParams method for dcmm object
march.dcmm.nbParams <- function(object){
  nbPPi <- 0
  # get the number of parameters for Pi :
  for( t in 1:object@orderHC ){
    nbPPi <- nbPPi + object@M^(t-1)*(object@M-1)
    for( i in 1:object@M^(t-1)){
      for( j in 1:object@M){
        if( object@Pi[t,i,j]==0 ){ nbPPi <- nbPPi-1 }
      }
    }
    if( all(object@Pi==0)){
      nbPPi <- nbPPi+1
    }
  }

  # get the number of parameters for A
  RA <- march.dcmm.h.compactA(d=object)
  nbPA <- object@M^object@orderHC*(object@M-1)-sum(RA==0)+sum(rowSums(RA)==0)

  # get the number of parameters for RB
  nbPRB <- 0
  for( i in 1:object@M ){
    nbPRB <- nbPRB+object@y@K^object@orderVC*(object@y@K-1)-sum(object@RB[1,,]==0)+sum(colSums(object@RB)==0)
  }
  nbPA+nbPPi+nbPRB
}

# This part create the generic method and describe how a call to this generic
# has to be redirected to the rigth method, according to the considerd object.
setGeneric(name="march.nbParams",def=function(object)march.model.nbParams(object))
setMethod(f="march.nbParams",signature="march.Indep",definition=march.indep.nbParams)
setMethod(f="march.nbParams",signature="march.Mc",definition=march.mc.nbParams)
setMethod(f="march.nbParams",signature="march.Mtd",definition=march.mtd.nbParams)
setMethod(f="march.nbParams",signature="march.Dcmm",definition=march.dcmm.nbParams)


###############################################################################
# thompson allows to compute confidence intervals according to a given model,
# using thompson's confidence interval method describe into: Thompson, S.K. (1987)
# "Sample size for estimating multinomial proportions," American Statistician, 41, 42-46.
# Adaptation to markov models is described into : TODO
###############################################################################

#' Thompson Confidence Intervals for a march.Model.
#'
#' Compute the confidence intervals using Thompson's formula on a march.Model
#' object. See Thompson SK (1987) Sample size for estimating multinomial proportions,
#' American Statistician 41:42-46, for details.
#'
#'
#' @param object the march.Model object on which compute the confidence intervals.
#' @param alpha the significance level among : 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.025, 0.02, 0.01, 0.005, 0.001, 0.0005, 0.0001.
#'
#' @return A list of half-length confidence intervals for each probability distribution of the considered model.
#' @author Ogier Maitre
#' @example tests/examples/march.thompson.example.R
#' 
#' @export
march.thompson <- function(object,alpha){}

march.model.thompson <- function(object,alpha){
  warning("Confidence interval with thompson formula cannot be computed for abstract class \"model\", check the parameters of the call to march.thompson")
}

march.mc.thompson <- function(object,alpha){
  d2n <- march.ci.h.d2n(alpha)
  d <- array(NA,object@y@K^object@order)

  for( i in 1:length(d)){
    s <- sum(object@RT[i,])
    if( s>0 ){
      d[i]<-sqrt(d2n/s)
    }

  }
  d
}

march.indep.thompson <- function(object,alpha){
  d2n <- march.ci.h.d2n(alpha)

  sqrt(d2n/object@dsL)
}

#
# TODO : Handle multi-sequences by computing the number of data
#  	 influencing each distribution and summing this number.
#
march.mtd.thompson <- function( object, alpha){
  d2n <- march.ci.h.d2n(alpha)

  s <- 0
  for( i in 1:object@order ){
    s<-s+march.mtd.h.l(object,object@y,i)
  }

  list(phi=sqrt(d2n/s),Q=sqrt(d2n/rowSums(march.mtd.h.n(object,object@y))))
}

march.dcmm.thompson <- function(object,alpha){

  ys <- march.dataset.h.extractSequence(object@y,1)
  alpha<-march.dcmm.forward(object,ys)
  beta<-march.dcmm.backward(object,ys)

  C
}

#This part create the generic method and describe how a call to this generic
#has to be redirected to the rigth method, according to the considered object.
setGeneric(name="march.thompson",def=function(object,alpha)march.model.thompson(object,alpha))

#' This method is called with the object "march.Indep" and the aplha "numeric" and
#' provides it to the march.thompson function.
#' @param object contains the name of the model.
#' @param alpha contains the Type I error
setMethod(f="march.thompson",signature=signature(object="march.Indep",alpha="numeric"),definition=march.indep.thompson)

#' This method is called with the object "march.Mc" and the aplha "numeric" and
#' provides it to the march.thompson function.
#'
#' @param object contains the name of the model.
#' @param alpha contains the Type I error
setMethod(f="march.thompson",signature=signature("march.Mc",alpha="numeric"),definition=march.mc.thompson)

#' This method is called with the object "march.Mtd" and the aplha "numeric" and
#' provides it to the march.thompson function.
#'
#' @param object contains the name of the model.
#' @param alpha contains the Type I error
setMethod(f="march.thompson",signature=signature("march.Mtd",alpha="numeric"),definition=march.mtd.thompson)

#' This method is called with the object "march.Dcmm" and the aplha "numeric" and
#' provides it to the march.thompson function.
#'
#' @param object contains the name of the model.
#' @param alpha contains the Type I error
setMethod(f="march.thompson",signature=signature("march.Dcmm",alpha="numeric"),definition=march.dcmm.thompson)

#' march.Model name.
#' 
#' Generate a name for the march model contained in the given \emph{object}.
#'
#' @param object contains the name of the model(Independence model, MTD,...). 
#' @author Ogier Maitre & Andre Berchtold
#' @example tests/examples/march.name.example.R
#' @export
march.name <- function(object){}

march.model.name <- function(object){ "abstract model" } # should never be called
march.indep.name <- function(object){ "Independence" }
march.mc.name <- function(object){ sprintf("MC(%d)",object@order) }
march.mtd.name <- function(object){
  if( dim(object@Q)[1]==1 ){
    sprintf("MTD(%d)",object@order)
  }
  else{
    sprintf("MTDg(%d)",object@order)
  }

}

march.dcmm.name <- function(object){
  if( object@orderVC==0 ){ sprintf("Hmm(%d)",object@orderHC) }
  else{ sprintf("Dcmm(%d,%d)",object@orderHC,object@orderVC) }
}


#This part create the generic method and describe how a call to this generic
#has to be redirected to the right method, according to the considered object.
setGeneric(name="march.name",def=function(object)march.model.name(object))

#' This method is called with the object "march.Indep" and provides it 
#' to the march.name function.
#'
#' @param object contains the name of the model.
setMethod(f="march.name",signature=signature(object="march.Indep"),definition=march.indep.name)

#' This method is called with the object "march.MC" and provides it 
#' to the march.name function.
#'
#' @param object contains the name of the model.
setMethod(f="march.name",signature=signature(object="march.Mc"),definition=march.mc.name)

#' This method is called with the object "march.Mtd" and provides it 
#' to the march.name function.
#'
#' @param object contains the name of the model.
setMethod(f="march.name",signature=signature(object="march.Mtd"),definition=march.mtd.name)

#' This method is called with the object "march.Dcmm" and provides it 
#' to the march.name function.
#'
#' @param object contains the name of the model.
setMethod(f="march.name",signature=signature(object="march.Dcmm"),definition=march.dcmm.name)

#quote=FALSE,digits = getOption("digits")

#' march.Model Summary.
#'
#' Print key values for the current model.
#' 
#' @param object can contain the results of any model computed  using march
#' @param ... should indicate any additional parameter passed to the function
#' @author Ogier Maitre & Andre Berchtold
#' @export
march.summary <- function(object,...){
  v <- array(NA,c(1,4))
  rownames(v) <- march.name(object)
  colnames(v) <- c(gettext("ll"),gettext("param"),gettext("BIC"),gettext("AIC"))
  v[1,]<-c(object@ll,march.nbParams(object),march.BIC(object),march.AIC(object))

  v
}

#setMethod('summary',"march.Model",march.summary)
#setMethod(f="summary",signature=signature(object="march.Mc"),definition=march.summary)
# setMethod(f="Summary",signature=signature(object="march.Mtd"),definition=march.summary)
# setMethod(f="Summary",signature=signature(object="march.Dcmm"),definition=march.summary)
