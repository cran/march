# TODO: Add comment
# 
# Author: Ogier Maitre
###############################################################################


#' Extract a sequence from a dataset.
#' @param y A sequence of integers.
#' @param i The number of observations to keep.
#' @author Ogier Maitre
#' @export
march.dataset.h.extractSequence <- function(y,i){
  new("march.Sequence",y=y@y[[i]],N=as.integer(y@T[i]),weight=y@weights[i])
}

#' Load a dataset from a file.
#' 
#' The function loads a dataset from a text file, where each row (resp. column) represents 
#' a data series when \emph{MARGIN} is 2 (resp. 1), using the character \emph{sep} as attribute separator. 
#' Each data sequence should be stored in a given column, (resp. row).
#' 
#' @param filename The complete path to the text file containing the dataset.
#' @param MARGIN The dimension of the extracted data.frame that contains the sequences (resp 1 for the column, 2 for the rows).
#' @param sep A caracter used as element separator on a line.
#' @param weights If specified, contains the weight of each sequence.
#' 
#' @return a \code{\link{march.Dataset-class}} object containing the data from the file found at \emph{filename}, using separator
#' \emph{sep}.
#' @author Ogier Maitre
#' #'
#' @export
march.dataset.loadFromFile <- function( filename, MARGIN=2,sep=",",weights=NA){
  r<- t(read.table(filename,sep=sep))
  
  march.dataset.loadFromDataFrame(r,MARGIN=MARGIN,weights)
}

#' Construct a dataset from a data.frame or a matrix.
#' 
#' The function creates a \code{\link{march.Dataset-class}} from a \emph{dataframe} or a \emph{matrix}, where each row (resp. column) represents 
#' an independent data series when \emph{MARGIN} is 2 (resp. 1).
#' 
#' @param dataframe A \code{\link{data.frame}} containing the dataset.
#' @param MARGIN The dimension of the matrix/data.frame that contains the sequences and of the covariates (resp 1 for the column, 
#' 2 for the rows).
#' @param weights If specified, contains the weight of each sequence.
#' @param missingDataRep If specified, the symbol representing a missing data.
#' @param covariates If specified, a three dimensional array of integers, representing the covariates. The data for the i-th covariates
#' should be in [, , i]. If the data are column-wise (respectively row-wise), each table of covariates should be column-wise (respectively row-wise).
#' If we only have one covariate, we can simply pass a two-dimensional array. The covariates should be coded as integers from 1 to the number of possible outputs.
#' @return A \code{\link{march.Dataset-class}} object containing the data contructed from the matrix or data.frame.
#' @example tests/examples/march.dataset.loadFromDataFrame.example.R
#' @author Ogier Maitre
#'
#' @export
march.dataset.loadFromDataFrame <- function( dataframe, MARGIN=2,weights=NA,missingDataRep=NA, covariates=NULL ){
  
  if( is.na(missingDataRep)==FALSE ){ dataframe[dataframe==missingDataRep]<-NA }
  
  if( MARGIN==1 ){
    r<- t(dataframe)
  }
  else{
    r<- dataframe
  }
  N <- dim(r)[1]
  T <- array(0,c(N))
  
  # Get the number of elements per sequence
  # read.table should return a matrix where non existing values
  # are represented by NA
  for( i in 1:N){
    isNa <- is.na(r[i,])
    if( sum(isNa)>0 ){ T[i] <- min(which(isNa))-1 }
    else{ T[i] <- dim(r)[2] }
  }
  
  # creating a factor representation and its dictionary
  yfactorV <- factor(as.vector(t(r)))
  dictionary <- levels(yfactorV)
  
  # determining the number of possible output
  K <- as.integer(max(as.numeric(yfactorV),na.rm=TRUE))
  
  
  cov=array()
  Ncov=as.integer(0)
  Kcov=vector()
  
  #
  if(is.null(covariates)==FALSE){
    if(length(dim(covariates))==3){
    Ncov=as.integer(dim(covariates)[3])
    } else{
      Ncov=as.integer(1)
    }
    cov=array(covariates,c(dim(covariates)[1],dim(covariates)[2],Ncov))
    Kcov=array(0,Ncov)
    if( MARGIN==1 ){
      tcov=array(0,c(dim(covariates)[2],dim(covariates)[1],Ncov))
      for(i in 1:Ncov){
        tcov[,,i]<-t(cov[,,i])
      }
      cov=tcov
    }
  
  
  for(i in 1:Ncov){
    s<-cov[,,i]
    cfactor<-factor(as.vector(t(s)))
    #dico<-levels(cfactor)
    Kcov[i]<-as.integer(max(as.numeric(cfactor),na.rm=TRUE))
  }

  }
  # create a list of values for each sequence
  y <- list(N)
  y[[1]] <- as.numeric(yfactorV[1:T[1]])
  for( i in march.h.seq(2,N) ){
    y[[i]] <- as.numeric(yfactorV[(dim(r)[2]*(i-1))+(1:(T[i]))])
  }
  

  
  # weights
  if( length(weights)==1 && is.na(weights) ){
    weights <- array(1,c(N))
    
  }
  if(length(weights)!=N){
    stop("The number of weights should be equal to the number of sequences")
  }
  new("march.Dataset",yRaw=as.matrix(r),T=T,y=y,dictionary=dictionary,N=N,K=K,weights=weights,Ncov=Ncov,cov=cov,Kcov=Kcov)
}

march.dataset.h.filtrateShortSeq <- function(y,limit){
  seqs <- list()
  nT <- list()
  nW <- list()
  
  count <- 0
  
  for( i in 1:y@N ){
    if( y@T[i]>=limit ){
      count <- count+1
      seqs[[count]] <- y@y[[i]]
      nT[[count]] <- y@T[i]
      nW[[count]] <- y@weights[i]
      
    }
  }
  
  new("march.Dataset",y=seqs,K=y@K,N=as.integer(count),T=unlist(nT),weights=unlist(nW),dictionary=y@dictionary,Ncov=y@Ncov,Kcov=y@Kcov,cov=y@cov)
}

march.dataset.h.cut <- function(y,start){
  seqs <- list()
  nT <- list()
  nW <- list()
  newCov<-y@cov
  
  for( i in 1:y@N ){
    seqs[[i]] <- y@y[[i]][(1+start):y@T[i]]
    nT[[i]] <- y@T[i]-start
    nW[[i]] <- y@weights[i]
  }
  if(y@Ncov>0){
    d<-dim(y@cov)[2]
    newCov<-array(0,c(y@N,d-start,y@Ncov))
    for(i in 1:y@Ncov){
      for(n in 1:y@N){
        newCov[n,,i]<-y@cov[n,(start+1):d,i]
      }
    }
  }
  
  new("march.Dataset",y=seqs,K=y@K,N=y@N,T=unlist(nT),weights=unlist(nW),dictionary=y@dictionary,Ncov=y@Ncov,Kcov=y@Kcov,cov=newCov)
}