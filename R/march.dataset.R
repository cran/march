# TODO: Add comment
# 
# Author: Ogier Maitre
###############################################################################

#' Extract a sequence from a dataset.
#' 
#' 
#' @param y is the dataset from which one sequence of data is extracted
#' @param i is the number of sequence to extract
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
#' @return a \code{\link[=march.Dataset-class]{march.Dataset}} object containing the data from the file found at \emph{filename}, using separator
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
#' The function creates a \code{\link[=march.Dataset-class]{march.Dataset}} from a \emph{dataframe} or a \emph{matrix}, where each row (resp. column) represents 
#' an independent data series when \emph{MARGIN} is 2 (resp. 1).
#' 
#' @param dataframe A \code{\link{data.frame}} containing the dataset.
#' @param MARGIN The dimension of the matrix/data.frame that contains the sequences (resp 1 for the column, 2 for the rows).
#' @param weights If specified, contains the weight of each sequence.
#' @param missingDataRep If specified, the symbol representing a missing data.
#' @return A \code{\link[=march.Dataset-class]{march.Dataset}} object containing the data contructed from the matrix or data.frame.
#' @example tests/examples/march.dataset.loadFromDataFrame.example.R
#' @author Ogier Maitre
#'
#' @export
march.dataset.loadFromDataFrame <- function( dataframe, MARGIN=2,weights=NA,missingDataRep=NA ){
  
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
  
  new("march.Dataset",yRaw=as.matrix(r),T=T,y=y,dictionary=dictionary,N=N,K=K,weights=weights)
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
  
  new("march.Dataset",y=seqs,K=y@K,N=as.integer(count),T=unlist(nT),weights=unlist(nW),dictionary=y@dictionary)
}

march.dataset.h.cut <- function(y,start){
  seqs <- list()
  nT <- list()
  nW <- list()
  
  for( i in 1:y@N ){
    seqs[[i]] <- y@y[[i]][(1+start):y@T[i]]
    nT[[i]] <- y@T[i]-start
    nW[[i]] <- y@weights[i]
  }
  
  new("march.Dataset",y=seqs,K=y@K,N=y@N,T=unlist(nT),weights=unlist(nW),dictionary=y@dictionary)
}