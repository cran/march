### 
# march.nbLLData <- function(rt,k,order){
# 	ones(1,k^order)%*%rt%*%ones(k,1);
# }

march.h.ones <- function(l,c){
  matrix(data=1,nrow=l,ncol=c);
}


#
# Generate a sequence of number from "from" to "to", incremented by "by".
# This function is basically the same as standard seq, except it retun
# a null sequence (numeric(0)) whenever from is less than to and by is positive
# or when from is more than to and by is negative.
# The function seq return an error in these cases.
#
# Parameters:
#	  from: the sequence start.
# 	to: the sequence end.
#	  by: the sequence step.
#
# Returns :
#	a sequence from "from" to "to", using step "by", if possible, numeric(0) otherwise.
#
march.h.seq <- function( from, to, by=1 ){
  if( by>0 ){
    if( to<from ){ numeric(0) }
    else{ seq(from,to,by) }
  }
  else if( from<to ){ numeric(0) }
  else{ seq(from,to,by) }
}

# Convert a parameter from real to integer type and warn the user if during the process something was lost.
march.h.paramAsInteger <- function( p ){
  pInt <- as.integer(p)
  if( pInt!=p ){
    warning(sprintf("%s parameter has been truncated to %d, as it should be an integer\n",toString(substitute(p)),pInt),call.=FALSE)
  }
  pInt
}

# Allow to print a matrix without row and col names.
march.h.printMatrix <- function( m ){
  write.table(format(m,justify="right"),row.names=F,col.names=F,quote=F)
}

#' Save a march.Model
#' 
#' Save a march.Model into a file pointed by \emph{filename}. The save will fails
#' if the file already exists unless force has been set to TRUE.
#' 
#' @param filename a path to the file where to write the model (absolute or relative to the current directory).
#' @param object the model to write.
#' @param force if TRUE and if the file pointed by the filename path already exists, overwrite it.
#' 
#'  @return invisible TRUE if the model has been written into the file pointed by filename, invisible FALSE otherwise.
#'
#' @export
march.write <- function(filename,object,force=FALSE){
  status <- FALSE
  c <- class(object)[1]
  if( c=="march.Indep"){ status <- TRUE}
  if( c=="march.Mc"){ status <- TRUE}
  if( c=="march.Mtd"){ status <- TRUE}
  if( c=="march.Dcmm"){ status <- TRUE}
  
  if( status ){
    if( file.exists(filename) ){
      if( force==FALSE ){
        warning(sprintf("This file already exists (\"%s\"), use force=TRUE to overwrite.",filename),call.=FALSE)
        return(invisible(FALSE))
      }
      else{
        warning(sprintf("This file already exists (\"%s\"), but will be overwritten.",filename),call.=FALSE)
      }
    }
    
    saveRDS(object,filename)
    return(invisible(TRUE))
  }
  else{
    warning(toString(sprintf("%s should be a march.Model",substitute(p))))
    return(invisible(FALSE))
  }
}

#' Load a march.Model.
#' 
#' Load a march.Model from a file pointed by \emph{filename} and check that the model is valid.
#' 
#' @param filename the path where load the mode
#' @return the march.Model contained into the file pointed by filename if it exists and contains a valid model.
#' 
#' @export
march.read <- function( filename ){
  if( file.exists(filename) ){
    object <- readRDS(filename)
    status <- FALSE
    c <- class(object)[1]
    if( c=="march.Indep"){ status <- TRUE}
    if( c=="march.Mc"){ status <- TRUE}
    if( c=="march.Mtd"){ status <- TRUE}
    if( c=="march.Dcmm"){ status <- TRUE}
    
    if( status==FALSE ){
      warning(sprintf("File (\"%s\") does not contain a march.Model object.",filename),call.=FALSE)
      return(NULL);
    }
    return(object)
  }
  else{
    warning(sprintf("File (\"%s\") does not exist.",filename),call.=FALSE)
    return(NULL);
  }
}

# Label rows and columns of a MC matrix (RT or RC), 
# with variable states (dictionary index), in order to print it.
march.h.mc.printableMatrix <- function(s,order,K){
  colnames(s) <- 1:K
  rn <- array(0,c(dim(s)[1]))
  for( i in 1:dim(s)[1]){
    n <- ""
    tmp <- i-1
    for( j in 1:order ){
      v <- tmp%%K
      tmp <- floor(tmp/K)
      n <- paste(n,as.character(v+1))
    }
    rn[i] <- n
  }
  rownames(s) <- paste(rn,":")
  s
}