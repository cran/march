#' Dataset for march package.
#' 
#' This class contains several discrete-valued time series, in a dataset. It contains for each sequence, 
#' its length and weights.
#' 
#' The internal representation uses factor-like representation. The integer values correspond to the words stored
#' into the dictionary vector. Therefor, they are in the interval [1,K].
#' 
#'  @section Slots:
#'  \describe{
#'    \item{\code{yRaw}:}{A matrix of \code{\link{character}} string, describing the content of 
#'    the original dataset or file, if any.}
#'    \item{\code{y}:}{A list of vector of \code{\link{integer}} representing the  each discrete-valued time series
#'    of the dataset, as can be used by the models.}
#'    \item{\code{T}:}{A vector of \code{\link{integer}} values representing the length of each sequence.}
#'    \item{\code{weights}:}{A vector of \code{\link{numeric}} values representing the weight of each sequence.}
#'    \item{\code{K}:}{A \code{\link{integer}} value representing the number of possible ouput and the number
#'    of words stored into the dictionary.}
#'    \item{\code{N}:}{A \code{\link{integer}} value representing the number of sequence.}
#'    \item{\code{Dictionary}:}{A vector of \code{\link{character}} string representing the translation between 
#'    the yRaw and y data. Each character string is stored according to the integer which represents it into y.}
#'  }
#'  @seealso \code{\link{march.dataset.loadFromFile}}, \code{\link{march.dataset.loadFromDataFrame}}
#'  @author Ogier Maitre
setClass("march.Dataset",
         representation(  yRaw="matrix",y="list","T"="vector",weights="vector",
                          K="integer",N="integer",dictionary="vector"))

# This class represents a discrete-valued time serie.
# 
# This class describes a discrete-valued time serie, as extracted from \code{\link{march.Dataset}}.
# It is not meant to be handled by user, but to be used internaly.
# 
# @section Slots:
#  \describe{
#    \item{\code{vector}:}{the discrete-valued time serie.}
#    \item{\code{weight}:}{the weight of the current discrete-valued time serie.}
#    \item{\code{N}:}{The number of samples contained in the serie.}
#  }
setClass("march.Sequence",representation(y="vector",weight="numeric",N="integer"))

#' A basic and virtual march model.
#' 
#' This class describe the basic and virtual model, that every model of the package will extend.
#' This is a virtual class, which is not meant to be handled by user directly.
#' 
#' @seealso The classes that inherit from march.Model are : \code{\link{march.Indep-class}}, \code{\link{march.Mc-class}}, \code{\link{march.Mtd-class}}, \code{\link{march.Dcmm-class}}.
#' 
#'  @section Slots:
#'    \describe{
#'      \item{\code{ll}:}{A \code{\link{numeric}} representing the log-likelihood for this model \emph{w.r.t} its 
#'      construction dataset.}
#'      \item{\code{y}:}{The \code{\link[=march.Dataset-class]{march.DataSet-class}} used to construct the model.}
#'      \item{\code{dsL}:}{A \code{\link{numeric}} representing the number of sample used to construct the model.}
#'      \item{\code{nbZeros}:}{A \code{\link{numeric}} representing the number of zeros created during model construction.}
#'      }
setClass("march.Model",representation(ll="numeric",y="march.Dataset",dsL="numeric",nbZeros="numeric","VIRTUAL"))


#' An independence model.
#' 
#' This class describes an independence model, represented by the probability distribution \emph{indP} of each event
#' and the number of data used to compute each member of the probability distribution. march.Indep inherits
#' from \code{\link{march.Model-class}} and therefore inherits its slots.
#' 
#' @section Slots:
#'  \describe{
#'    \item{\code{indP}:}{A vector of \code{\link{numeric}} representing the model probability distribution.}
#'    \item{\code{indC}:}{A vector of \code{\link{integer}} representing the number of data used to compute each
#'    member of the probability distribution.}
#'    }
#' @seealso \code{\link{march.indep.construct}}, \code{\link{march.Model-class}}.
setClass("march.Indep",contains="march.Model",representation(indP="vector",indC="vector"))

#' A Markov chain of order >= 1.
#' 
#' This class describes a Markov chain of order \emph{order}, represented by matricess RC (transition matrix in reduced form)
#' and RT (number of data points used to compute each transition). march.Mc extends \code{\link{march.Model-class}} class and therefore
#' inherits its slots.
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{RC}:}{A matrix of \code{\link{numeric}} representing the reduced form of the 
#'    transition matrix of the current Markov Chain.}
#'    \item{\code{order}:}{An \code{\link{integer}} representing the order of the current Markov Chain.}
#'    \item{\code{RT}:}{A matrix of \code{\link{integer}} representing the number of sample used to compute each 
#'    transition row of the current RC matrix.}
#'  }
#'  
#' @seealso \code{\link{march.mc.construct}}, \code{\link{march.Model-class}}.
setClass("march.Mc",contains="march.Model",representation(RC="array",order="integer",RT="array"))


#' A Mixture Transition Distribution (MTD) model.
#' 
#' This class describes a Mixture Transition Distribution (MTD) model, represented by its transition matrix Q,
#' its vector phi of lag parameters and its order. march.Mtd extends \code{\link{march.Model-class}} class and therefore
#' inherits its slots.
#' march.Mtd extends \code{\link{march.Model-class}} class and therefore inherits its slots.
#' 
#' The model used here is described into :
#' \itemize{
#'  \item Raftery, A. E. A Model for High-Order Markov Chains. 
#'  JRSS B 47(1985), pp. 528-539.
#'  \item Berchtold, A. Estimation in the mixture transition distribution model.
#'  Journal of Time Series Analysis, 22 (4) (2001), pp. 379-397
#' }

#' 
#'  @section Slots:
#'  \describe{
#'    \item{\code{Q}:}{A matrix of \code{\link{numeric}} representing the transition matrix associated with the 
#'    current MTD model.}
#'    \item{\code{phi}:}{A vector of \code{\link{numeric}} representing the vector of lag parameters.}
#'    \item{\code{order}:}{An \code{\link{integer}} representing the order of the model.}
#'  }
#'    
#' @seealso \code{\link{march.mtd.construct}}, \code{\link{march.Model-class}}.
setClass("march.Mtd",contains="march.Model",
         representation(Q="array",phi="vector",order="integer"))


#' A Double Chain Markov Model (DCMM).
#' 
#' This class describes a Double Chain Markov Model (DCMM) represented by Pi, the probability distributions of the first
#' hidden states; by A, the transition matrix between hidden states; by RB, the transition matrix between sucessive output.
#' march. Dcmm extends \code{\link{march.Model-class}} class and therefore inherits its slots.
#' 
#' The model used here is described in :
#' \itemize{
#'  \item Berchtold, A.: The Double Chain Markov Model. Commun. Stat., Theory Methods 28 (1999), pp. 2569-2589  
#'  \item Berchtold, A.: High-order extensions of the Double Chain Markov Model. Stochastic Models 18 (2002), pp. 193-227.  
#' }
#' 
#' @section Slots:
#'  \describe{
#'    \item{\code{Pi}:}{A 3D matrix of \code{\link{numeric}} representing the probability distribution 
#'    of the first hidden state.}
#'    \item{\code{A}:}{A matrix of \code{\link{numeric}} representing the transition matrix between hidden states.}
#'    \item{\code{RB}:}{A 3D matrix of \code{\link{numeric}} representing the transition matrix between successive 
#'    output, in a reduced form.}
#'    \item{\code{M}:}{An \code{\link{integer}} value representing the number of hidden state.}
#'    \item{\code{orderVC}:}{An \code{\link{integer}} value representing the order of the visible Markov chain.}
#'    \item{\code{orderHC}:}{An \code{\link{integer}} value representing the order of the hidden Markov chain.}
#'  }
#'    
#' @seealso \code{\link{march.dcmm.construct}}, \code{\link{march.Model-class}}.
setClass( "march.Dcmm", 
          representation(Pi="array",A="array",RB="array",M="integer",orderVC="integer",orderHC="integer"),
          contains="march.Model");


# Those are abstract classes used into abstract EA.
# Real EA algorithms should all use implementations of
# these classes.
setClass("march.ea.InitParameters",representation(fct="function"))
setClass("march.ea.EvalParameters",representation(fct="function"))
setClass("march.ea.CrossoverParameters",representation(fct="function"))
setClass("march.ea.MutationParameters",representation(fct="function"))
setClass("march.ea.OptimizingParameters",representation(fct="function"))


setClass("march.ea.Parameters",representation(
  crossoverProb="numeric",
  optimizing="logical",
  
  initParameters="march.ea.InitParameters",
  evalParameters="march.ea.EvalParameters",
  mutationParameters="march.ea.MutationParameters",
  crossoverParameters="march.ea.CrossoverParameters",
  optimizingParameters="march.ea.OptimizingParameters",
  populationSize="integer",
  generation="integer")
)


# those are implementation of the abstract classes defined in march.ea.R
# init, eval and mutation parameters :
setClass("march.dcmm.ea.InitParameters",contains="march.ea.InitParameters",
         representation( AConst="logical",
                         M="integer",
                         K="integer",
                         orderVC="integer",
                         orderHC="integer",
                         y="march.Dataset"
         )
)

setClass("march.dcmm.ea.EvalParameters",contains="march.ea.EvalParameters",
         representation(ds="march.Dataset")		
)

setClass("march.dcmm.ea.MutationParameters",contains="march.ea.MutationParameters",
         representation(pMut="numeric")
)

setClass("march.dcmm.ea.OptimizingParameters",contains="march.ea.OptimizingParameters",
         representation(ds="march.Dataset",iterBw="integer",stopBw="numeric")
)
