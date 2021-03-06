% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/march.tools.R
\name{march.write}
\alias{march.write}
\title{Save a march.Model}
\usage{
march.write(filename, object, force = FALSE)
}
\arguments{
\item{filename}{a path to the file where to write the model (absolute or relative to the current directory).}

\item{object}{the model to write.}

\item{force}{if TRUE and if the file pointed by the filename path already exists, overwrite it.

 @return invisible TRUE if the model has been written into the file pointed by filename, invisible FALSE otherwise.}
}
\description{
Save a march.Model into a file pointed by \emph{filename}. The save will fails
if the file already exists unless force has been set to TRUE.
}
