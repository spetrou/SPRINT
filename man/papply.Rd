\name{papply}
\alias{papply}
\title{Parallel Apply}
\description{ 
Parallel apply function used to perform the same operation over all the 
elements of data objects like matrices, data frames or lists. This function 
provides a parallel implementation of both the apply() and lapply() functions
from the core of the R programming language. This parallel implementation 
papply() has been extended to accept an ff object as input which allows the 
processing of much larger data stored on disk. 
}
\usage{
papply(data, fun, margin = 1, out_filename = NULL)
}
\arguments{
    \item{data}{array, list or ff object}
    \item{fun}{function to be applied}
    \item{margin}{vector indicating which elements of the matrix the function
                  will be applied to. The default value is 1 and indicates the
                  rows, 2 indicates the columns and the parameter is ignored 
                  if data is a list.}
    \item{out_filename}{string, name of the result file when input is an ff
                        object}
}
\details{
The function to be applied can be supplied to papply() either as a function 
name or as a function definition. When only the function name is provided, 
the package implementing the function has to be loaded before the SPRINT 
library is initialised in order to ensure that the name is recognised by all 
the processes involved in the computation.
}
\seealso{
    \code{\link{apply}}
    \code{\link{lapply}}
    \code{\link[ff]{ff}}
    \code{\link{SPRINT}}
}
\author{
University of Edinburgh SPRINT Team
    \email{sprint@ed.ac.uk}
    \url{www.r-sprint.org}
    }
\keyword{utilities}
\keyword{interface}