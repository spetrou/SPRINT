\name{pcor}
\alias{pcor}
\title{Parallel Correlation}
\description{
Sister function to "cor" implementing a parallel Pearson correlation using a 
replicated data strategy.
}
\usage{
pcor(data_x, data_y = NULL, distance = FALSE, caching_ = "mmeachflush", 
     filename_ = NULL)
}
\arguments{
    \item{data_x}{double precision 2D array of data}
    \item{data_y}{NULL or second double precision 2D array of data}
    \item{distance}{boolean, whether the distance or correlation coefficient
                    matrix is returned}
    \item{caching_}{string, either "mmeachflush" or "mmnoflush" select the back-end 
                    caching scheme}
    \item{filename_}{string, name of the result file}
}
\seealso{
   \code{\link{cor}}
   \code{\link{sprint}}
}
\keyword{utilities}
\keyword{interface}