\name{SPRINT}
\alias{sprint-package}
\docType{package}
\title{Overview of SPRINT}
\description{
SPRINT (Simple Parallel R INTerface) is a parallel framework for R. It 
provides a High Performance Computing (HPC) harness which allow R scripts to 
run on HPC clusters. SPRINT contains a library of selected R functions that 
have been parallelized. Functions are named after the original R function with
the added prefix 'p', i.e. the parallel version of cor() in SPRINT is called 
pcor(). These parallelized functions are written in C and MPI. Call to these 
functions are included directly in standard R scripts. 

The following functions are implemented in SPRINT beta 0.3.0:
 - ptest
 - pterminate
 - pcor
 - pmaxT
 - ppam
 
See the User Guide and Release Notes in the sprint folder or the SPRINT web 
page: http://www.r-sprint.org for more information.
}
\keyword{utilities}
\keyword{interface}