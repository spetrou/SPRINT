##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

# The R stub for the sprint.setOpenMPThreads function.
# This does some rudimentary argument type checking and then hands off to
# the C stub.

## consistent error / warning messages; could use for internationalization
..msg <- list(error =
              c(non.integer = "Input vector must contain integer (consider use of as.integer())",
                non.numeric = "Input vector must be a numeric matrix",
                no.vector = "Input argument must be a vector of integers",
                no.argument = "The input vector is set to NULL",
                no.negative = "The integers need to be positive numbers greater than 0"
                ), warn = c()
              )

sprint.setOpenMPThreads <- function(
  threadPerMPIProcess = NULL   # input numerical vector
)
  {

    # Check that a vector is given (size does not matter for now)
    if(is.null(threadPerMPIProcess)) {
        stop(..msg$error["no.argument"])
    }

    # Check if values are numeric
    if (!is.numeric(threadPerMPIProcess)) {
        stop(..msg$error["non.numeric"])
    }

    # Check that input argument is a vector
    # (scalars are vector of length 1)
    if(!is.vector(threadPerMPIProcess)) {
        stop(..msg$error["no.vector"])
    }

    # Check that arguments are integers
    if (!is.integer(threadPerMPIProcess)) {
        stop(..msg$error["non.integer"])
    }

    # Check that the integers are positive numbers
    for(i in c(1:length(threadPerMPIProcess))) {
        if(threadPerMPIProcess[i] <= 0) {
            stop(..msg$error["no.negative"])
        }
    }

    # Call C interface
    return_val <- .Call("setSprintOpenMPThreads", threadPerMPIProcess,
                        length(threadPerMPIProcess))

    # Return values from the interface have meaning.
    #  0    -->     success
    # -1    -->     MPI is not initialized
    # -2    -->     Only the master process exists, no workers
    if ( return_val == 0 ) {
        result <- TRUE
    } else {

        if ( return_val == -1 )
            warning(paste("MPI is not initialized. Function is aborted.\n"))
        if ( return_val == -2 )
            warning(paste("No worker processes exist. Function setOpenMPThreads() is aborted.\n"))
        result <- FALSE
    }

    return(result)
}

