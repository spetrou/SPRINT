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

# The R stub for the sprint.MPI function family.
# These functions are implemented in order to expose MPI information to the
# R environment.

## consistent error / warning messages; could use for internationalization
..msg <- list(error =
              c(not.initialized = "MPI environment not initialized"),
              warn = c()
              )

sprint.MPI_WORLD_SIZE <- function()
{
    # Call C interface
    return_val <- .Call("mpi_getWorldSize")

    # Return values from the interface have meaning.
    # -1    -->     MPI is not initialized
    if ( return_val == -1 ) {
        warning(..msg$error["non.numeric"])
    }

    return(return_val)
}

