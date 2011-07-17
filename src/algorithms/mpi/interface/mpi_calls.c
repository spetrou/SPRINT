/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008,2009 The University of Edinburgh                     *
 *                                                                        *
 *  This program is free software: you can redistribute it and/or modify  *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  any later version.                                                    *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful,       *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this program. If not, see <http://www.gnu.or/licenses/>.   *
 *                                                                        *
 **************************************************************************/

#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../../functions.h"

/* ******************************************************** *
 *  The stub for the R side of the setSprintOpenMPThreads   *
 *  command. This command will broadcast the number of      *
 *  threads that each MPI process should use for            *
 *  computations.
 * ******************************************************** */

SEXP mpi_getWorldSize(void)
{
    SEXP result;
    int response, worldSize;

    MPI_Initialized(&response);
    if (response) {
        DEBUG("MPI is init'ed in mpi_getWorldSize\n");
    } else {

        DEBUG("MPI is NOT init'ed in mpi_getWorldSize\n");
        PROTECT(result = NEW_INTEGER(1));
        INTEGER(result)[0] = -1;
        UNPROTECT(1);

        return result;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    // Create response integer and set it to the world size value
    PROTECT(result = allocVector(INTSXP, 1));
    INTEGER(result)[0] = worldSize;;
    UNPROTECT(1);

    return result;
}

