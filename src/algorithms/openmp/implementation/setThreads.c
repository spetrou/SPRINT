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

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

/**
 * A command used to set the number of OpenMP threads in all
 * MPI hosts.
 **/

int setThreads(int n,...)
{
    int worldSize;
    int worldRank;
    va_list ap;

    // An array of integers (as many elements as MPI processes)
    // The sanity checks are perfmormed by the "interface".
    int* numThreads = NULL;
    int scatterNumThreads = 1;

    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

    // master processor gathers results form slaves
    if (worldRank == 0) {

        // Get input variables
        va_start(ap, n);
        numThreads = va_arg(ap, int*);
        va_end(ap);

        // Scatter the number of threads each MPI process
        // should be configured for.
        MPI_Scatter(numThreads, 1, MPI_INTEGER,
                    &scatterNumThreads, 1, MPI_INTEGER,
                    0, MPI_COMM_WORLD);
    } else {

        MPI_Scatter(NULL, 1, MPI_INTEGER,
                    &scatterNumThreads, 1, MPI_INTEGER,
                    0, MPI_COMM_WORLD);
    }

    omp_set_num_threads(scatterNumThreads);

    printf("\nProcess %d will spawn %d threads", worldRank, scatterNumThreads);
    return 0;
}

