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
extern int setThreads(int n,...);

/* Local function to check the size of the input vector.
 * There are three cases:
 *
 * (1) Vector has more integers than the number of MPI processes.
 * --------------------------------------------------------------
 * If the vector has more integers that are MPI processes
 * warn the user that the extra numbers will be ignored.
 *
 * (2) Vector has less integers than the number of MPI processes.
 * --------------------------------------------------------------
 * If the vector has less elements than MPI processes we extend
 * the vector and fill in the missing elements with 1s.
 *
 * (3) Vector has exactly as many integers as MPI processes.
 * ---------------------------------------------------------
 * Do nothing.
 *
 */
static int* checkAndInitVector(int *numThreadsVec, const int vecLength,
                               const int worldSize, int* freeMemory);

/* ******************************************************** *
 *  The stub for the R side of the setSprintOpenMPThreads   *
 *  command. This command will broadcast the number of      *
 *  threads that each MPI process should use for            *
 *  computations.
 * ******************************************************** */

SEXP setSprintOpenMPThreads(SEXP numThreads, SEXP vecLength)
{
    SEXP result;
    int response, worldSize, freeMemory = 0;
    enum commandCodes commandCode;
    int* inputVectorPtr = NULL;

    MPI_Initialized(&response);
    if (response) {
        DEBUG("MPI is init'ed in setSprintOpenMPThreads\n");
    } else {

        DEBUG("MPI is NOT init'ed in setSprintOpenMPThreads\n");
        PROTECT(result = NEW_INTEGER(1));
        INTEGER(result)[0] = -1;
        UNPROTECT(1);

        return result;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    // broadcast command to other processors
    commandCode = SETSPRINTOPENMPTHREADS;
    MPI_Bcast(&commandCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Find the length of the vector. It needs to have AT LEAST
    // as many integers as they are MPI processes or else the
    // scatter will fail.
    // If we have more element we can ignore the extra integers,
    // if we have fewer then we need to pad the array with 1s.
    inputVectorPtr = checkAndInitVector(INTEGER_POINTER(numThreads),
                                        INTEGER(vecLength)[0],
                                        worldSize, &freeMemory);

    if(inputVectorPtr == NULL)
    {
        DEBUG("Error occured in checkAndInitVector function");
        response = -3;
    }
    else
    {
        response = setThreads(1, inputVectorPtr);

        if(freeMemory)
        {
            free((void *)inputVectorPtr);
        }
    }

    // Create response integer and set it to the "response" value
    PROTECT(result = allocVector(INTSXP, 1));
    INTEGER(result)[0] = response;

    UNPROTECT(1);
    return result;
}


static int* checkAndInitVector(int *numThreadsVec, const int vecLength,
                               const int worldSize, int* freeMemory)
{

    int i = 0;

    // If same size then return.
    // Integers sanity checks are performed in R level.
    // Checks for: (a) non-negative, (b) integers
    if(vecLength == worldSize)
    {
        return numThreadsVec;
    }
    else if(vecLength > worldSize)
    {
        printf("\nInput array contains more integers than available MPI proceses\n"
               "The first %d integers are going to be used.\nThe extra values are"
               " going to be ignored", worldSize);
        return numThreadsVec;

    }
    else if(vecLength < worldSize)
    {

        int* newVector = NULL;

        newVector = (int *)malloc(worldSize * sizeof(int));
        if(newVector == NULL) {
            printf("\nMemory allocation error. Function returning.");
            return NULL;
        }

        // Copy values to new vector.
        for(i = 0; i < vecLength; ++i)
        {
            newVector[i] = numThreadsVec[i];
        }

        // Fill in remaining values with 1s.
        // Index variable "i" points to the right element.
        for(; i < worldSize; ++i)
        {
            newVector[i] = 1;
        }

        // Raise flag to trigger memory free.
        *freeMemory = 1;

        return newVector;
    }
    else
    {
        // There is really no need for this... just for safety
        return NULL;
    }
}

