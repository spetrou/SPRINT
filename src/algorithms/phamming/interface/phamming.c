/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2012 The University of Edinburgh                          *
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
#include <string.h>
#include "../../../sprint.h"
#include "../../../functions.h"

extern int hamming(int n,...);

/* **************************************************************** *
 *  Accepts information from R gets response and returns it.        *
 *                                                                  *
 * **************************************************************** */
int phamming(char **x,
             char **out_filename,
             int *sample_width,
             int *number_of_samples)
{

  int response;
  int worldSize, worldRank;
  
  enum commandCodes commandCode;
  
  // Check that MPI is initialized
  MPI_Initialized(&response);
  if (response) {
    DEBUG("\nMPI is init'ed in phamming\n");
  } else {
    DEBUG("\nMPI is NOT init'ed in phamming\n");
    return(-1);
  }
  
  // Get size and rank from communicatorx
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  
  // Broadcast command to other processors
  commandCode = PHAMMING;
  MPI_Bcast(&commandCode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // Call the partitioning around medoids function
  response = hamming(4, *x, *sample_width, *number_of_samples, *out_filename);
  
  return response;

}
