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

#ifndef _CORRELATION_KERNEL_H
#define _CORRELATION_KERNEL_H

enum CORRELATION_METHOD {
    PEARSON = 0,
    SPEARS  = 1,
    KENDALL = 2
};

// FIXME:   Move this in 'sprint.h'?? It should really be shared
//          to everyone. 'static' should not be needed, not it needs it
//          to avoid multiple definition.
static const int MASTER_PROCESS_RANK = 0;   /*! MPI rank of 'master' */

/*!
 *  @brief      The correlation execution kernel.
 *  @details    Implementation of correlation kernel for the 'pcor' function.
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 *
 *
 *  The idea of the kernel
 *  ----------------------
 *
 *  The kernel is implemented based on the 'task farm' parallel concept.
 *  The 'master' process dispatches 'work' to the 'worker' processes and
 *  waits for them to finish before sending them the next piece of 'work'.
 *  When all 'work' pieces are finished the 'master' gathers the results
 *  from the 'worker' processes. The gathering is implemented as a collective
 *  parallel file write operation.
 *
 *  The steps performed by the 'master' process are:
 *
 *  Initialization
 *  --------------
 *  At this step the master process sends the first 'work' piece to each
 *  'worker'. In order to do this it compares the number of available 'work'
 *  pieces to the number of available 'worker' processes:
 *      (a) If we have more 'worker' processes than 'work', then a few of
 *          the processes will not receive any work. The 'master' sends
 *          all active 'worker' processes their first piece and informs
 *          the rest that they will need to just wait.
 *      (b) If we have equal or more pieces of 'work' than 'worker' processes
 *          then the 'master' sends each one their first piece of 'work'.
 *
 *  Main loop
 *  ---------
 *  The 'master' process waits for 'worker' processes to finish and dispatch
 *  more work to them, if available, or inform them that the 'work' is done
 *  and they should move on to the gathering step. When a 'worker' process
 *  informs the 'master' that a 'work' is done, it also sends back an error
 *  code. In case of an error, the 'master' stops the dispatch and moves
 *  onto a 'cleanup' stage.
 *
 *
 *  The steps performed by the 'worker' process are:
 *
 *  Initialization
 *  --------------
 */

/*!
 *  @brief  The main computation kernel imlpementation.
 *  @param rank
 *  @param size
 *  @param matrix
 *  @param dataMatrixX
 *  @param dataMatrixY
 *  @param width
 *  @param height
 *  @param out_filename
 *  @param distance_flag
 *  @param method_flag
 *  @return On success 0, otherwise non-zero value.
 */
int correlationKernel(int rank, int size, double *matrix, double* dataMatrixY,
                      int width, int height, char *out_filename, int distance_flag,
                      int method_flag);

#endif

