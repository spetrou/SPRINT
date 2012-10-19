#ifndef _KENDALL_H
#define _KENDALL_H

/*!
 *  @brief      Header file for the Kendall algorithm functions.
 *  @details    The functions of this file are used for computing the Kendall
 *              correlation algorithmm.
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */

#include <stdio.h>
#include <math.h>

#include <omp.h>

/*!
 *  @brief  Implementation of the Kendall correlation algorithm.
 *  @param  data        Pointer to the first data array.
 *  @param  first_row   Index of the row to compute correlation for.
 *  @param  result      Pointer to the results array.
 *  @param  rows        Count of rows.
 *  @param  columns     Count of columns.
 */
void kendall(double *data, const int first_row, double *result,
             const int rows, const int columns);

#endif

