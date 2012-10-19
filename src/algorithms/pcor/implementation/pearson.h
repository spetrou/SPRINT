#ifndef _PEARSON_H
#define _PEARSON_H

/*!
 *  @brief      Header file for the Pearson algorithm functions.
 *  @details    The functions of this file are used for computing the Pearson
 *              correlation algorithmm.
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */

#include <stdio.h>
#include <math.h>

#include <omp.h>

/*!
 *  @brief      Implementation of the Pearson correlation algorithm.
 *  @details    The single input array implementation of the Pearson
 *              correlation algorithm.
 *  @param  data        Pointer to the first data array.
 *  @param  row_x       Index of row to compute correlation for.
 *  @param  result      Pointer to the results array.
 *  @param  rows        Count of rows.
 *  @param  columns     Count of columns.
 *  @param  Sxx_vector  Pointer to the expected values array.
 *
 */
void pearson(double *data, const int row_x, double *result,
             const int rows, const int columns, const double *Sxx_vector);

#endif

