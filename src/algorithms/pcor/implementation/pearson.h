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
 *  @return The number of coefficients computed. Negative on error.
 */
int pearson_row(const double *data, const int row_x, double *result,
                const int rows, const int columns,
                const double *Sxx_vector);


/*!
 *  @brief      Pearson correlation function for one input arrays.
 *  @details    Function to compute the pearson correlation coefficient
 *              between rows of the same array. One coefficiant
 *              is computed and returned.
 *  @param  data        Pointer to the beginning of the array.
 *  @param  row_x       Index to the first row of the array.
 *  @param  row_y       Index to the second row of the array.
 *  @param  size        Length of the rows.
 *  @param  Sxx_vector  Pointer to the expected values vector.
 *  @return The correlation coefficient between the two rows.
 */
double pearson(const double *data,
               const int row_x,
               const int row_y,
               const int size,
               const double *Sxx_vector);

/*!
 *  @brief      Pearson correlation function for two input arrays.
 *  @details    Function to compute the pearson correlation coefficient
 *              between rows of two distinct arrays. One coefficiant
 *              is computed and returned.
 *  @param  dataMatrixX     Pointer to the beginning of the first array.
 *  @param  dataMatrixX     Pointer to the beginning of the second array.
 *  @param  row_x           Index of the row of first array.
 *  @param  row_x           Index of the row of second array.
 *  @param  size            Length of the rows. The length should be
 *                          identical for the two rows.
 *  @param  Sxx_vector      Pointer to the expected values vector for
 *                          first array.
 *  @param  Sxx_vector      Pointer to the expected values vector for
 *                          second array.
 *  @return The correlation coefficient between the two rows.
 */
double pearson_XY(const double *dataMatrixX,
                  const double *dataMatrixY,
                  const int row_x,
                  const int row_y,
                  const int size,
                  const double *Sxx_vector,
                  const double *Syy_vector);

#endif

