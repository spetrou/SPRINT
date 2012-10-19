#ifndef _TRANSFORM_H
#define _TRANSFORM_H

/*!
 *  @brief      Header file for functions used to transform arrays.
 *  @details    The functions of this file transform the data of an input array
 *              to the format needed by the correlation algorithm. Two functions
 *              are implemented. One for transforming an array to 'expected'
 *              values and one for transforming it to ranks.
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Inlcude the header files for the two parallel libraries
#include <omp.h>

/*!
 *  @brief      Compute parameters used by Pearson method
 *  @details    It will transform the input array to "expected values"
 *              and also create vectors for the mean values and variances
 *  @param  data
 *  @param  rows
 *  @param  columns
 *  @param  Sxx_vector
 */
void transform_to_expected(double *data, const int rows, const int columns,
                           double *Sxx_vector);


/*!
 *  @brief  Function to transform the input dataset into ranks
 *          
 *  @param  data
 *  @param  rows
 *  @param  columns
 *  @return 
 */
void transform_to_ranks(double *data, const int rows, const int columns);

#endif

