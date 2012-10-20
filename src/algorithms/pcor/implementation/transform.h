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
 *  @brief  Return codes for functions in this file
 */
#define RCODE_SUCCESS                       0
#define RCODE_INPUT_PARAMETER_CHECKS_FAIL   1
#define RCODE_MEMORY_ALLOC_FAIL             2

/*!
 *  @brief      Compute expected values of an array.
 *  @details    It will transform the input array to "expected values".
 *              It also populates the vectors for the mean values and variances.
 *              These arrays need to be passed to the function (pre-allocated).
 *              If any of the pre-processing checks fail, the values of the
 *              array and the vectors are not altered.
 *  @see    transform_to_expected_with_alloc
 *  @param  data                Pointer to the biginning on the input array.
 *  @param  rows                Number of rows of input array.
 *  @param  columns             Number of columns of input array.
 *  @param  Sxx_vector          Vector to store the expected value of each row.
 *  @param  mean_value_vector   Vector to store mean values of each row.
 *  @return Value 0 on success, non-zero otherwise.
 */
unsigned int transform_to_expected(double *data, const int rows,
                                   const int columns, double *Sxx_vector,
                                   double *mean_value_vector);

/*!
 *  @brief      Compute expected values of an array.
 *  @details    It will transform the input array to "expected values".
 *              It also populates the vectors for the mean values and variances.
 *              These arrays will be allocated by the function itself. If they
 *              are not initialized to NULL the function will try to free the
 *              previous memory and allocate the necessary memory according
 *              to the input array dimensions.
 *              ***CAUTION***: Is the caller's responsibility to free the
 *              memory. If any of the pre-processing checks fail, the values
 *              of the array are not altered and the vectors are not touched
 *              in any way. If memory was already attached to them will remain
 *              allocated.
 *  @see    transform_to_expected_with_alloc
 *  @param  data                Pointer to the biginning on the input array.
 *  @param  rows                Number of rows of input array.
 *  @param  columns             Number of columns of input array.
 *  @param  Sxx_vector          Vector to store the expected value of each row.
 *  @param  mean_value_vector   Vector to store mean values of each row.
 *  @return Value 0 on success, non-zero otherwise.
 */
unsigned int transform_to_expected_with_alloc(double *data, const int rows,
                                              const int columns, double *Sxx_vector,
                                              double *mean_value_vector);

/*!
 *  @brief  Function to transform the input array into ranks
 *          
 *  @param  data
 *  @param  rows
 *  @param  columns
 *  @return Value 0 on success, non-zero otherwise.
 */
unsigned int transform_to_ranks(double *data, const int rows,
                                const int columns);

#endif

