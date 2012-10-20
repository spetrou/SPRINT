
/*!
 *  @brief      Implementation of correlation transform functions.
 *  @details    This file contains the implementation of the correlation
 *              transformation functions described in 'cor_func.h'.
 *  @see        cor_func.h
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */

#include "transform.h"
#include "../../../sprint.h"

/*!
 *  @brief  Comparison function for QuickSort.
 *  @param  v1  Pointer to the index of the first value to be compared.
 *  @param  v2  Pointer to the index of the second value to be compared.
 *  @return If v1 < v2 then -1, if v1 > v2 then 1, otherwise 0.
 */
static int func_cmp(const void *v1,
                    const void *v2);

/*! @brief Pointer used in QuickSort */
static double *gp_arr;

unsigned int transform_to_expected_with_alloc(double *data, const int rows,
                                              const int columns, double *Sxx_vector,
                                              double *mean_value_vector)
{
    double *new_expected = NULL;    /*! Temporary variable for memory allocation */
    double *new_mean = NULL;        /*! Temporary variable for memory allocation */
    unsigned int rcode = 0;         /*! Return code */

    // Check dimensions are valid
    if ( columns <= 0 || rows <= 0) {
        ERR("\nRow/column sizes in \"%s\" invalid. rows[%d], columns[%d]\n",
            "transform_to_expected_with_alloc", rows, columns);
        return RCODE_INPUT_PARAMETER_CHECKS_FAIL;
    }

    if(data == NULL) {
        ERR("\nInput array pointer is null in \"%s\"",
            "transform_to_expected_with_alloc");
        return RCODE_INPUT_PARAMETER_CHECKS_FAIL;
    }

    // Allocate memory for vectors.
    new_mean = (double *)malloc(sizeof(double) * rows);
    if(new_mean == NULL) {
        ERR("\nMemory allocation in \"%s\", for \"%s\" failed.\n", 
            "transform_to_expected_with_alloc", "mean_value_vector");
        return RCODE_MEMORY_ALLOC_FAIL;
    }

    new_expected = (double *)malloc(sizeof(double) * rows);
    if(new_expected == NULL) {
        ERR("\nMemory allocation in \"%s\", for \"%s\" failed.\n", 
            "transform_to_expected_with_alloc", "Sxx_vector");
        // Don't forget to free the mean vector's memory
        free(new_mean);
        return RCODE_MEMORY_ALLOC_FAIL;
    }

    Sxx_vector = new_expected;
    mean_value_vector = new_mean;

    rcode = transform_to_expected(data, rows, columns,
                                  Sxx_vector, mean_value_vector);

    // Free memory on failure.
    if(rcode)
    {
        free(Sxx_vector);
        free(mean_value_vector);
        Sxx_vector = NULL;
        mean_value_vector = NULL;
    }

    return rcode;
}


unsigned int transform_to_expected(double *data, const int rows,
                                   const int columns, double *Sxx_vector,
                                   double *mean_value_vector)
{
    int i = 0;  /*! Temporary ndex for rows */
    int j = 0;  /*! Temporary index for columns */

    // Check dimensions are valid
    if ( columns <= 0 || rows <= 0) {
        ERR("\nRow/column sizes in \"%s\" invalid. rows[%d], columns[%d]\n",
            "transform_to_expected", rows, columns);
        return RCODE_INPUT_PARAMETER_CHECKS_FAIL;
    }

    if(data == NULL) {
        ERR("\nInput array pointer is null in \"%s\"", "transform_to_expected");
        return RCODE_INPUT_PARAMETER_CHECKS_FAIL;
    }

    if(mean_value_vector == NULL) {
        ERR("\nMemory allocation in \"%s\", for \"%s\" failed.\n", 
            "transform_to_expected", "mean_value_vector");
        return RCODE_INPUT_PARAMETER_CHECKS_FAIL;
    }

    if(Sxx_vector == NULL) {
        ERR("\nMemory allocation in \"%s\", for \"%s\" failed.\n", 
            "transform_to_expected", "Sxx_vector");
        return RCODE_INPUT_PARAMETER_CHECKS_FAIL;
    }

#pragma omp parallel shared(mean_value_vector, Sxx_vector, data) \
                     private(i, j) \
                     default(none)
{
    // Compute the mean values of each row
#pragma omp for
    for(i=0; i < rows; i++) {
        mean_value_vector[i] = 0.0;
        for(j=0; j < columns; j++) {
            mean_value_vector[i] += data[(i*columns)+j];
        }
        mean_value_vector[i] /= columns;
    }

    // Transform the data array
#pragma omp for
    for(i=0; i < rows; i++) {
        Sxx_vector[i] = 0.0;
        for(j=0; j < columns; j++) {
            data[(i*columns)+j] -= mean_value_vector[i];
            Sxx_vector[i] += (data[(i*columns)+j] * data[(i*columns)+j]);
        }
        Sxx_vector[i] /= (columns-1);
    }
} // end of parallel region

    return RCODE_SUCCESS;
}


unsigned int transform_to_ranks(double *data, const int rows, const int columns)
{
    int i = 0;      /*! Temporary ndex for rows */
    int j = 0;      /*! Temporary index for columns */
    int *R = NULL;  /*! Index vector for QuickSort algorithm */

    // Make sure value is valid
    if ( columns <= 0 || rows <= 0) {
        ERR("\nRow/column sizes in \"%s\" invalid. rows[%d], columns[%d]\n",
            "transform_to_ranks", rows, columns);
        return RCODE_INPUT_PARAMETER_CHECKS_FAIL;
    }

    // Set the pointer to the beginning of the data array
    gp_arr = data;

#pragma omp parallel shared(gp_arr) private(i, j, R) default(none)
{
    // Allocate space for R
    if ( (R = (int *)malloc(columns * sizeof(int))) == NULL ) {
        ERR("\nMemory allocation in \"transform_to_ranks\" function failed.\n");
        ERR("\nValues remain as they were before the call.\n");
        exit(EXIT_FAILURE);
    }

    // Transform the values of each row into ranks
#pragma omp for
    for(i=0; i<rows; i++) {

        // Initialize the index vector
        for(j=0; j<columns; j++)
            R[j] = j + (i * columns);

        // QuickSort the row
        qsort(R, (size_t)columns, sizeof(int), func_cmp);

        // Replace values with ranks
        for(j=0; j<columns; j++)
            gp_arr[R[j]] = (double)(j+1);

    }

    // Free each thread's memory before joining
    free(R);

} // end of parallel region

    return RCODE_SUCCESS;
}


int func_cmp(const void *v1, const void *v2) {

    double f1 = *(gp_arr+*(int *)v1);
    double f2 = *(gp_arr+*(int *)v2);

    if (f1<f2)
        return -1;
    else if (f1>f2)
        return 1;
    else
        return 0;
}

