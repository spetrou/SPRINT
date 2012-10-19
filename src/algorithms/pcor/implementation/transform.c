
/*!
 *  @brief      Implementation of correlation transform functions.
 *  @details    This file contains the implementation of the correlation
 *              transformation functions described in 'cor_func.h'.
 *  @see        cor_func.h
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */

#include "transform.h"

/*!
 *  @brief  Comparison function for QuickSort.
 *  @param  v1  Pointer to the index of the first value to be compared.
 *  @param  v2  Pointer to the index of the second value to be compared.
 *  @return If v1 < v2 then -1, if v1 > v2 then 1, otherwise 0.
 */
static int func_cmp(const void *v1,
                    const void *v2);

/* Pointer used with QuickSort */
static double *gp_arr;

void transform_to_expected(double *data, const int rows, const int columns,
                           double *Sxx_vector)
{
    int i, j;
    double *mean_value_vector;

    // Make sure value is valid
    if ( columns <= 0 || rows <= 0) {
        printf("\nRow/column size passed to \"transform_to_ranks\" function invalid.\n");
        printf("\nValues remain as they were before the call.\n");
        return;
    }

    mean_value_vector = (double *)malloc(sizeof(double) * rows);

    if( mean_value_vector == NULL ) {
        printf("\nMemory allocation in \"transform_to_expected\" function failed.\n");
        printf("\nValues remain as they were before the call.\n");
        return;
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

    free(mean_value_vector);
}



/* ********************************************* *
 *  Function to transform input values to ranks  *
 * ********************************************* */
void transform_to_ranks(double *data, const int rows, const int columns)
{
    int i, j, *R = NULL;

    // Make sure value is valid
    if ( columns <= 0 || rows <= 0) {
        printf("\nRow/column size passed to \"transform_to_ranks\" function invalid.\n");
        printf("\nValues remain as they were before the call.\n");
        return;
    }

    // Set the pointer to the beginning of the data array
    gp_arr = data;

#pragma omp parallel shared(gp_arr) private(i, j, R) default(none)
{
    // Allocate space for R
    if ( (R = (int *)malloc(columns * sizeof(int))) == NULL ) {
        printf("\nMemory allocation in \"transform_to_ranks\" function failed.\n");
        printf("\nValues remain as they were before the call.\n");
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

