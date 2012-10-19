
/*!
 *  @brief      Implementation of Pearson correlation functions.
 *  @details    This file contains the implementation of the functions
 *              related to the Pearson correlation.
 *  @see        pearson.h
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */

#include "pearson.h"

void pearson_row(const double *data, const int row_x, double *result,
                 const int rows, const int columns, const double *Sxx_vector)
{
    int i, j;
    double sxy = 0.0;
    const double *x, *y;

    // Make sure value is valid
    if ( columns <= 0 || rows <= 0) {
        printf("\nRow/column size passed to \"pearson\" function invalid.\n");
        printf("\nValues remain as they were before the call.\n");
        return;
    }

    // Set pointer to the row the array is compared against
    x = &data[row_x*columns];

    // Compute the correlation coefficients for the given row
    for (i=0; i < rows; i++) {
        // Set main diagonal to 1
        if ( i == row_x ) {
            result[i] = 1.0;
            continue;
        }

        y = &data[i*columns];

        // Corvariance
        sxy = 0.0;
        for (j=0; j<columns; j++) {
            sxy += x[j] * y[j];   
        }
        sxy /= (columns-1);

        // This *should* set r=NAN if sxx*syy = 0 *and*
        // sxy = 0 - this replicates R's behaviour
        result[i] = sxy / sqrt(Sxx_vector[row_x]*Sxx_vector[j]);
    }
}


double pearson(const double *data,
               const int row_x,
               const int row_y,
               const int size,
               const double *Sxx_vector)
{
    int i = 0;
    double sxy = 0.0;
    double r = 0.0;
    const double *x = NULL, *y = NULL;

    x = &data[row_x*size];
    y = &data[row_y*size];

    // Correlation
    for (i=0; i<size; i++) {
        sxy += x[i] * y[i];   
    }

    sxy /= (size-1);

    // This *should* set NAN if sxx*syy = 0 *and*
    // sxy = 0 - this replicates R's behaviour
    return sxy / sqrt(Sxx_vector[row_x]*Sxx_vector[row_y]);
}


double pearson_XY(const double *dataMatrixX,
                  const double *dataMatrixY,
                  const int row_x,
                  const int row_y,
                  const int size,
                  const double *Sxx_vector,
                  const double *Syy_vector)
{
    int i;
    double sxy = 0.0;
    double r;
    const double *x, *y;

    x = &dataMatrixX[row_x*size];
    y = &dataMatrixY[row_y*size];

    // correlation
    for (i=0; i<size; i++) {
        sxy += x[i] * y[i];
    }

    sxy /= (size-1);

    // This *should* set NAN if sxx*syy = 0 *and*
    // sxy = 0 - this replicates R's behaviour
    return sxy / sqrt(Sxx_vector[row_x]*Syy_vector[row_y]);
}

