
/*!
 *  @brief      Implementation of Pearson correlation functions.
 *  @details    This file contains the implementation of the functions
 *              related to the Pearson correlation.
 *  @see        pearson.h
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */

#include "pearson.h"

void pearson(double *data, const int row_x, double *result,
             const int rows, const int columns, const double *Sxx_vector)
{
    int i, j;
    double sxy = 0.0;
    double *x, *y;

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

//double pearson_XY(double *dataMatrixX,
//                  double *dataMatrixY,
//                  int row_x,
//                  int row_y,
//                  int size)
//{
//    int i;
//    double sxy = 0.0;
//    double r;
//    double *x, *y;
//
//    x = &dataMatrixX[row_x*size];
//    y = &dataMatrixY[row_y*size];
//
//    // correlation
//    for (i=0; i<size; i++) {
//        sxy += x[i] * y[i];
//    }
//
//    sxy /= (size-1);
//
//    // This *should* set r=NAN if sxx*syy = 0 *and*
//    // sxy = 0 - this replicates R's behaviour
//    r = sxy / sqrt(Sxx_vector[row_x]*Syy_vector[row_y]);
//
//    return r;
//}

