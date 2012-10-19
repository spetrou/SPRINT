
/*!
 *  @brief      Implementation of Kendall correlation functions.
 *  @details    This file contains the implementation of the functions
 *              related to the Kendall correlation.
 *  @see        kendall.h
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */

#include "kendall.h"

void kendall(double *data, const int row_x, double *result,
             const int rows, const int columns) {

    int i, j, z, k=0, X_sign, Y_sign;
    double *x, *y, denominator;

    // Make sure value is valid
    if ( columns <= 0 || rows <= 0 ) {
        printf("\nRow/column size passed to \"pearson\" function invalid.\n");
        printf("\nValues remain as they were before the call.\n");
        return;
    }

    // Compute the denominator of the final corellation coefficient
    denominator = (double)((columns * (columns - 1)) / 2);

    // Select the row to be compared against the dataset
    x = &data[row_x*columns];

#pragma omp parallel for shared(result, x, data, denominator) \
                         private(j, y, z, k, Y_sign, X_sign) \
                         default(none)
    for (i=0; i < rows; i++) {
        // Set main diagonal to 1
        if ( i == row_x ) {
            result[i] = 1.0;
            continue;
        }

        y = &data[i*columns];
        k = 0;

        // Execute Kendall's rank correlation algorithm for each pair of rows.
        for(j=0; j<columns; j++) {
            for(z=j+1; z<columns; z++) {
                X_sign = (x[j] > x[z]? 1:-1);
                Y_sign = (y[j] > y[z]? 1:-1);
                k += (X_sign * Y_sign);
            }
        }
        
        result[i] = (double)k / denominator;
    }
}

