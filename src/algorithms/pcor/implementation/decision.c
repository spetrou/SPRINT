
#include "cor_func.h"
/*@-compdef@*/
/*@-boolops@*/
/*@-mustdefine@*/

extern int _cor_MPI_procs;
extern int _cor_OMP_procs;


/* **************************************************************** *
 *  Depending on the number of MPI processes and OMP threads this   *
 *  function will decide which approach to follow in order to       *
 *  get best performance.                                           *
 * **************************************************************** */

int decision(int method, double *values, int rows, int columns, double *results) {

    int worldSize=1, worldRank=0;
    int init_flag=-1, i;
    static int fake_argc=1;
    double *Sxx;

    omp_set_num_threads(_cor_OMP_procs);

    // Execute correlation algorithm for each pair of rows...
    if ( method == 1 || method == 3 ) {

        Sxx = (double *)malloc(sizeof(double) * rows);
        if( Sxx == NULL ) {
            ERR("Error allocating Sxx vector.\n\n");
            return(-1);
        }

        // If Spearman's algorithm then transform to ranks first
        if ( method == 3 )
            transform_to_ranks(values, rows, columns);

        // Compute expected values
        transform_to_expected(values, rows, columns, Sxx);

        for(i=0; i<rows; i++) {
            pearson(values, i, &results[i*rows], rows, columns, Sxx);
        }

        // Free memory for expected values
        free(Sxx);
    }

    if ( method == 2 ) {

        transform_to_ranks(values, rows, columns);

        for(i=0; i<rows; i++) {
            kendall(values, i, &results[i*rows], rows, columns);
        }
    }

    return(0);
}


