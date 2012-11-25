#ifndef _COR_FUNC_H
#define _COR_FUNC_H

/*!
 *  @brief      Header file for functions used 
 *  @details    
 *              
 *              
 *  @author     Savvas Petrou <savvaspetrou@hotmail.com>
 */



#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

// Inlcude the header files for the two parallel libraries
#include <omp.h>
#include <mpi.h>


/* ************************************************************ *
 *  Declare the global variables that are going to be used for  *
 *  controlling the parallelism in all functions.               *
 * ************************************************************ */
//#ifdef _MAIN_FUNC_
//int _cor_MPI_procs = 1;
//int _cor_OMP_procs = 2;
//#endif


/* ******************** *
 *  1  -->  Kendall's   *
 *  2  -->  Pearson's   *
 *  3  -->  Spearman's  *
 * ******************** */
#define COR_METHOD 1

int mpi_boot(void);

/*!
 *  @brief  
 *          
 *  @param  method
 *  @param  values
 *  @param  rows
 *  @param  columns
 *  @param  results
 *  @return 
 */
int decision(int     method,
             double *values,
             int     rows,
             int     columns,
             double *results);


/* ****************************** *
 *  Kendall's correlation method  *
 * ****************************** */
void kendall(                  /*@in@*/ double *data,
                                        int     first_row,
              /*@notnull@*/ /*@out@*/   double *result,
                                        int     rows,
                                        int     columns);


/* ****************************** *
 *  Pearson's correlation method  *
 * ****************************** */
void pearson (                /*@in@*/ double *data,
                                       int     row_x,
               /*@notnull@*/ /*@out@*/ double *result,
                                       int     rows,
                                       int     columns,
                    /*@in@*/ /*@out@*/ double *Sxx_vector);


/* ******************************* *
 *  Spearman's correlation method  *
 * ******************************* */
void spearman (                /*@in@*/ double *data,
                                        int     row_x,
                /*@notnull@*/ /*@out@*/ double *result,
                                        int     rows,
                                        int     columns,
                     /*@in@*/ /*@out@*/ double *Sxx_vector);



// Timer function
double timer(void);


#endif

