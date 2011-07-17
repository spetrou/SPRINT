/**************************************************************************
 *  Copyright © 2008,2009 The University of Edinburgh                     *
 *                                                                        *
 *  This program is free software: you can redistribute it and/or modify  *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  any later version.                                                    *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful,       *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this program. If not, see <http://www.gnu.or/licenses/>.   *
 *                                                                        *
 **************************************************************************/

#include <Rdefines.h>
#include "../../../sprint.h"
#include "../../../functions.h"
extern int boot(int n,...);
void Carray2Rmatrix(double * array, SEXP matrix, int r, int c);
void Rmatrix2Carray(SEXP matrix, int * array, int r, int c); 
void Rmatrix2CDBLarray(SEXP matrix, double * array, int r, int c); 

/* ******************************************************** *
 *  The stub for the R side of a very simple test command.  *
 *  Simply issues the command and returns.                  *
 * ******************************************************** */

//SEXP pboot(SEXP data, SEXP statistic, SEXP ind, SEXP lt0, SEXP varg){
SEXP pboot(SEXP scenario,...){
   
  SEXP result;
  double *func_results;
  int response, worldSize;
  enum commandCodes commandCode;
  int scene = asInteger(scenario);
  va_list ap;
  int c;
  // get the function arguments common to all scenarios 
  // 1, R, lt0, vargs, strdata, strstatistic
  va_start(ap, scenario); 
  int r = asInteger(va_arg(ap, SEXP)); // the number of replications to perform
  int ltn = asInteger(va_arg(ap, SEXP)); // the number of results the statistical function returns
  SEXP varg = va_arg(ap, SEXP);
  SEXP data = va_arg(ap, SEXP);
  SEXP statistic = va_arg(ap, SEXP);

  MPI_Initialized(&response);
  if (response) {
      DEBUG("MPI is init'ed in ptest\n");
  } else {

      DEBUG("MPI is NOT init'ed in ptest\n");
      PROTECT(result = NEW_INTEGER(1));
      INTEGER(result)[0] = -1;
      UNPROTECT(1);

      return result;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
 
  // broadcast command to other processors
  commandCode = PBOOT;
  MPI_Bcast(&commandCode, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // intialise the memory to store the results
  func_results = (double *)malloc(sizeof(double) * r * ltn); 

  int * find;
  int * wind;
  SEXP f, w, Spred;
  int m; 
  int * pred;

  switch(scene) {   
    case 1:;
      SEXP rangen, mle;
      rangen = va_arg(ap, SEXP);
      mle = va_arg(ap, SEXP);
      response = boot(1, func_results, r, ltn, varg, CHAR(STRING_ELT(data,0)), translateChar(PRINTNAME(statistic)), rangen, mle);
      break; 
    case 2:;
      f = va_arg(ap, SEXP); 
      c = ncols(f); // number of columns in the index
      find = (int *)malloc(sizeof(int) * r * c);
      Rmatrix2Carray(f, find, r, c);
      response = boot(2, func_results, r, ltn, varg, CHAR(STRING_ELT(data,0)), translateChar(PRINTNAME(statistic)), c, find);
      free(find);
      break; 
    case 3:;
      f = va_arg(ap, SEXP); 
      Spred = va_arg(ap, SEXP); 
      c = ncols(f); // number of columns in the index
      m = ncols(Spred); 
      find = (int *)malloc(sizeof(int) * r * c);
      pred = (int *)malloc(sizeof(int) * r * m);
      Rmatrix2Carray(f, find, r, c);
      Rmatrix2Carray(Spred, pred, r, m);
      response = boot(3, func_results, r, ltn, varg, CHAR(STRING_ELT(data,0)), translateChar(PRINTNAME(statistic)), c, find, pred, m);
      free(find);
      free(pred);
      break; 
    case 4:
      w = va_arg(ap, SEXP);
      c = ncols(w); // number of columns in the index
      wind = (int *)malloc(sizeof(int) * r * c);
      Rmatrix2Carray(w, wind, r, c);
      response = boot(4, func_results, r, ltn, varg, CHAR(STRING_ELT(data,0)), translateChar(PRINTNAME(statistic)), c, wind);
      free(wind);
      break; 
    case 5:;
      w = va_arg(ap, SEXP);
      Spred = va_arg(ap, SEXP);
      c = ncols(w); // number of columns in the index
      m = ncols(Spred);
      wind = (int *)malloc(sizeof(int) * r * c);
      pred = (int *)malloc(sizeof(int) * r * m);
      Rmatrix2Carray(w, wind, r, c);
      Rmatrix2Carray(Spred, pred, r, m);
      response = boot(5, func_results, r, ltn, varg, CHAR(STRING_ELT(data,0)), translateChar(PRINTNAME(statistic)), c, wind, pred, m);
      free(wind);
      free(pred);
      break; 
    case 8:
      ;// work around for gcc bug
      // retrieve function arguments 
      SEXP ind = va_arg(ap, SEXP);
      c = ncols(ind); // replications are the number of columns in the index
      // convert the ind from (horrible) SEXP format to C array
      int * cind;
      cind = (int *)malloc(sizeof(int) * r * c);
      Rmatrix2Carray(ind, cind, r, c);
      
      // sending everything to the implementation function
      response = boot(8, func_results, r, ltn, varg, CHAR(STRING_ELT(data,0)), translateChar(PRINTNAME(statistic)), c, cind);
      free(cind);
      break; // end of scenario 8
    default:
      break; 
  }// end of switch
  va_end(ap);
  // Turn the array passed back from the implementation into 
  // a SEXP object that can be returned to R.
  PROTECT(result = allocMatrix(REALSXP,r ,ltn)); // t.star <- matrix(NA, sum(R), lt0)
  Carray2Rmatrix(func_results, result, r, ltn);
  free(func_results);
  UNPROTECT(1); 
  return result;
}

void Rmatrix2Carray(SEXP matrix, int * array, int r, int c){
  int i,j;
  int count = 0;
  int indx = 0; 
  for(i=0;i<r;i++){
    indx = i;
    for(j=0;j<c;j++){
      array[count] = INTEGER(matrix)[indx];
      count++;
      indx += r;
    }
  }
}

void Rmatrix2CDBLarray(SEXP matrix, double * array, int r, int c){
  int i,j;
  int count = 0;
  int indx = 0;
  for(i=0;i<r;i++){
    indx = i;
    for(j=0;j<c;j++){
      array[count] = REAL(matrix)[indx];
      count++;
      indx += r;
    }
  }
}


void Carray2Rmatrix(double * array, SEXP matrix, int r, int c){
  int i, j; 
  int count = 0;
  int indx = 0;
  for(i=0; i < c; i++) {
    indx = i;
    for (j=0; j <r; j++){
      REAL(matrix)[count]= array[indx];
      count++;
      indx += c;
    }
  }
}
