/**************************************************************************
 *                                                                        *
 *  SPRINT: Simple Parallel R INTerface                                   *
 *  Copyright © 2008,2009,2010 The University of Edinburgh                *
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

#include <stdarg.h>
#include "../../../sprint.h"
#include <R_ext/Parse.h>

/* These match up with the R interface.  So that nthcdr(args,
 * arguments_t) gives the relevant SEXP value. */
typedef enum _arguments_t {
    X = 0,
    Y,
    XTEST,
    YTEST,
    NTREE,
    MTRY,
    REPLACE,
    CLASSWT,
    CUTOFF,
    STRATA,
    SAMPSIZE,
    NODESIZE,
    MAXNODES,
    IMPORTANCE,
    LOCALIMP,
    NPERM,
    PROXIMITY,
    OOB_PROX,
    NORM_VOTES,
    DO_TRACE,
    KEEP_FOREST,
    COR_BIAS,
    KEEP_INBAG
} arguments_t;

/* When the serialize.c interface is made available to third party
 * packages this could be pulled entirely into the C layer.  As it is
 * we have to go through the interpreter. */
SEXP serialize_form(SEXP form)
{
    SEXP thunk;
    SEXP ret;
    PROTECT(thunk = allocVector(LANGSXP, 3));
    SETCAR(thunk, install("serialize"));
    SETCADR(thunk, form);
    SETCADDR(thunk, R_NilValue);

    ret = eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
    return ret;
}

/* As above. */
SEXP unserialize_form(SEXP form)
{
    SEXP thunk;
    SEXP ret;
    PROTECT(thunk = allocVector(LANGSXP, 2));
    SETCAR(thunk, install("unserialize"));
    SETCADR(thunk, form);

    ret = eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
    return ret;
}

SEXP getListElement(SEXP list, char *str)
{
  SEXP elmt = R_NilValue;
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

void setListElement(SEXP list, char *str, SEXP value)
{
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      SET_VECTOR_ELT(list, i, value);
      break;
    }
}

/* Call the serial randomForest call with the args we were passed in
 * parallel.  We use the size and rank arguments to calculate the
 * subset of the total trees we're going to compute. */
SEXP serial_randomForest(SEXP args, int rank, int size)
{
    SEXP thunk;
    SEXP ret;
    SEXP tmp;
    int length;
    int i;
    int ntree;
    int chunksize;
    int remainder;
#ifdef _PROF
    double t;
    t = MPI_Wtime();
#endif
    length = length(args) + 1;  /* Space for function name too */

    PROTECT(thunk = allocVector(LANGSXP, length));
    SETCAR(thunk, install("randomForest"));

    ntree = INTEGER(coerceVector(CAR(nthcdr(args, NTREE)), INTSXP))[0];

    if ( ntree < size ) {
        UNPROTECT(1);
        return NULL;
    }
    chunksize = ntree/size;
    remainder = ntree - (chunksize * size);
    if ( rank < remainder ) {
        ++chunksize;
    }
    /* We modify this value in place but it's come from R, so we'll need to
     * restore the value afterwards */
    REAL(CAR(nthcdr(args, NTREE)))[0] = (double)chunksize;

    /* Need to keep pointer to head of args list, so just modify a
     * pointer copy here. */
    PROTECT(tmp = args);
    for ( i = 1; i < length; i++ ) {
        /* Because args is a pairlist, but thunk is a linked list, Nil
         * values in thunk confuse the R->C call that the evaluation
         * below carries out.  Basically the Nil value is seen as the
         * end of the list, rather than just a null argument.  To fix
         * this, if we got any Nil arguments through the R call, we
         * replace them with "missing" values.  These get fixed up in
         * the serial randomForest call. */
        if ( CAR(tmp) == R_NilValue ) {
            SET_MISSING(CAR(tmp), 1);
            SETCAR(tmp, R_MissingArg);
        }
        SETCAR(nthcdr(thunk, i), CAR(tmp));
        tmp = CDR(tmp);
    }
    UNPROTECT(1);
    PROTECT(ret = eval(thunk, R_GlobalEnv));

    /* Don't need to send the data back over, so set the "call" part of
     * the list to NULL, we repopulate in the R interface. */
    setListElement(ret, "call", R_NilValue);
    /* The xlevels and ncat variables are the same for a given
     * dataset, so we don't need to transfer these over either if
     * we're not rank 0. */
    if ( rank != 0 ) {
        PROTECT(tmp = getListElement(ret, "forest"));
        setListElement(tmp, "xlevels", R_NilValue);
        setListElement(tmp, "ncat", R_NilValue);
        UNPROTECT(1);
    }
    /* Need to restore saved number of trees value */
    REAL(CAR(nthcdr(args, NTREE)))[0] = (double)ntree;
    UNPROTECT(2);
#ifdef _PROF
    Rprintf("[%d: subforest-gen] %g (%d trees)\n", rank, MPI_Wtime() - t, chunksize);
    fflush(stdout);
#endif
    return ret;
}

/* Load a package.  We need this so that we can load the randomForest
 * package on the slaves. */
static void require_package(char *package_name)
{
    SEXP thunk;
    PROTECT(thunk = allocVector(LANGSXP, 2));
    SETCAR(thunk, install("require"));
    SETCADR(thunk, install(package_name));
    eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
}

/* Combining forests uses the R function pcombine we've just
 * defined which calls out to the random forest library function
 * combine to do most of the heavy lifting. */
SEXP combine_forests(SEXP a, SEXP b)
{
    SEXP thunk;
    SEXP ret;
    PROTECT(thunk = allocVector(LANGSXP, 3));
    SETCAR(thunk, install("pcombine"));
    SETCADR(thunk, a);
    SETCADDR(thunk, b);
    ret = eval(thunk, R_GlobalEnv);
    UNPROTECT(1);
    return ret;
}

/* Do a tree reduction across the processes in COMM combining IN
 * using COMBINE_FN (an associate, commutative, side-effect-free
 * function), writing the result to OUT on ROOT process.  OUT need
 * only be defined on the root process (you can pass NULL on other
 * processes). */
void reduce_combine(SEXP in, SEXP *out, SEXP (*combine_fn)(SEXP, SEXP),
                    int root, MPI_Comm comm)
{
    MPI_Group g;
    MPI_Group doing_comms;
    MPI_Group sitting_out;
    MPI_Comm new;
    MPI_Status status;
    MPI_Request request;
    SEXP ret;
    SEXP tmp;
    int size;
    int rank;
    int neighbour;
    int nbytes;
    int *ranks;
    int i;

    if ( comm == MPI_COMM_NULL ) {
        return;
    }

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    if ( root >= size ) {
        Rprintf("Invalid root rank specified in reduce, aborting\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    /* Rotate ranks so that root is at rank zero, this makes the logic
     * simpler. */
    MPI_Comm_split(comm, 0, (rank - root + size) % size, &new);
    MPI_Comm_dup(new, &comm);
    MPI_Comm_free(&new);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_group(comm, &g);

    PROTECT(ret = in);
    /* We've reached the root of the tree, assign to the output
     * buffer, which should be defined since we're on the root
     * process, and return. */
    if ( size == 1 ) {
        *out = ret;
        UNPROTECT(1);
        return;
    }

    /* Create group containing the ranks that will participate in this
     * round of communication.  If size is even, this is all of them,
     * otherwise it misses the last one off. */
    ranks = calloc((size_t)size, sizeof(int));
    for ( i = 0; i < size - (size % 2); i++ ) {
        ranks[i] = i;
    }
    MPI_Group_incl(g, size - (size % 2), ranks, &doing_comms);
    /* Create group containing the ranks not participating in the
     * current communication round. */
    MPI_Group_difference(g, doing_comms, &sitting_out);

    MPI_Comm_create(comm, doing_comms, &new);
    /* If we're in the comms round */
    if ( new != MPI_COMM_NULL ) {
        MPI_Comm_rank(new, &rank);
        /* Even rank receives, odd rank sends. */
        if ( rank % 2 == 0 ) {
            neighbour = rank + 1;
            MPI_Recv(&nbytes, 1, MPI_INT, neighbour, 0, new, &status);
            PROTECT(tmp = allocVector(RAWSXP, nbytes));
            MPI_Recv(RAW(tmp), nbytes, MPI_BYTE, neighbour, 0, new, &status);
            ret = (*combine_fn)(in, unserialize_form(tmp));
            UNPROTECT(1);
        } else {
            neighbour = rank - 1; 
            PROTECT(tmp = serialize_form(in));
            nbytes = length(tmp);
            MPI_Isend(&nbytes, 1, MPI_INT, neighbour, 0, new, &request);
            MPI_Isend(RAW(tmp), nbytes, MPI_BYTE, neighbour, 0, new, &request);
            MPI_Wait(&request, &status);
            UNPROTECT(1);
        }
    }

    /* Now create a group containing all the ranks that will
     * participate in the next comms round.  All the even ones from
     * this round... */
    for ( i = 0; i < size/2; i++ ) {
        ranks[i] = i * 2;
    }

    MPI_Group_free(&doing_comms);
    MPI_Group_incl(g, size/2, ranks, &doing_comms);
    MPI_Group_free(&g);
    free(ranks);

    /* ...plus the one that was left out in this round. */
    MPI_Group_union(doing_comms, sitting_out, &g);
    MPI_Group_free(&doing_comms);
    MPI_Group_free(&sitting_out);

    if ( new != MPI_COMM_NULL ) {
        MPI_Comm_free(&new);
    }
    /* Create communicator */
    MPI_Comm_create(comm, g, &new);
    MPI_Group_free(&g);
    MPI_Comm_free(&comm);
    /* Do next level of reduction.  Now the root is zero, because we
     * rotated ranks in the top-level call. */
    reduce_combine(ret, out, combine_fn, 0, new);

    if ( new != MPI_COMM_NULL ) {
        MPI_Comm_free(&new);
    }

    UNPROTECT(1);
    return;
}

/* Workhorse function.  This is where we jump into as a slave. */
int random_forest_driver(int n, ...)
{
    int size;
    int rank;
    MPI_Comm comm;
    va_list ap;
    SEXP args;               /* args to randomForest */
    SEXP tmp;                /* serialized data */
    SEXP *ret;               /* result passed back (on root) */
    int length;
#ifdef _PROF
    double t;
#endif
    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if ( 0 == rank ) {
        va_start(ap, n);
        args = va_arg(ap, SEXP);
        ret = va_arg(ap, SEXP *);
        va_end(ap);

#ifdef _PROF
        t = MPI_Wtime();
#endif
        /* Throw away name of function call */
        PROTECT(args = CDR(args));
        /* Serialize input call.  The data come through here too. */
        PROTECT(tmp = serialize_form(args));

        length = length(tmp);

        /* Send this serialized chunk to everyone else. */
        MPI_Bcast(&length, 1, MPI_INT, 0, comm);
        MPI_Bcast(RAW(tmp), length, MPI_BYTE, 0, comm);
        UNPROTECT(1);           /* tmp */
#ifdef _PROF
        Rprintf("[0: setup-and-send] %g (%d bytes)\n", MPI_Wtime() - t, length);
        fflush(stdout);
#endif
        /* Do some work on master too */
        PROTECT(tmp = serial_randomForest(args, rank, size));
        if ( NULL == tmp ) {
            Rprintf("Can't call prandomForest with more processes than trees, returning NULL\n");
            *ret = R_NilValue;
            UNPROTECT(2);
            return 0;
        }

#ifdef _PROF
        t = MPI_Wtime();
#endif
        /* Get the results back */
        reduce_combine(tmp, ret, &combine_forests, 0, comm);
#ifdef _PROF
        Rprintf("[%d: combine-forests] %g\n", rank, MPI_Wtime() - t);
        fflush(stdout);
#endif
        UNPROTECT(2);           /* tmp and args */
    } else {
#ifdef _PROF
        t = MPI_Wtime();
#endif
        /* Make sure randomForest library is loaded on slave */
        require_package("randomForest");
        /* Get data */
        MPI_Bcast(&length, 1, MPI_INT, 0, comm);
        PROTECT(tmp = allocVector(RAWSXP, length));
        MPI_Bcast(RAW(tmp), length, MPI_BYTE, 0, comm);
        PROTECT(args = unserialize_form(tmp));
        UNPROTECT_PTR(tmp);
#ifdef _PROF
        Rprintf("[%d: setup-and-recv] %g (%d bytes)\n", rank, MPI_Wtime() - t, length);
        fflush(stdout);
#endif
        /* Do our bit of work */
        PROTECT(tmp = serial_randomForest(args, rank, size));
        if ( NULL == tmp ) {
            UNPROTECT(2);
            return 0;
        }

#ifdef _PROF
        t = MPI_Wtime();
#endif
        reduce_combine(tmp, NULL, &combine_forests, 0, comm);
#ifdef _PROF
        Rprintf("[%d: combine-forests] %g\n", rank, MPI_Wtime() - t);
        fflush(stdout);
#endif
        UNPROTECT(2);           /* tmp and args */
    }
    return 0;
}
