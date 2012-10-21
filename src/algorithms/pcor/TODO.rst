=========
PCOR TODO
=========

#. Support executing pcor with only one thread.

    -   Add a new function that does the same as the parallel but without MPI ?
    -   Use the same but pass a flag to skip the MPI calls? Will make the existing implementation
        tricky to manage...probably not a good idea.
    -   What will happen to the final result? pcor creates a 'ff' object to the file created by
        'MPI_File', I'll need to flush them in a file to be consistent with what the parallel
        implementation does.

#. Add the other two algorithms to the pcor.R file and pass the method flag down to C interface.

    -   The UI needs to do some parameter checking for the two new functions.
    -   I could possibly borrow the code from cor.R. 

#. Change the call in the 'kernel.c' file to process an entire line at a time

    -   Current call computes a single coefficient at a time and it's been called
        'rows' number of times. More efficient to minimize the function calls. I already
        have an implementation of the Pearson algorithm to do this.

