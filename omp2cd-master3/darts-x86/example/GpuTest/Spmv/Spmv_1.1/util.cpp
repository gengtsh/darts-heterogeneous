#include "util.h"


// ****************************************************************************
// Function initRandomMatrix
//
// Purpose:
//   Assigns random positions to a given number of elements in a square
//   matrix, A.  The function encodes these positions in compressed sparse
//   row format.
//
// Arguments:
//   cols:          array for column indexes of elements (size should be = n)
//   rowDelimiters: array of size dim+1 holding indices to rows of A;
//                  last element is the index one past the last element of A
//   n:             number of nonzero elements in A
//   dim:           number of rows/columns in A
//
// Programmer: Kyle Spafford
// Creation: July 28, 2010
// Returns: nothing
//
// ****************************************************************************
void initRandomMatrix(int *cols, int *rowDelimiters, const int n, const int dim)
{
    int nnzAssigned = 0;

    // Figure out the probability that a nonzero should be assigned to a given
    // spot in the matrix
    double prob = (double)n / ((double)dim * (double)dim);

    // Seed random number generator
    srand48(8675309L);

    // Randomly decide whether entry i,j gets a value, but ensure n values
    // are assigned
    bool fillRemaining = false;
    for (int i = 0; i < dim; i++)
    {
        rowDelimiters[i] = nnzAssigned;
        for (int j = 0; j < dim; j++)
        {
            int numEntriesLeft = (dim * dim) - ((i * dim) + j);
            int needToAssign   = n - nnzAssigned;
            if (numEntriesLeft <= needToAssign) {
                fillRemaining = true;
            }
            if ((nnzAssigned < n && drand48() <= prob) || fillRemaining)
            {
                // Assign (i,j) a value
                cols[nnzAssigned] = j;
                nnzAssigned++;
            }
        }
    }
    // Observe the convention to put the number of non zeroes at the end of the
    // row delimiters array
    rowDelimiters[dim] = n;
    assert(nnzAssigned == n);
}
