# MPI_SLE
Gaussian elimination algorithm parallelized using MPI standard

## Command line arguments:

### n - matrix size
A size of a matrix to read or create. A positive integer.

### printed submatrix size
A bound for the number of rows and columns to be printed (works both for the initial matrix and the solution). Ranges from -1 up to the matrix size.
* -1    the full solution with no auxillary data
* 0     only auxillary data
* 1-n   the first n by n submatrix of the initial matrix, the first n elements of the solution and some auxillary data

### initialization method
An integer in the range 0-4.
* 0     reading a matrix from a file, requires a filename (the last parameter)
* 1..4  different test formulas for filling the matrix, a filename is not required

### filename (optional)
A name of a file to read a matrix from

## run examples
mpirun -np 4 ./run 1000 7 3
mpirun -np 3 ./run 50 5 0 data.txt
