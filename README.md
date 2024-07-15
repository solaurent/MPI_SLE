# MPI_SLE
Gaussian elimination algorithm parallelized using MPI interface

## Command line arguments:

### n - matrix size
Size of a matrix to read or to create. A positive integer.

### m - printed submatrix size
A bound for a number of rows and columns to be printed (works both for the initial matrix and a solution). Ranges from -1 up to the matrix size.
* -1    the full solution with no auxillary data
* 0     only auxillary data
* 1..n  the first n by n submatrix of the initial matrix, the first n elements of the solution and some auxillary data

### k - initialization method
An integer in range 0-4.
* 0     reading the matrix from a file, requires a filename (which is given as the last arguement)
* 1..4  different test formulas for filling the matrix, the filename is not required

### filename (optional)
A name of the file to read the matrix from

## run examples
* mpirun -np 4 ./run 1000 7 3
* mpirun -np 3 ./run 50 5 0 data.txt
