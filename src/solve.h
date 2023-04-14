#pragma once
#include <mpi.h>
#include <math.h>
#include <string.h>

#define BUF		(block.matrix[block.rows])
#define RHS(i)	(block.matrix[i][block.n])
#define OFFS(i)	(i / mpi.size + \
(mpi.rank < i % mpi.size))

#define IS_ROOT_PROC (mpi.rank == mpi.root)
#define IS_MAIN_PROC (mpi.rank == globalMain.proc)
#define IS_CURR_PROC (mpi.rank == i % mpi.size)

#define SWAP_ROWS(A, B) \
{ double *tmp; tmp = A; A = B; B = tmp; }

typedef struct MpiData {
	int rank, size, root;
	MPI_Status status;
} MpiData;

typedef struct Block {
	double **matrix;
	double *solution;
	int rows, n;
} Block;

typedef struct Element {
	double value;
	int proc;
} Element;

int solve(
	MpiData mpi,
	Block block
);
