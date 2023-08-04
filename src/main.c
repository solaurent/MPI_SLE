#include "aux.h"

int main(int argc, char **argv) {

	MpiData mpi = { 0 };
	Block block = { 0 };

	const int n = atoi(argv[1]);
	const int m = atoi(argv[2]);
	const int k = atoi(argv[3]);
	const char *filename = k ? NULL : argv[4];

	int err, i;
	double discr, error, time;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi.size);

	if (n < 1 || m < -1 || m > n || k < 0 || k > 4) {
		LOG("\nSome parameters are out of range\n\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (n < mpi.size) {
		LOG("\nSome processes have no payload\n\n");
		MPI_Abort(MPI_COMM_WORLD, 2);
	}

	block.rows = OFFS(n);
	block.n = n;

	// (n+1) x (n+1) array
	// the last column - rhs
	// the last row - buffer

	block.matrix = (double **)malloc(
		sizeof(double *) * (block.rows + 1) +
		sizeof(double) * (block.rows + 1) * (block.n + 1)
	);

	block.solution = (double *)malloc(sizeof(double) * block.n);

	if (!block.matrix || !block.solution) {
		printf("\nProcess %d: memory error\n\n", mpi.rank);
		MPI_Abort(MPI_COMM_WORLD, 3);
	}

	// setting pointers

	block.matrix[0] = (double *)(block.matrix + (block.rows + 1));
	for (i = 1; i < block.rows + 1; ++i)
		block.matrix[i] = block.matrix[i - 1] + (block.n + 1);
	
	// initialization

	err = init(mpi, block, k, filename);

	if (err == 1) {
		LOG("\nUnable to open the file\n\n");
		MPI_Abort(MPI_COMM_WORLD, 4);
	}

	if (err == 2) {
		LOG("\nUnexpected symbols in the file\n\n");
		MPI_Abort(MPI_COMM_WORLD, 5);
	}

	// output

	if (m > 0) {
		LOG("Matrix:\n");
		print(mpi, block.matrix, BUF, n, m);
	}

	// solving

	MPI_Barrier(MPI_COMM_WORLD);
	time = MPI_Wtime();
	err = solve(mpi, block);
	time = MPI_Wtime() - time;

	if (err == 1) {
		LOG("\nThe given matrix is singular\n\n");
		MPI_Abort(MPI_COMM_WORLD, 6);
	}

	if (m > 0) {
		LOG("Solution:\n");
		print(mpi, &block.solution, BUF, 1, m);
	}
	if (m == -1) print(mpi, &block.solution, BUF, 1, n);

	// discr and error calculation

	if (m != -1) {
		init(mpi, block, k, filename);
		discr = calcDiscr(mpi, block);
		error = calcError(block);

		if (IS_ROOT_PROC) {
			printf("Discr:\t%16.8e\n", discr);
			printf("Error:\t%16.8e\n", error);
			printf("Time:\t%16.3f\n", time);
		}
	}

	free(block.matrix);
	free(block.solution);
	MPI_Finalize();
	return 0;
}
