#include "solve.h"

int equal(const float a, const float b);

int equal(const float a, const float b) {
	return fabs(a - b) < 1e-16;
}

int solve(
	MpiData mpi,
	Block block
) {
	Element localMain = { 0 }, globalMain = { 0 };
	int mainRow, i, j, k;

	for (i = 0; i < block.n; ++i) {

		// searching for the main element

		for (j = mainRow = OFFS(i); j < block.rows; ++j)
			if (fabs(block.matrix[j][i]) > fabs(block.matrix[mainRow][i]))
				mainRow = j;

		localMain.value = (OFFS(i) < block.rows) ? fabs(block.matrix[mainRow][i]) : 0.0;
		localMain.proc = mpi.rank;

		MPI_Allreduce(
			&localMain,
			&globalMain,
			1,
			MPI_DOUBLE_INT,
			MPI_MAXLOC,
			MPI_COMM_WORLD
		);

		if (equal(globalMain.value, 0.0)) return 1;

		// normalization

		if (IS_MAIN_PROC) {
			SWAP_ROWS(BUF, block.matrix[mainRow]);
			for (j = block.n; j >= i; --j)
				BUF[j] /= BUF[i];
		}

		// broadcasting the main row

		MPI_Bcast(
			BUF + i,
			block.n + 1 - i,
			MPI_DOUBLE,
			globalMain.proc,
			MPI_COMM_WORLD
		);

		// triangulation

		for (j = OFFS(i); j < block.rows; ++j)
			for (k = block.n; k >= i; --k)
				block.matrix[j][k] -= BUF[k] * block.matrix[j][i];

		// swapping/exchanging rows

		if (IS_MAIN_PROC && IS_CURR_PROC) {
			SWAP_ROWS(BUF, block.matrix[mainRow]);
			SWAP_ROWS(block.matrix[mainRow], block.matrix[i / mpi.size]);
		}

		if (!IS_MAIN_PROC && IS_CURR_PROC) {
			MPI_Send(
				block.matrix[i / mpi.size] + i,
				block.n + 1 - i,
				MPI_DOUBLE,
				globalMain.proc,
				0,
				MPI_COMM_WORLD
			);
			SWAP_ROWS(BUF, block.matrix[i / mpi.size]);
		}

		if (IS_MAIN_PROC && !IS_CURR_PROC)
			MPI_Recv(
				block.matrix[mainRow] + i,
				block.n + 1 - i,
				MPI_DOUBLE,
				i % mpi.size,
				MPI_ANY_TAG,
				MPI_COMM_WORLD,
				&mpi.status
			);
	}
	
	// diagonalization

	for (i = block.n - 1; i >= 0; --i) {

		if (IS_CURR_PROC) SWAP_ROWS(BUF, block.matrix[i / mpi.size]);

		MPI_Bcast(
			BUF + block.n,
			1,
			MPI_DOUBLE,
			i % mpi.size,
			MPI_COMM_WORLD
		);

		for (j = 0; j < OFFS(i); ++j)
			RHS(j) -= block.matrix[j][i] * BUF[block.n];

		if (IS_CURR_PROC) SWAP_ROWS(BUF, block.matrix[i / mpi.size]);
	}

	// writing a solution fragment to the buffer

	memset(BUF, 0, sizeof(double) * block.n);

	for (i = 0; i < block.rows; ++i)
		BUF[i * mpi.size + mpi.rank] = RHS(i);

	// obtaining the entire solution in every process

	MPI_Allreduce(
		BUF,
		block.solution,
		block.n,
		MPI_DOUBLE,
		MPI_SUM,
		MPI_COMM_WORLD
	);
	
	return 0;
}
