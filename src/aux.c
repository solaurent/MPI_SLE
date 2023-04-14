#include "aux.h"

double f(
    const int k,
    const int n,
    int i,
    int j
) {
    ++i; ++j;

    switch (k) {
    case 1: return (double)n - fmax(i, j) + 1.0;
    case 2: return fmax((double)i, (double)j);
    case 3: return fabs((double)i - (double)j);
    case 4: return 1.0 / (double)(i + j - 1);
    default: return 0.0;
    }
}

int readFile(
    MpiData mpi,
    Block block,
    const char *filename
) {
    FILE *file;
    int i, j;

    if (!IS_ROOT_PROC) {
        for (i = 0; i < block.rows; ++i)
            MPI_Recv(
                block.matrix[i],
                block.n,
                MPI_DOUBLE,
                mpi.root,
                MPI_ANY_TAG,
                MPI_COMM_WORLD,
                &mpi.status
            );

        return 0;
    }

    if (!(file = fopen(filename, "r"))) return 1;

    for (i = 0; i < block.n; ++i) {
        for (j = 0; j < block.n; ++j)
            if (fscanf(file, "%lf", &BUF[j]) != 1)
                return 2;

        if (IS_CURR_PROC) SWAP_ROWS(BUF, block.matrix[i / mpi.size])
        else
            MPI_Send(
                BUF,
                block.n,
                MPI_DOUBLE,
                i % mpi.size,
                0,
                MPI_COMM_WORLD
            );
    }

    fclose(file);
    return 0;
}

int init(
    MpiData mpi,
    Block block,
    const int k,
    const char *filename
) {
    int err, i, j;

    if (k)
        for (i = 0; i < block.rows; ++i)
            for (j = 0; j < block.n; ++j)
                block.matrix[i][j] = f(k, block.n, mpi.rank + mpi.size * i, j);
    else
        if ((err = readFile(mpi, block, filename)))
            return err;

    for (i = 0; i < block.rows; ++i) {
        RHS(i) = 0.0;
        for (j = 0; j < block.n; j += 2)
            RHS(i) += block.matrix[i][j];
    }

    return 0;
}

void print(
    MpiData mpi,
    double **mem,
    double *buf,
    const int l,
    const int m
) {
    int i, j, b = l < m ? l : m;
    double *ptr;

    if (!IS_ROOT_PROC) {
        for (i = 0; i < OFFS(b); ++i) {
            MPI_Send(
                mem[i],
                m,
                MPI_DOUBLE,
                mpi.root,
                0,
                MPI_COMM_WORLD
            );
        }

        return;
    }

    for (i = 0; i < b; ++i) {
        if (IS_CURR_PROC) ptr = mem[i / mpi.size];
        else {
            MPI_Recv(
                buf,
                m,
                MPI_DOUBLE,
                i % mpi.size,
                MPI_ANY_TAG,
                MPI_COMM_WORLD,
                &mpi.status
            );

            ptr = buf;
        }

        for (j = 0; j < m; ++j)
            printf("%16.8e ", ptr[j]);
        printf("\n");
    }
}

double calcDiscr(
    MpiData mpi,
    Block block
) {
    int i, j;
    double localNorms[2] = { 0 }, globalNorms[2] = { 0 };

    for (i = 0; i < block.rows; ++i) {

        localNorms[1] += RHS(i) * RHS(i);

        for (j = 0; j < block.n; ++j)
            RHS(i) -= block.matrix[i][j] * block.solution[j];

        localNorms[0] += RHS(i) * RHS(i);
    }

    MPI_Reduce(
        localNorms,
        globalNorms,
        2,
        MPI_DOUBLE,
        MPI_SUM,
        mpi.root,
        MPI_COMM_WORLD
    );

    return IS_ROOT_PROC ? sqrt(globalNorms[0]) / sqrt(globalNorms[1]) : 0.0;
}

double calcError(
    Block block
) {
    int i;
    double tmp, sum = 0.0;

    for (i = 0; i < block.n; ++i) {
        tmp = block.solution[i] - (double)(i % 2 == 0);
        sum += tmp * tmp;
    }

    return sqrt(sum);
}
