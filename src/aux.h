#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "solve.h"

#define LOG(msg) \
{ if (IS_ROOT_PROC) printf(msg); }

double f(
	const int k,
	const int n,
	int i,
	int j
);

int readFile(
	MpiData mpi,
	Block block,
	const char *filename
);

int init(
	MpiData mpi,
	Block block,
	const int k,
	const char *filename
);

void print(
	MpiData mpi,
	double **mem,
	double *buf,
	const int l,
	const int m
);

double calcDiscr(
	MpiData mpi,
	Block block
);

double calcError(
	Block block
);
