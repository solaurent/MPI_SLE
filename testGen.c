#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main(int argc, char **argv) {

	int i, j;
	const char *testName = argv[1];
	const int testSize = atoi(argv[2]);
	double sample;

	FILE *file = fopen(testName, "w");
	srand(time(NULL));

	for (i = 0; i < testSize; ++i) {
		for (j = 0; j < testSize; ++j) {
			sample = (double)rand() / (double)RAND_MAX * 100.0;
			fprintf(file, "%2.2f\t", sample);
		}
		fprintf(file, "\n");
	}

	return 0;
}
