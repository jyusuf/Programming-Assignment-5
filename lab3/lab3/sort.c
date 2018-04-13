#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Prints a vector without decimal places...
 */
void print_vector(int n, double vector[]) {
	int i;
	printf("[%.0f", vector[0]);
    for (i = 1; i < n; i++) {
  		printf(", %.0f", vector[i]);
    }
	printf("]\n");
}

/**
 * Just checks sequentially if everything is in ascending order.
 *
 * Note: we don't care about stability in this sort since we
 * have no data attached to the double value, see
 * http://en.wikipedia.org/wiki/Sorting_algorithm#Stability
 * for more detail.
 */
void test_correctness(int n, double v[]) {
	int i;
	for (i = 1; i < n; i++) {
		if (v[i] < v[i-1]) {
			printf("Correctness test found error at %d: %.4f is not < %.4f but appears before it\n",
					i, v[i-1], v[i]);
		}
	}
}

/**
 * Initialize random vector.
 *
 * You may not parallelize this (even though it could be done).
 */
void init_random_vector(int n, double v[]) {
    int i, j;
	for (i = 0; i < n; i++) {
		v[i] = rand() % n;
	}
}

double *R;

/**
 * Merges two arrays, left and right, and leaves result in R
 */
void merge(double *left_array, double *right_array, int leftCount, int rightCount) {
	int i,j,k;

	// i - to mark the index of left aubarray (left_array)
	// j - to mark the index of right sub-raay (right_array)
	// k - to mark the index of merged subarray (R)
	i = 0; j = 0; k =0;

	while (i < leftCount && j < rightCount) {
		if(left_array[i]  < right_array[j])
			R[k++] = left_array[i++];
		else
			R[k++] = right_array[j++];
	}
	while (i < leftCount)
		R[k++] = left_array[i++];
	while (j < rightCount)
		R[k++] = right_array[j++];
}

/**
 * Recursively merges an array of n values using group sizes of s.
 * For example, given an array of 128 values and starting s value of 16
 * will result in 8 groups of 16 merging into 4 groups of 32, then recursively
 * calling merge_all which merges them into 2 groups of 64, then once more
 * recursively into 1 group of 128.
 */
void merge_all(double *v, int n, int s) {
  if (s < n) {
    int i;
    for (i = 0; i < n; i += 2*s) {
      merge(v+i, v+i+s, s, s);
      // result is in R starting at index 0
      memcpy(v+i, R, 2*s*sizeof(double));
    }
    merge_all(v, n, 2*s);
  }
}

void merge_sort(int n, double *v) {
  merge_all(v, n, 1);
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		// () means optional
		printf("usage: %s n (seed)\n", argv[0]);
		return 0;
	}
	int n = atoi(argv[1]);
	int seed = 0;
	if (argc > 2) {
		seed = atoi(argv[2]);
	}
  srand(seed);

	int mpi_p, mpi_rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_p);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // init the temporary array used later for merge
  R = malloc(sizeof(double)*n);

	// all MPI processes will allocate a vector of length n, but only
	// process 0 will initialize the values.  You must distribute the
	// values to processes.  You will also likely need to allocate
	// additional memory in your processes, make sure to clean it up
	// by adding a correct free at the end of each process.
  double *v = malloc(sizeof(double)*n);
  if (mpi_rank == 0) {
  	init_random_vector(n, v);
	}

	double start = MPI_Wtime();

	// do all the work ourselves! (you should make a better algorithm here!)
  if (mpi_rank == 0) {
	   merge_sort(n, v);
	}

	double end = MPI_Wtime();

	if (mpi_rank == 0) {
		printf("Total time to solve with %d MPI Processes was %.6f\n",
				mpi_p, (end-start));

		test_correctness(n, v);
	}

	MPI_Finalize();

	free(v);
  free(R);
  return 0;
}
