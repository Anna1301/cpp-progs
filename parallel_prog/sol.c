#include <math.h>
#include <mpi.h>

#include "sol.h"

void Solution(int n, double *a, double *b, double *x, int rank, int p) {
	int rows;

	if (rank >= n % p) {
        rows = n / p;
    } else {
	    rows = n / p + 1;
	}

	for (int i = 0; i < n; i++) {
		if (rank == i % p) {
			double buf = 1.0 / a[i/p * n + i];
			for (int j = i; j < n; j++) {
                a[(i / p) * n + j] *= buf;
            }
			b[i / p] *= buf;

			for (int j = i; j < n; j++) {
                x[j] = a[i / p * n + j];
            }
			x[n] = b[i / p];

			MPI_Bcast(x, n + 1, MPI_DOUBLE, i%p, MPI_COMM_WORLD);
			for (int j = (i / p) + 1; j < rows; j++) {
				buf = a[j * n + i];
				for (int k = i; k < n; k++) {
                    a[j * n + k] -= buf * a[(i / p) * n + k];
                }
				b[j] -= buf * b[i/p];
			}
			for (int j = 0; j < i / p; j++) {
				buf = a[j * n + i];
				for (int k = i; k < n; k++) {
                    a[j * n + k] -= buf * a[(i / p) * n + k];
                }
				b[j] -= buf * b[i / p];
			}
		} else {
			MPI_Bcast(x, n + 1, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
			for (int j = 0; j < rows; j++) {
				double buf = a[j * n + i];
				for (int k = i; k < n; k++) {
                    a[j * n + k] -= buf * x[k];
                }
				b[j] -= buf * x[n];
			}
		}
	}
}
