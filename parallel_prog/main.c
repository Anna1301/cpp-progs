#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "sol.h"

#define MAX_OUTPUT_SIZE 5

double f(int i, int j) {
    return 1.0 / (i + j + 1.0);
}

void ScanMatrix(int n, double *a, double *b, int rank, int p) {
    int rows;
    if (rank >= n % p) {
        rows = n / p;
    } else {
        rows = n / p + 1;
    }

    for (int i = 0; i < rows; i++) {
        double x = 0.0;
        for (int j = 0; j < n; j++) {
            a[i * n + j] = f(rank + p * i, j);
            if (!(j & 1)) {
                x += a[i * n + j];
            }
        }
        b[i] = x;
    }
}

void PrintMatrix(int n, double *a, double *b, double *x, int rank, int p) {
    MPI_Status status;

    int m = (n < MAX_OUTPUT_SIZE) ? n : MAX_OUTPUT_SIZE;

    for (int i = 0; i < m; i++) {
        if (rank == 0) {
            if (rank == i % p) {
                printf("| ");
                for (int j = 0; j < m; j++) {
                    printf("%10.3g ", a[(i / p) * n + j]);
                }
                printf("|   %10.3g\n", b[i / p]);
            } else {
                MPI_Recv(x, m + 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &status);
                printf("| ");
                for (int j = 0; j < m; j++) {
                    printf("%10.3g ", x[j]);
                }
                printf("|   %10.3g\n", x[m]);
            }
        } else if (rank == i % p) {
            for (int j = 0; j < m; j++) {
                x[j] = a[(i / p) * n + j];
            }
            x[m] = b[i / p];
            MPI_Send(x, m + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
}

void PrintVector(int n, double *b, double *x, int rank, int p) {
    MPI_Status status;

    int m = (n < MAX_OUTPUT_SIZE) ? n : MAX_OUTPUT_SIZE;

    for (int i = 0; i < m; i++) {
        if (rank == 0) {
            if (rank == i % p) {
                printf("%10.3g ", b[i / p]);
            } else {
                MPI_Recv(x, 1, MPI_DOUBLE, i % p, 0, MPI_COMM_WORLD, &status);
                printf("%10.3g ", x[0]);
            }
        } else if (rank == i % p) {
            x[0] = b[i / p];
            MPI_Send(x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char **argv)
{
    int rank, p;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n;
	if (argc == 2) {
	    n = atoi(argv[1]);
	} else {
		if (rank == 0) {
		    printf("Not correct usage.\n");
		}

		MPI_Finalize();
		return 1;
	}

    int rows;
	if (rank >= n % p) {
	    rows = n / p;
	} else {
	    rows = n / p + 1;
	}

    double* a = (double*)malloc(rows * n * sizeof(*a));
	double* b = (double*)malloc(rows * sizeof(*b));
	double* x = (double*)malloc((n + 1)* sizeof(*x));

    int err1 = 0, err2 = 0;
	if (!(a && b && x)) {
	    err1 = 1;
	}

	MPI_Allreduce(&err1, &err2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if (err2) {
		if (rank == 0) {
		    printf("Not enough memory.\n");
		}

		if (a) {
		    free(a);
		}
		if (b) {
		    free(b);
		}
		if (x) {
		    free(x);
		}

		MPI_Finalize();

		return 1;
	}

	ScanMatrix(n, a, b, rank, p);

	if (rank == 0) {
	    printf("Matrix A:\n\n");
	}
	PrintMatrix(n, a, b, x, rank, p);

	MPI_Barrier(MPI_COMM_WORLD);
    double t = MPI_Wtime();

	Solution(n, a, b, x, rank, p);

	MPI_Barrier(MPI_COMM_WORLD);
	t = MPI_Wtime() - t;

	if (rank == 0) {
	    printf("\nSolution:\n");
	}
	PrintVector(n, b, x, rank, p);

	if (rank == 0) {
	    printf("\n\nSolution time = %e\n", t);
	}

	free(a);
	free(b);
	free(x);

	MPI_Finalize();

	return 0;
}
