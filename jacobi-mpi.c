#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "util.h"

double residual(const double *u, const double *f, long int n)
{
	if (n == 1)
	{
		return abs(2*(n+1)*(n+1)*u[0] -f[0]);
	}
	double rss=0;
	rss += pow((n+1)*(n+1)*(2*u[0]-u[1]) - f[0],2);
	long int i = 1;
	while ( i < n-1)
	{
		rss += pow((n+1)*(n+1)*(2*u[i]-u[i+1]-u[i-1]) - f[i],2);
		i++;
	}
	rss += pow((n+1)*(n+1)*(2*u[n-1]-u[i-2]) - f[n-1],2);
	return sqrt(rss);
}

void jacobi(double *u, const double *f, long int n)
{
	if (n == 1)
	{
		u[0] = f[0]/(2*(n+1)*(n+1));
		return;
	}
	double u_old = 0;
	long int i;
	for (i = 0; i < n-1; i++)
	{
		double tmp = f[i]/(2*(n+1)*(n+1))+0.5*u_old+0.5*u[i+1];
		u_old = u[i];
		u[i] = tmp;
	}
	u[n-1] = f[i]/(2*(n+1)*(n+1))+0.5*u_old;
	return;
}

int main (int argc, char **argv)
{
	if (argc != 3)
	{
		fprintf(stderr,"Function needs vector length and number of iterations as input arguments!\n");
		abort();
	}
	
	long int n = atoi(argv[1]), k = atoi(argv[2]);
	
	if (n <=0)
	{
		fprintf(stderr,"Vector length must be positive!\n");
		abort();
	}
	else if (k <= 0)
	{
		fprintf(stderr, "Number of iteration must be positive!\n");
		abort();
	}

	int rank, size;
	MPI_Status status;

	timestamp_type time0, time1;
	get_timestamp(&time0);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (n % size)
	{
		if (rank == 0)
			fprintf(stderr, "Vector length is not a multiply of process number!");
		abort();
	}

	int after = rank + 1, before = rank - 1;
	int tag = 99;
	
	double *u, *f;
	int length = n / size + 2;	//vector length in each process.

	u = (double *) malloc(sizeof(double)*length);
	f = (double *) malloc(sizeof(double)*length);
	
	int i;
	for ( i = 0; i < length; i++)
	{
		u[i] = 0;
		f[i] = 1;
	}
	
	for (i = 0; i < k; i++)
	{
		jacobi(u, f, length);

		if (rank == 0)
		{
			MPI_Send(u + length - 2, 1, MPI_DOUBLE, after, tag, MPI_COMM_WORLD);
			MPI_Recv(u + length - 1, 1, MPI_DOUBLE, after, tag, MPI_COMM_WORLD, &status);
			u[0] = 0;
		}
		else if (rank == size - 1)
		{
			MPI_Send(u + 1, 1, MPI_DOUBLE, before, tag, MPI_COMM_WORLD);
			MPI_Recv(u , 1, MPI_DOUBLE, before, tag, MPI_COMM_WORLD, &status);
			u[length-1] = 0;
		}
		else
		{
			MPI_Send(u + 1, 1, MPI_DOUBLE, before, tag, MPI_COMM_WORLD);
			MPI_Send(u + length - 2, 1, MPI_DOUBLE, after, tag, MPI_COMM_WORLD);
			MPI_Recv(u, 1, MPI_DOUBLE, before, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(u + length - 1, 1, MPI_DOUBLE, after, tag, MPI_COMM_WORLD, &status);
		}

	}
	
	free(u);
	free(f);

        get_timestamp(&time1);

	if (rank == 0) {
		double elapsed = timestamp_diff_in_seconds(time0, time1);
		printf("Time elapsed: %f seconds.\n", elapsed);
	}

	MPI_Finalize();

	return 0;
	
}
