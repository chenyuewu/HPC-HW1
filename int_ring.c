#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "util.h"

int main(int argc, char **argv)
{
  if (argc != 2) {
    fprintf(stderr, "Need parameter to specify how often the message is sent!");
    return 1;
  }

  int N = atol(argv[1]);

  if (N <= 0) {
    fprintf(stderr, "Number of loops must be positive!");
    return 2;
  }

  int rank, size;
  MPI_Status status;

  timestamp_type time0, time1;
  get_timestamp(&time0);

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int message;

  if (rank == 0)
    message = 0;

  int destination = (rank + 1) % size;
  int origin = (rank - 1) % size;
  origin = origin < 0? origin + size : origin;
  int tag = 99;

  //printf("process %d has destination %d and origin %d.\n", rank, destination, origin);

  int i = 0;

  for (i = 0; i < N; i++)
  {
    if (rank == 0) {
      MPI_Send(&message, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
      //printf("At loop %d, process %d send %d to process %d.\n", i, rank, message, destination);
      MPI_Recv(&message, 1, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
    }
    else {
      MPI_Recv(&message, 1, MPI_INT, origin, tag, MPI_COMM_WORLD, &status);
      message = message + rank;
      MPI_Send(&message, 1, MPI_INT, destination, tag, MPI_COMM_WORLD);
      //printf("At loop %d, process %d send %d to process %d.\n", i, rank, message, destination);
    }
  }

  if (rank == 0) {
    printf("After %d loops on %d processes, the integer grows to %d.\n", N, size, message);
    get_timestamp(&time1);
    double elapsed = timestamp_diff_in_seconds(time0,time1);
    printf("Total communication: %d.\nTime elapsed: %f seconds.\nLatency: %f seconds.\n", size*N, elapsed, elapsed/size/N);
  }
  
  MPI_Finalize();
  return 0;

}
