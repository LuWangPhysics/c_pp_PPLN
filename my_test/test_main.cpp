#include <stdio.h>
#include "mpi.h"
#include <omp.h>
#include <iostream>
#include <sstream>
int main(int argc, char *argv[]) {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
omp_set_nested(1);
omp_set_dynamic(0);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  #pragma omp parallel default(shared) num_threads(16)
  { int np = omp_get_thread_num(); //get the first level thread number
 
     #pragma omp parallel num_threads(2)
    {   
     int iam = omp_get_thread_num();
    printf("Hello from thread %d out of %d from process %d out of %d on %s\n",
           iam, np, rank, numprocs, processor_name);
     }
  }

  MPI_Finalize();
}