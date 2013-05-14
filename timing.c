/*
 * Copyright (C) 2013 Edward B. Duffy
 * License: http://www.gnu.org/licenses/gpl.html GPL version 2 or higher
 */

#include <stdio.h>
#include <mpi.h>
#include "divvy.h"

extern int quiet;
extern int verbose;

void tic(double *clock)
{
   *clock = MPI_Wtime();
}

void toc(double *clock, const char *message)
{
   int rank;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   if(rank == 0) {
      printf("%-13s %12.2lf sec\n", message, MPI_Wtime() - *clock);
   }
   *clock = MPI_Wtime();
}
