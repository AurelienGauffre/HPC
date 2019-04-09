/* Standard headers */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Math */
#include <math.h>

/* Time */
#include <time.h>

/* Simulation */
#include "kernel.h"
#include "common.h"
#include "functions.h"

#define INDEX2D(i,j) (i + problemInfo.nx*j)

int main(int argc, char* argv[]) {
  /*  Variables */
  BOOL analytical = TRUE;
  clock_t startCPUTime = 0, endCPUTime = 0, totalCPUTime = 0;

  /* Problem Definition */
  ProblemInfo problemInfo;
  if (!setupProblem(argc, argv, &problemInfo)) return 0;
  problemInfo.source = f1;

  /* Memory Allocation */
  double *u   = createVector(problemInfo.nx * problemInfo.ny);
  double *f   = createVector(problemInfo.nx * problemInfo.ny);
  double *f_n = createVector(problemInfo.nx * problemInfo.ny);

  /* Kernel Parameters */
  SolverParameters param;

  param.tol = 1.e-10;
  param.it_max= 20000;
  param.io_interval = 100;

  /* Simulation */
  int iter = 0;

  startCPUTime = clock();
  /* Problem Initialisation */
  initializeField(u ,f , problemInfo);
  while(iter++*problemInfo.dt <= problemInfo.T) {
    printf("\n ================== Time %f ================= \n\n",iter*problemInfo.dt);
    /* Apply time scheme : update f */
    for ( int j = 1; j < problemInfo.ny-1; j++ )
    {
      for ( int i = 1; i < problemInfo.nx-1; i++ )
      {
        f_n[INDEX2D(i,j)] = problemInfo.dt * f[INDEX2D(i,j)] + u[INDEX2D(i,j)];
      }
    }

    /* Solve linear system */
    jacobi(u, f_n, problemInfo , param);

    /* Solution validation */
    if (analytical == TRUE & iter%param.io_interval == 1)
    {
      computeError( u, problemInfo);
    }
  }
  endCPUTime = clock();

  /* Compute effective time */
  totalCPUTime= endCPUTime-startCPUTime;
  printf("Grid %d x %d  (dx=%12.5e,dy=%12.5e) et dt %12.5e \n",
      problemInfo.nx,problemInfo.ny,problemInfo.dx,problemInfo.dy,problemInfo.dt);
  printf("Execution Time %lf s\n",(double) totalCPUTime/CLOCKS_PER_SEC);

  /* Save solution */
  sauvegarde(u, problemInfo,"sol.dat");

  /* Free memory */
  free(u);
  free(f);
  free(f_n);

  return 0;
}

#undef INDEX2D
