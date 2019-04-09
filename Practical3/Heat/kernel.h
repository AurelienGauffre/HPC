#ifndef KERNEL_H
#define KERNEL_H

#include "common.h"

/* Computation parameters */
typedef struct {
  int it_max;
  int io_interval;
  double tol;
} SolverParameters;

/* Jacobi solver */
void jacobi(double u[], double f[], ProblemInfo problemInfo, SolverParameters param);

#endif
