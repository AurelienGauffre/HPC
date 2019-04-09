/* Standard headers */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Math */
#include <math.h>

/* Heat tools*/
#include "common.h"
#include "kernel.h"

/* Access Macro */
#define INDEX2D(i,j)    ((i)+(j) * nx)


void jacobi(double u[], double f[], ProblemInfo problemInfo, SolverParameters param) {
  /* Loop indices */
  int i, j, it, nx;

  /* Stencil coefficients */
  double ax, ay, d;

  /* Compute norm */
  double norm = 0., diff = 0.;

  /* Backup intermediate steps */
  double *u_old;

  /* Initialise coefficients */
  nx = problemInfo.nx;
  d  = problemInfo.eps * problemInfo.dt;
  ax = d /problemInfo.dx2;
  ay = d /problemInfo.dy2;

  u_old = createVector(problemInfo.nx*problemInfo.ny);

  for ( it = 1; it <= param.it_max; it++ ) {
    norm = 0.0;

    /* Duplicate new solution to old one */
    memcpy(u_old, u, (unsigned) problemInfo.nx * (unsigned) problemInfo.ny * sizeof( double ));

    /* Compute the value of the stencil and update it */
    for ( j = 1; j < problemInfo.ny-1; j++ ) {
      for ( i = 1; i < problemInfo.nx-1; i++ ) {
        u[INDEX2D(i,j)] = (f[INDEX2D(i,j)]   +
            ( ax * ( u_old[INDEX2D(i-1,j)] + u_old[INDEX2D(i+1,j)] ) +
              ay * ( u_old[INDEX2D(i,j-1)] + u_old[INDEX2D(i,j+1)] ) ) ) / ( 1 + 2*(ax+ay) );

        diff = fabs(u[INDEX2D(i,j)] - u_old[INDEX2D(i,j)]);

        /* Calcul de la norme 2 */
        norm = norm + diff*diff;
      }
    }

#ifdef VERBOSE
    /* Display evolution of error */
    if (0 == it% param.io_interval)
    {
      printf ( "iteration  %5d norme %12.3e\n", it, norm );
    }
#endif
    /* Check error in the method  */
    if ( norm <= param.tol )
    {
      break;
    }
  }

#ifdef VERBOSE
  /* Display data at the end of the iteration */
  printf ( "== Nombre total d'iteration %d avec l'erreur %12.3e \n", it, norm );
#endif

  /* Free allocated memory */
  free ( u_old );

  return;
}


/* Cleanup macros */
#undef INDEX2D
