#include <stdlib.h>
#include <stdio.h>

#include <math.h>

// Check for input arguments
#include <unistd.h>

#include "common.h"

/* Access macro */
#define INDEX2D(i,j)    ((i)+(j) * nx)

/* Create vector */
double* createVector(int size) {
  double* val;
  val = (double*) malloc((unsigned) size*sizeof(double));
  return val;
}

/* Define parameters of equation to be solved and the discretisation parameters */
bool setupProblem(int argc, char *argv[], ProblemInfo *problemInfo) {
  int c;
  /* Set default values */
  problemInfo->T = 0.1;

  /* Physical parameters */
  problemInfo->eps = 100;
  problemInfo->a = 0.;

  /* Discretisation parameters */
  problemInfo->nx = 32;
  problemInfo->ny = 32;

  problemInfo->CFL = problemInfo->nx*problemInfo->nx * problemInfo->nx *problemInfo->eps;

  /* Parse command line arguments */
  while ((c = getopt(argc , argv , "?hm:n:t:e:a:c:")) != -1 )
  {
    switch (c){
      case 'm' :
        problemInfo->nx = atoi(optarg);
        if (problemInfo->nx <= 0 )
        {
          printf("Illegal parameter\n");
          usage();
        };
        break;
      case 'n' :
        problemInfo->ny = atoi(optarg);
        if (problemInfo->ny <= 0 )
        {
          printf("Illegal parameter\n");
          usage();
        };
        break;
      case 't' :
        problemInfo->T = atof(optarg);
        if (problemInfo->T<= 0 )
        {
          printf("Illegal parameter\n");
          usage();
        };
        break;
      case 'e' :
        problemInfo->eps = atof(optarg);
        if (problemInfo->eps<= 0 )
        {
          printf("Illegal parameter\n");
          usage();
        };
        break;
      case 'a' :
        problemInfo->a = atof(optarg);
        break;
      case 'c' :
        problemInfo->CFL = atof(optarg);
        if (problemInfo->CFL<= 0 )
        {
          printf("Illegal parameter\n");
          usage();
        };
        break;
      case '?' :
        usage();
      default :
        usage();
    }
  }

  /* Initialise discretisation parameters */
  problemInfo->dx =  1.0/problemInfo->nx;
  problemInfo->dy =  1.0/problemInfo->ny;

  problemInfo->dx2 = problemInfo->dx*problemInfo->dx;
  problemInfo->dy2 = problemInfo->dy*problemInfo->dy;

	problemInfo->dt = problemInfo->CFL;

	return true;
}

/* */
void usage(void) {
  printf("Usage: ./chaleur [-a nu] [-c dt] [-e eps] [-m nx] [-n ny] [-t nsecs] \n");
  printf("\t -h -? This message \n");
	printf("\t -a (float) value of convection coefficient\n");
	printf("\t -c (float) time step\n");
	printf("\t -e (float) valie of diffusion coefficient\n");
  printf("\t -m (int) number of points in x direction\n");
  printf("\t -n (int) number of points in y direction\n");
  printf("\t -t (double) physical time of simulation\n");
  exit (8);
}

/* Save data to ASCII file */
void sauvegarde(double u[], ProblemInfo problem, char* fileName) {
	int nx = problem.nx;
	int ny = problem.ny;
	double deltaX =(double) 1.0/problem.nx;
	double deltaY =(double) 1.0/problem.ny;
	FILE* dataFile;
	dataFile = fopen(fileName,"w");
	for (int ix = 0; ix < nx; ix++)
  {
		for (int jy = 0; jy < ny;jy++)
    {
			fprintf(dataFile,"%lf    %lf     %lf\n",ix*deltaX,jy*deltaY,u[INDEX2D(ix,jy)]);
		}
		fprintf(dataFile,"\n");
	}
	fclose(dataFile);
}

/* Define value at the boundary with Dirichlet boundary conditions */
void initializeField(double u[] ,double f[] , ProblemInfo problemInfo) {
  int    i, j;
  double x, y;

	int nx = problemInfo.nx;
  /* Initialize at the domain boundary */
  j = 0;
  y = j * problemInfo.dy;

  for ( i = 0; i < problemInfo.nx; i++ ) {
    x = i * problemInfo.dx;
    u[INDEX2D(i,j)] = problemInfo.source ( x, y, 0 );
    f[INDEX2D(i,j)] =  -problemInfo.eps*(problemInfo.source( x, y, 2 ) + problemInfo.source( x, y ,1));
  }

  j = problemInfo.ny - 1;
  y = j * problemInfo.dy;

  for ( i = 0; i < problemInfo.nx; i++ ) {
    x = i * problemInfo.dx;
    u[INDEX2D(i,j)] = problemInfo.source ( x, y, 0 );
    f[INDEX2D(i,j)] =  -problemInfo.eps*(problemInfo.source( x, y, 2 ) + problemInfo.source( x, y ,1));
  }

  i = 0;
  x = i * problemInfo.dx;

  for ( j = 0; j < problemInfo.ny; j++ ) {
    y =  j * problemInfo.dy;
    u[INDEX2D(i,j)] = problemInfo.source ( x, y, 0 );
    f[INDEX2D(i,j)] =  -problemInfo.eps*(problemInfo.source( x, y, 2 ) + problemInfo.source( x, y ,1));
  }

  i = problemInfo.nx - 1;
  x = i * problemInfo.dx;

  for ( j = 0; j < problemInfo.ny; j++ ) {
    y = j * problemInfo.dy;
    u[INDEX2D(i,j)] = problemInfo.source ( x, y, 0 );
    f[INDEX2D(i,j)] =  -problemInfo.eps*(problemInfo.source( x, y, 2 ) + problemInfo.source( x, y ,1));
  }

  /* Initialize inside the domain */
  for ( j = 1; j < problemInfo.ny - 1; j++ ){
    y = j * problemInfo.dy;
    for ( i = 1; i < problemInfo.nx - 1; i++ ) {
      x = i * problemInfo.dx;
      u[INDEX2D(i,j)] = 0.;
      f[INDEX2D(i,j)] = -problemInfo.eps*(problemInfo.source( x, y, 2 ) + problemInfo.source( x, y ,1));
    }
  }

  return ;
}

void computeError(double u[], ProblemInfo problemInfo) {
    /* Error variables */
    double  error_max, error_l2;

    /* Coordinates for source */
    double  x, y, term;

    int nx = problemInfo.nx;

    /* Discretisation */
    double dx = problemInfo.dx;
    double dy = problemInfo.dy;

    /* Error initialisation */
    error_max   = 0.0;
    error_l2    = 0.0;

    for (int j = 0; j < problemInfo.ny; j++ ) {
      y = j * dy;
      for (int i = 0; i < problemInfo.nx; i++ ) {
            x =  i * dx;
            term   =  u[INDEX2D(i,j)] - problemInfo.source( x, y, 0);
            error_l2  = error_l2 + term*term;
            error_max = fmax(fabs(term), error_max);
        }
      }

    error_l2 = sqrt(dx*dy*error_l2);

    printf( "Linfinite Error : %12.5e \n", error_max );
    printf( "L2 Error        : %12.5e \n", error_l2 );

    return;
}

#undef INDEX2D
