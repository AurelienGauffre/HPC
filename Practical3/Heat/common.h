#ifndef COMMON_H
#define COMMON_H

typedef enum {
  FALSE, TRUE
} BOOL;

/* Problem definition */
typedef struct {
	int nx;
	int ny;
	double T;
	double eps;
	double CFL;

  double dx;
  double dy;

  double dx2;
  double dy2;
  double dt;

  /* Scheme normalisation coefficients */
  double a;
  double b;

  double (*source)(double, double, int);
} ProblemInfo ;


/* Manage arrays */
double* createVector(int arraySize);

void initializeField(double u[] ,double f[] , ProblemInfo problemInfo);
void computeError(double u[], ProblemInfo problemInfo);

/* Input/Output */
bool setupProblem(int argc, char *argv[], ProblemInfo *problemInfo);
void usage(void) __attribute__ ((noreturn));
void sauvegarde(double u[], ProblemInfo problem, char* filename);

#endif
