#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <omp.h>

double f(double x);
double g(double x);
void timestamp();

int main()
{
  // Number of discretization point in space
  int m = 100;

  // Number of time step
  int n = 100;

  double alpha;
  double c = 0.25;
  std::ofstream command_unit;
  std::ofstream outfile;

  // Time interval
  double t;
  double t1 = 0.0;
  double t2 = 3.0;

  // Space interval
  double x;
  double x1 = 0.0;
  double x2 = 1.0;

  // Discretization steps
  double dt;
  double dx;

  dx = (x2 - x1) / (double)n;
  dt = (t2 - t1) / (double)m;
  alpha = pow(c * dt / dx, 2);

  // Define array holding the amplitude of the string
  double *u = new double[n + 1];
  double *u_new = new double[n + 1];
  double *u_old = new double[n + 1];
  //  Write data file.
  outfile.open("string_data.dat");

#pragma omp parallel
  {
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    // Define the size allocated to each threads
    int nloc = n / nthreads;
    // Need to take into account overflow if mod(n, nthreads) is not 0
    if (tid < n - nloc * nthreads)
      nloc += 1;

    // Need to consider storage of boundary conditions on each threads
    if (tid != 0 || tid != nthreads - 1)
      nloc += 2;
    else
      nloc += 1;

    // Each thread initialize its own part
    // Initial condition
    if (tid == 0)
      u_old[0] = 0.0;
#pragma omp for
    for (int j = 1; j < n; j++)
    {
      x = j * dx;
      u_old[j] = f(x);
    }
    if (tid != nthreads - 1)
      u_old[n] = 0.0;

    // Bootstrap for the first time step
    if (tid == 0)
      u[0] = 0.0;
#pragma omp for
    for (int j = 1; j < n; j++)
    {
      x = j * dx;
      u[j] = (alpha / 2.0) * (u_old[j - 1] + u_old[j + 1]) + (1.0 - alpha) * u_old[j] + dt * g(x);
    }
    if (tid != nthreads - 1)
      u[n] = 0.0;

    // Loop in time
    for (int i = 2; i <= m; i++)
    {
      if (tid == 0)
        u_new[0] = 0.0;
#pragma omp for
      for (int j = 1; j < n; j++)
      {
        u_new[j] =
            alpha * (u[j - 1] + u[j + 1]) + 2.0 * (1.0 - alpha) * u[j] - u_old[j];
      }
      if (tid != nthreads - 1)
        u_new[n] = 0.0;
#pragma omp master
      {
        t += dt;
        x = 0.;
        for (int j = 0; j <= n; j++)
        {
          outfile << "  " << t << "  " << x
                  << "  " << u_new[j] << "\n";
          x += dx;
        }
        outfile << "\n";
        std::swap(u, u_old);
        std::swap(u, u_new);
      }
    }
#pragma omp barrier
  }
  outfile.close();

  return 0;
}

// Initial condition
double f(double x)
{
  double value;

  if (0.25 <= x && x <= 0.50)
  {
    value = (x - 0.25) * (0.50 - x);
  }
  else
  {
    value = 0.0;
  }

  return value;
}

// Derivative of the initial condition
double g(double x)
{
  double value;

  value = 0.0;
  return value;
}
