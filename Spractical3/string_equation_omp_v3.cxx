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

#define SAVEFILE 0
int main()
{
  // Number of time step
  int m = 100;

  // Number of discretization point in space
  int n = 101;

  double c = 0.25;

  // Time interval
  double t;
  double t1 = 0.0;
  double t2 = 3.0;

  // Space interval
  double x1 = 0.0;
  double x2 = 1.0;

  // Discretization steps
  double dt;
  double dx;

  double alpha;
  dx = (x2 - x1) / (double)(n - 1);
  dt = (t2 - t1) / (double)m;
  alpha = pow(c * dt / dx, 2);

  // We need to store boundary conditions for each thread in order
  // to exchange them with their neighbors
  double *u_l, *u_r;

// Start parallel region
#pragma omp parallel shared(u_l, u_r)
  {
#if (SAVEFILE == 1)
    // Each thread will save its own data
    std::ofstream outfile;
#endif
    // Get thread information
    int nthreads = omp_get_num_threads();
    int tid = omp_get_thread_num();

    // Define the size to allocate for each thread
    int nx = n / nthreads;
    int inc = n - nx * nthreads;

    // We need to take into account overflow if modulo(n, nthreads) is not 0
    if (tid < inc)
    {
      nx += 1;
    }

    // Set starting and ending indices for each thread in global
    // computational domain
    int j0 = (n / nthreads) * tid + (tid < inc ? tid : inc);
    int j1 = j0 + nx - 1;

    // Starting and ending incdices for each thread in local computational domain
    int j0l = 1;
    int j1l = nx;
    j1l = j0l + nx - 1;

    // Consider overlapping domains between thread
    int nx_s = nx;
    if (tid != 0)
    {
      nx_s += 1;
    }
    if (tid != nthreads - 1)
    {
      nx_s += 1;
    }

// Allocation for left and right boundary exchange arrays
#pragma omp single
    {
      u_l = new double[nthreads];
      u_r = new double[nthreads];
    }

#pragma omp barrier

    // Initialize boundary values
    u_l[tid] = 0.;
    u_r[tid] = 0.;

    //  Arrays storing the amplitude of the string
    double *u = new double[nx_s];
    double *u_new = new double[nx_s];
    double *u_old = new double[nx_s];

#if (SAVEFILE == 1)
    // Open data file for specfic thread
    outfile.open("string_data_" + std::to_string(tid) + ".dat");
#endif
    // Starting point of thread domain
    double x0 = j0 * dx;
    double x;

#pragma omp barrier
    for (int j = j0l, k = j0; j <= j1l; j++, k++)
    {
      x = k * dx;
      u_old[j] = f(x);
    }

#pragma omp barrier
    // initialize boundary arrays
    if (tid != nthreads - 1)
    {
      u_l[tid + 1] = u_old[j1l];
    }
    if (tid != 0)
    {
      u_r[tid - 1] = u_old[1];
    }
#pragma omp barrier

    // Finalize initialization of boundaries on current thread
    u_old[0] = u_l[tid];
    u_old[nx_s - 1] = u_r[tid];

    for (int j = j0l, k = j0; j <= j1l; j++, k++)
    {
      x = k * dx;

      u[j] = (alpha / 2.0) * (u_old[j - 1] + u_old[j + 1]) + (1.0 - alpha) * u_old[j] + dt * g(x);
    }
#pragma omp barrier

    // Set boundary exchange array
    if (tid != nthreads - 1)
    {
      u_l[tid + 1] = u[j1l];
    }
    if (tid != 0)
    {
      u_r[tid - 1] = u[1];
    }

#pragma omp barrier
    u[0] = u_l[tid];
    u[nx_s - 1] = u_r[tid];

    // Loop in time
    for (int i = 2; i <= m; i++)
    {

      for (int j = j0l; j <= j1l; j++)
      {
        u_new[j] =
            alpha * (u[j - 1] + u[j + 1]) + 2.0 * (1.0 - alpha) * u[j] - u_old[j];
      }
      // Reset boundary exchange array
      if (tid != nthreads - 1)
      {
        u_l[tid + 1] = u_new[j1l];
      }
      if (tid != 0)
      {
        u_r[tid - 1] = u_new[1];
      }
// Set boundary for next iteration
#pragma omp barrier

      u_new[0] = u_l[tid];
      u_new[nx_s - 1] = u_r[tid];

#pragma omp master
      t += dt;

#if (SAVEFILE == 1)
      // Save data to file
      x = x0;
      for (int j = 0; j < nx_s; j++)
      {
        outfile << "  " << t << "  " << x
                << "  " << u_new[j] << "\n";
        x += dx;
      }
      outfile << "\n";
#endif

      // Shuffle array in time
      std::swap(u, u_old);
      std::swap(u, u_new);
    }

#if (SAVEFILE == 1)
    // Close file
    outfile.close();
#endif
  }

  return 0;
}

// Initial condition
double f(double x)
{
  double value;

  if (0.4 <= x && x <= 0.60)
  {
    value = (x - 0.4) * (0.60 - x);
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
