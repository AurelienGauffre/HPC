#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <ctime>
#include <cstdlib>

int main();
double f(double x);
double g(double x);
void timestamp();

int main()
{
  // Number of discretization point in space
  int m = 100;

  // Number of time step
  int n = 100;

  // Define array holding the amplitude of the string
  double u[m + 1][n + 1];

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

  // Initial condition
  u[0][0] = 0.0;
  for (int j = 1; j < n; j++)
  {
    x = j * dx;
    u[0][j] = f(x);
  }
  u[0][n] = 0.0;

  // Bootstrap for the first time step
  u[1][0] = 0.0;
  for (int j = 1; j < n; j++)
  {
    x = j * dx;
    u[1][j] = (alpha / 2.0) * (u[0][j - 1] + u[0][j + 1]) + (1.0 - alpha) * u[0][j] + dt * g(x);
  }
  u[1][n] = 0.0;

  // Loop in time
#pragma omp parallel
  {
    for (int i = 2; i <= m; i++)
    {
      u[i][0] = 0.0;
#pragma omp for
      for (int j = 1; j < n; j++)
      {
        u[i][j] =
            alpha * (u[i - 1][j - 1] + u[i - 1][j + 1]) + 2.0 * (1.0 - alpha) * u[i - 1][j] - u[i - 2][j];
      }
      u[i][n] = 0.0;
    }
  }
  //  Write data file.
  outfile.open("string_data.dat");

  t = 0.;
  x = 0.;

  for (int i = 0; i <= m; i++)
  {
    for (int j = 0; j <= n; j++)
    {
      outfile << "  " << x << "  " << t
              << "  " << u[i][j] << "\n";
      x += dx;
    }
    t += dt;
    x = 0.;
    outfile << "\n";
  }
  outfile.close();

  return 0;
#undef m
#undef n
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

// Derivative of the intial condition
double g(double x)
{
  double value;

  value = 0.0;
  return value;
}
