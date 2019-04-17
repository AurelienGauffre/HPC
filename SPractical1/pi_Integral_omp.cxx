/*!
 * \file  pi_Integral_omp.cxx
 * \brief Compute pi using numerical integration
 */

// Display to std
#include <iostream>
// Value of PI
#include <cmath>
// Mesure time
#include <chrono>

int main()
{
  // Number of integration interval
  int num_steps = 10000;
  // Evaluation point
  double x;
  // Approximate value of pi
  double pi;
  // Approximate value of the integral
  double sum = 0.0;
  // Integration step
  double step = 1.0 / (double)num_steps;

  auto start = omp_get_wtime();
#pragma omp parallel for reduction(sum \
                                   : x)
  for (int i = 1; i <= num_steps; i++)
  {
    x = (i - 0.5) * step;
    sum = sum + 4.0 / (1.0 + x * x);
  }
  auto stop = omp_get_wtime();

  pi = step * sum;
  auto elapsed = stop - start;
  std::cout << "The value of PI after " << num_steps << "steps is";
  std::cout << pi << " with error ";
  std::cout << fabs(pi - M_PI) << " and was obtained in ";
  std::cout << elapsed << " seconds." << std::endl;

  return EXIT_SUCCESS;
}
