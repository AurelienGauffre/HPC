/*!
 * \file  pi_Integral.cxx
 * \brief Compute pi using numerical integration
 */

// Display to std
#include <iostream>
// Value of PI
#include <cmath>
// Mesure time
#include <chrono>
#include "omp.h"

int main ()
{
  int nb_thread= 8;
  // Number of integration interval
  int num_steps = 10000;
  // Evaluation point
  double x;
  // Approximate value of pi
  double pi;
  // Approximate value of the integral
  double sum = 0.0;
  // Integration step
  double step = 1.0/(double) num_steps;

  auto start = std::chrono::high_resolution_clock::now();
  for (int i=1;i<= num_steps/nb_thread; i++)
  {
    #pragma omp parallel
    {
      std::cout<<omp_get_thread_num()<<std::endl;
    sum = sum + 4.0/(1.0+(pow((8*i+omp_get_thread_num()-0.5)*step,2)));
    }


}

    // for (int i=1;i<= num_steps; i++)
    // {
    // x = (i-0.5)*step;
    // sum = sum + 4.0/(1.0+x*x);
    //
    // }
  auto stop = std::chrono::high_resolution_clock::now();

  pi = step * sum;
  auto elapsed =  std::chrono::duration<double>(stop-start).count();
  std::cout << "The value of PI after "<< num_steps << " steps is";
  std::cout << pi << " with error ";
  std::cout << fabs(pi - M_PI)<< " and was obtained in ";
  std::cout << elapsed << " seconds."<<std::endl;

  return EXIT_SUCCESS;
}
