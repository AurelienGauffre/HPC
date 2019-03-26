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
using namespace std;
int main ()

{
  int const nb_thread= 8;
  // Number of integration interval
  int num_steps = 100000;
  // Evaluation point
  double x;
  // Approximate value of pi
  double pi;
  // Approximate value of the integral
  double sum = 0.0;
  // Integration step
  double step = 1.0/(double) num_steps;

  double xtab[nb_thread] ;
 auto start = std::chrono::high_resolution_clock::now();

//   for (int i=1;i<= num_steps/nb_thread; i++)
//   {
//     #pragma omp parallel
//     {
//       int nt=omp_get_thread_num();
//       xtab[nt]=4.0/(1.0+ ((num_steps/nb_thread*nt+i-0.5)*step)*(num_steps/nb_thread*nt+i-0.5)*step);
//
//     }
//
//     for(int j = 0; j<8;j++)
//     {
//       sum = sum + xtab[j];
//     }
// }
#pragma omp for
for (int i=1;i<= num_steps; i++)
{

    int nt=omp_get_thread_num();
    xtab[nt]+=(double)(4.0/(1.0+ ((i-0.5)*step)*(i-0.5)*step));


}
for(int j = 0; j<8;j++)
{
  sum = sum + xtab[j];
}

    for (int i=1;i<= num_steps; i++)
    {
    x = (i-0.5)*step;
    sum = sum + 4.0/(1.0+x*x);

    }
 auto stop = std::chrono::high_resolution_clock::now();

  pi = step * sum;
 auto elapsed =  std::chrono::duration<double>(stop-start).count();
  std::cout << "The value of PI after "<< num_steps << " steps is";
  std::cout << pi << " with error ";
  std::cout << fabs(pi/2 - M_PI)<< " and was obtained in ";
 std::cout << elapsed << " seconds."<<std::endl;

  return 0;
}
