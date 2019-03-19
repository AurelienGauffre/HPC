/*!
 * \file  pi_MonteCarlo.cxx
 * \brief Compute pi using Monte Carlo
 */


// Display to std
#include <iostream>
// Value of PI
#include <cmath>
// Mesure time
#include <chrono>
// Generate random number
#include <random>

int main()
{
  // Number of Monte-Carlo tries
  int niter = 1000000;
  // x,y coordinates from the random tries
  double x,y;
  // Counter for number of good tries
  int count=0;
  // Approximate value of pi
  double pi;
  // Time to Initialize random number generator
  time_t t;

  // Initialize random generator
  srand48(((unsigned)time(&t)));
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);

  auto start = std::chrono::high_resolution_clock::now();
  for (int i=0; i<=niter; ++i)
  {
    // Get random x and y coordinates
    x = dis(gen);
    y = dis(gen);

    // Checks if point lies inside unit circle
    if (((x*x)+(y*y)) <=1)
    {
      //if it is, consider it a valid random point
      ++count;
    }
  }
  auto stop = std::chrono::high_resolution_clock::now();

  // Pi is equal to the ratio of hits to the total number of tries
  pi = ((double)count/(double)niter)*4.0;
  auto elapsed =  std::chrono::duration<double>(stop-start).count();
  std::cout << "The value of PI is " << pi << " with error ";
  std::cout << fabs(pi - M_PI)<< " and was obtained in ";
  std::cout << elapsed << " seconds."<<std::endl;

  return EXIT_SUCCESS;
}
