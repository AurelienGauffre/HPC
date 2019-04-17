/*!
 * \file   hello_omp.cxx
 * \brief  Displays "Hello world" in parallel
 *
 * Usage: ./hello
 */

#include <iostream>
#include <omp.h>

int main()
{
  int nthreads = 0;

#pragma omp parallel
  {
    int id = omp_get_thread_num();
#pragma omp master
    {
      int nthreads = omp_get_num_threads();
    }
    std::cout << "Hello World from thread = " << id;
    std::cout << " with " << nthreads << " threads" << std::endl;
  }

  std::cout << "all done, with hopefully " << nthreads << " threads." << std::endl;
}
