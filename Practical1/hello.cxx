/*!
In order to modify the number of threads used :
  export OMP_NUM_THREADS=8
 * \file   hello.cxx
 * \brief  Make each thread diplays its number
 *
 * Usage: ./hello
 */
#include "omp.h"
#include <iostream>
#include <chrono>

int main ()
  omp_set_num_threads(8);
{

  #pragma omp parallel
  {
    std::cout<<omp_get_thread_num()<<std::endl;
  }
  std::cout<<omp_get_num_threads()<<std::endl; //
  return 0;
}
