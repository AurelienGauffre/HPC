/*!
export OMP_NUM_THREADS = 8
 * \file   hello.cxx
 * \brief  Displays "Hello world"
 *
 * Usage: ./hello
 */
#include "omp.h"
#include <iostream>
#include <chrono>
int main ()
{
  int a = 0 ;


  for(int i = 0; i<3;i++)
  {
  #pragma omp parallel
  {
    std::cout<<omp_get_thread_num()<<std::endl;

  }
  }
  std::cout<<omp_get_num_threads()<<std::endl;
  return 0;
}
