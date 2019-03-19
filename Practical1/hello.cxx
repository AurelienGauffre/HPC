/*!
export OMP_NUM_THREADS = 8
 * \file   hello.cxx
 * \brief  Displays "Hello world"
 *
 * Usage: ./hello
 */
#include "omp.h"
#include <iostream>

int main ()
{
  int a = 0 ;
  #pragma omp parallel
  {

  //std::cout << a<< std::endl ;
  std::cout<<"Hello World!"<<std::endl;
  std::cout<< omp_get_thread_num()<<std::endl;
  }
  return EXIT_SUCCESS;
}
