#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include <omp.h>

int main()
{

  int n = 10;

  /// nxn Matrix filled with random numbers between (-1,1)
  MatrixXd A = MatrixXd::Random(n,n);
  MatrixXd B = MatrixXd::Random(n,n);
  MatrixXd C = MatrixXd::Zero(n,n);

#pragma parallel omp for collapse(2)
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      double sum = 0;
      for(int k=0; k<n; k++) {
        sum = A(i,k) * B(k,j);
      }
      C(i,j) = sum;
    }
  }
}
