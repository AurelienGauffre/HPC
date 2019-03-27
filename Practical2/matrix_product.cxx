#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include "omp.h"

using namespace Eigen;
using namespace std;
static int nb_threads = 8;

MatrixXd matrix_product(MatrixXd A, MatrixXd B, int n, int p, int q);
// void matrix_product_parallel(MatrixXd A, MatrixXd B, int n,  int p, int q);

int main()
{
  int n = 10;
  /// nxn Matrix filled with random numbers between (-1,1)
  MatrixXd A = MatrixXd::Random(n,n);
  MatrixXd B = MatrixXd::Random(n,n);
  MatrixXd C = MatrixXd::Zero(n,n);
  MatrixXd D = MatrixXd::Zero(n,n);

  auto start = std::chrono::high_resolution_clock::now();
  C = matrix_product(A, B, n, n, n);
  auto stop = std::chrono::high_resolution_clock::now();
  auto elapsed =  std::chrono::duration<double>(stop-start).count();
  cout <<"Sequential : " <<elapsed << endl;

  // start = std::chrono::high_resolution_clock::now();
  // D = matrix_product_parallel(A, B, n, n, n);
  // stop = std::chrono::high_resolution_clock::now();
  // elapsed =  std::chrono::duration<double>(stop-start).count();
  // cout <<"Parallel : " <<elapsed << endl;
  //
  // cout << C-D << endl;
}

MatrixXd matrix_product(MatrixXd A, MatrixXd B, int n, int p, int q){
  MatrixXd C = MatrixXd::Zero(n,n);
  #pragma omp parallel for
  for(int i=0; i<n; i++) {
    for(int j=0; j<q; j++) {
      double sum = 0;
      for(int k=0; k<p; k++) {
        sum = A(i,k) * B(k,j);
      }
      C(i,j) = sum;
    }
  }
  return C;
}

// void matrix_product_parallel(MatrixXd A, MatrixXd B, int n, int p, int q){
//   int m1 = n/2;
//   int m2 = n - m1;
//   matrixXD A0 = A.block(0,0,m1,m1);
//   matrixXD A1 = A.block(0,m1,m1,m2);
//   matrixXD A2 = A.block(m1,0,m2,m1);
//   matrixXD A3 = A.block(m1,m1,m2,m2);
//   matrixXD B0 = B.block(0,0,m1,m1);
//   matrixXD B1 = B.block(0,m1,m1,m2);
//   matrixXD B2 = B.block(m1,0,m2,m1);
//   matrixXD B3 = B.block(m1,m1,m2,m2);
//   matrixXD C0 = C.block(0,0,m1,m1);
//   matrixXD C1 = C.block(0,m1,m1,m2);
//   matrixXD C2 = C.block(m1,0,m2,m1);
//   matrixXD C3 = C.block(m1,m1,m2,m2);
//   matrix_product_parallel(A0,B0,C0,m1);
//   #pragma omp for
//   for(int i=0; i<n; i++) {
//     for(int j=0; j<n; j++) {
//       double sum = 0;
//       for(int k=0; k<n; k++) {
//         sum = A(i,k) * B(k,j);
//       }
//       C(i,j) = sum;
//     }
//   }
// }
