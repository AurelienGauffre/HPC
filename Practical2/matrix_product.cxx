#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include "omp.h"

using namespace Eigen;
using namespace std;
static int nb_threads = 8;

MatrixXd matrix_product(MatrixXd A, MatrixXd B, int n);
MatrixXd matrix_product_parallel(MatrixXd A, MatrixXd B, int n);
double norm(MatrixXd A, int n);

int main()
{
  int n = 400;
  /// nxn Matrix filled with random numbers between (-1,1)
  MatrixXd A = MatrixXd::Random(n,n);
  MatrixXd B = MatrixXd::Random(n,n);
  MatrixXd C = MatrixXd::Zero(n,n);
  MatrixXd D = MatrixXd::Zero(n,n);

  auto start = std::chrono::high_resolution_clock::now();
  C = matrix_product(A, B, n);
  auto stop = std::chrono::high_resolution_clock::now();
  auto elapsed =  std::chrono::duration<double>(stop-start).count();
  cout <<"Sequential : " <<elapsed << endl;


  start = std::chrono::high_resolution_clock::now();
  D = matrix_product_parallel(A, B, n);
  stop = std::chrono::high_resolution_clock::now();
  elapsed =  std::chrono::duration<double>(stop-start).count();
  cout <<"Parallel : " <<elapsed << endl;

  MatrixXd diff = C-D;
  double norm_diff = norm(diff,n);
  cout << "Norme de C-D : " << norm_diff << endl;
}

MatrixXd matrix_product(MatrixXd A, MatrixXd B, int n){
  MatrixXd C = MatrixXd::Zero(n,n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      double sum = 0;
      for(int k=0; k<n; k++) {
        sum += A(i,k) * B(k,j);
      }
      C(i,j) = sum;
    }
  }
  return C;
}

MatrixXd matrix_product_parallel(MatrixXd A, MatrixXd B, int n){
  int nb_blocks = 4;
  int m = n/nb_blocks;
  MatrixXd C = MatrixXd::Zero(n,n);
  #pragma omp parallel for 
  for(int i=0; i<nb_blocks; i++) {
    for(int j=0; j<nb_blocks; j++) {
      MatrixXd sum = MatrixXd::Zero(m,m);
      for(int k=0; k<nb_blocks; k++) {
        sum += matrix_product(A.block(i*m,k*m,m,m),B.block(k*m,j*m,m,m),m);
      }
      C.block(i*m,j*m,m,m) = sum;
    }
  }
  return C;
}

double norm(MatrixXd A, int n){
  double sum = 0;
  for (int i = 0; i< n; i++){
    for (int j = 0; j< n; j++){
      sum += A(i,j)*A(i,j);
    }
  }
  return sum;
}
