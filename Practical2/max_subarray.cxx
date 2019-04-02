#include <iostream>
#include <limits>
#include <Eigen/Dense>
#include <chrono>
using namespace Eigen;

void kadane(const VectorXd & array, double &maxSum, int &l, int &r)
{
  maxSum = -std::numeric_limits<double>::infinity();
  l = 0;
  r = 0;
  double sum = 0;
  int currentStartIndex = 0;
  for(int i = 0; i<array.size(); ++i){
    sum += array(i);
    if (sum > maxSum){
      maxSum = sum;
      l = currentStartIndex;
      r = i;
    }
    if(sum < 0) {
      sum= 0;
      currentStartIndex = i + 1 ;
    }
  }
}



void maxSubarray2D(const MatrixXd & array,
  double &maxSum, int &left, int &right, int &top, int &bottom){

  maxSum = -std::numeric_limits<double>::infinity();
  left = -1;
  right = -1;
  top = -1;
  bottom = -1;
  double sum = 0;
  int start, finish;

  for(int i = 0; i<array.rows(); ++i){
    VectorXd temp = VectorXd::Zero(array.cols());
    for(int j = i; j<array.rows(); ++j){
      for(int k = 0; k<array.cols(); ++k){
        temp(k) += array(j,k);
      }

      kadane(temp, sum, start, finish);
      if(sum > maxSum){
        maxSum = sum;
        left = i;
        right = j;
        top = start;
        bottom = finish;

      }
    }
  }
}

void maxSubarray2DparallelNaive(const MatrixXd & array,
  double &maxSum, int &left, int &right, int &top, int &bottom){

  maxSum = -std::numeric_limits<double>::infinity();
  left = -1;
  right = -1;
  top = -1;
  bottom = -1;
  double sum = 0;
  int start, finish;

  #pragma omp parallel for
  for(int i = 0; i<array.rows(); ++i){
    VectorXd temp = VectorXd::Zero(array.cols());
    for(int j = i; j<array.rows(); ++j){
      for(int k = 0; k<array.cols(); ++k){
        temp(k) += array(j,k);
      }
      #pragma omp critical
      {
      kadane(temp, sum, start, finish);
      if(sum > maxSum){
        maxSum = sum;
        left = i;
        right = j;
        top = start;
        bottom = finish;
      }
      }
    }
  }
}
void maxSubarray2DparallelOpti(const MatrixXd & array,
  double &maxSum, int &left, int &right, int &top, int &bottom){

  maxSum = -std::numeric_limits<double>::infinity();
  left = -1;
  right = -1;
  top = -1;
  bottom = -1;
  double sum = 0;
  int start, finish;

  int local_left[8];
  int local_right[8];
  int local_top[8];
  int local_bottom[8];
  double local_maxSum[8];


  #pragma omp parallel  //firstprivate(left,right,top,bottom,sum,start, finish)
{
  int tid = omp_get_thread_num();
  double local_sum = 0 ;
  int local_start;
  int local_finish;
  local_left[tid]= -1;
  local_right[tid]= -1;
  local_top[tid]= -1;
  local_bottom[tid]= -1;
  local_maxSum[tid]=0;
  #pragma omp for
  for(int i = 0; i<array.rows(); ++i){
    VectorXd temp = VectorXd::Zero(array.cols());
    for(int j = i; j<array.rows(); ++j){

      for(int k = 0; k<array.cols(); ++k){

        temp(k) += array(j,k);
      }
      kadane(temp, local_sum, local_start, local_finish);
      if(local_sum > local_maxSum[tid]){
        local_maxSum[tid] = local_sum;
        local_left[tid] = i;
        local_right[tid] = j;
        local_top[tid] = local_start;
        local_bottom[tid] = local_finish;
      }
    }
  }
  }
  double maxi = local_maxSum[0];
  int maxtid = 0 ;
  for(int p=1; p<8;p++){
       if (local_maxSum[p]>maxi){
         maxi = local_maxSum[p];
         maxtid = p ;
       }
  }

  maxSum = local_maxSum[maxtid] ;
  left = local_left[maxtid] ;
  right = local_right[maxtid] ;
  bottom = local_bottom[maxtid] ;
  top = local_top[maxtid] ;
}

int main()
{
  /// Size of the matrix
  int n = 200;

  /// nxn Matrix filled with random numbers between (-1,1)
  srand((unsigned int) time(0));
  MatrixXd m = MatrixXd::Random(n,n);
  double maxSum;
  int left, right, top, bottom;


  auto start = std::chrono::high_resolution_clock::now();
  maxSubarray2D(m, maxSum, left, right, top, bottom);
  auto stop = std::chrono::high_resolution_clock::now();
  auto elapsed =  std::chrono::duration<double>(stop-start).count();

  std::cout <<" Sequential :     "<< maxSum<<" "<<" "<< left<<" " <<right<< "  "<<top<<" "<< bottom <<" Time "<< elapsed<< std::endl;


  start = std::chrono::high_resolution_clock::now();
  maxSubarray2DparallelNaive(m, maxSum, left, right, top, bottom);
  stop = std::chrono::high_resolution_clock::now();
  elapsed =  std::chrono::duration<double>(stop-start).count();

  std::cout <<" Naive Parallel : "<< maxSum<<" "<<" "<< left<<" " <<right<< "  "<<top<<" "<< bottom <<" Time "<< elapsed<< std::endl;





  start = std::chrono::high_resolution_clock::now();
  maxSubarray2DparallelOpti(m, maxSum, left, right, top, bottom);
  stop = std::chrono::high_resolution_clock::now();
  elapsed =  std::chrono::duration<double>(stop-start).count();

  std::cout <<" Opti  Parallel : "<< maxSum<<" "<<" "<< left<<" " <<right<< "  "<<top<<" "<< bottom <<" Time "<< elapsed<< std::endl;

}
