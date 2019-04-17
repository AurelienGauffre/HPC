#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <ctime>
#include <cstdlib>
#include "omp.h"
#include <chrono>
#define NUM_THREADS 8

int main ( );
double f ( double x );
double g ( double x );
void timestamp ( );

int main (){

    // Number of dsicretisation point in time
    # define m 200

    // Number of space step
    # define n 800

    // Define array holding the amplitude of the string
    double u[m+1][n+1];
    double uparal[m+1][n+1];

    double alpha;
    double c = 0.25;
    std::ofstream command_unit;
    std::ofstream outfile;

    // Time interval
    double t;
    double t1 = 0.0;
    double t2 = 3.0;

    // Space interval
    double x;
    double x1 = 0.0;
    double x2 = 1.0;

    // Discretization steps
    double dt;
    double dx;

    dx = ( x2 - x1 ) / ( double ) n;
    dt = ( t2 - t1 ) / ( double ) m;
    alpha = pow ( c * dt / dx, 2 );

    // Initial condition
    u[0][0] = 0.0;
    for (int j = 1; j < n; j++ ) {
        x = j * dx;
        u[0][j] = f ( x );
    }
    u[0][n] = 0.0;

    // Bootstrap for the first time step
    u[1][0] = 0.0;
    for ( int j = 1; j < n; j++ ) {
        x = j * dx;
        u[1][j] = ( alpha / 2.0 ) * (u[0][j-1] + u[0][j+1])
        + ( 1.0 - alpha ) * u[0][j] + dt * g ( x );
    }
    u[1][n] = 0.0;

    for (int i = 0; i<=1; i++) {
        for (int j = 0; j<n+1; j++) {
            uparal[i][j] = u[i][j];
        }
    }
////////////////////////// SEQUENTIAL VERSION ////////////////////////
    auto start1 = std::chrono::high_resolution_clock::now();
    // Loop in time
    for ( int i = 2; i <= m; i++ ) {
        u[i][0] = 0.0;
        for ( int j = 1; j < n; j++ ) {
            u[i][j] =
            alpha   * (u[i-1][j-1] + u[i-1][j+1])
            + 2.0 * ( 1.0 - alpha ) * u[i-1][j]
            - u[i-2][j];
        }
        u[i][n] = 0.0;
    }
    auto stop1 = std::chrono::high_resolution_clock::now();
    auto elapsed1 =  std::chrono::duration<double>(stop1-start1).count();
    std::cout <<"Sequential : "<< elapsed1 << '\n';



////////////////////////// PARALLEL  VERSION ////////////////////////

    double edges[2*NUM_THREADS];
    int n_loc = n/NUM_THREADS;

    auto start2 = std::chrono::high_resolution_clock::now();

      #pragma omp parallel
      {
      int tid = omp_get_thread_num();

      double uparal_loc_old[n_loc];
      double uparal_loc[n_loc];
      double uparal_loc_new[n_loc];

      int j0 = tid*n_loc;
      int j1 = j0 + n_loc - 1;
      // Initial condition
      for (int j = j1; j < j2; j++ ) {
          x = j * dx;
          uparal_loc_old[j] = f ( x );
      }

      //Update edges
      edges[2*tid]=uparal_loc_old[0];
      edges[2*tid+1]=uparal_loc_old[n_loc-1];

      edges[0] = 0.0;
      edges[NUM_THREADS-1] = 0.0;


      // Bootstrap for the first time step
      for ( int j = j0+1; j < j1-1; j++ ) {
          x = j * dx;
          uparal_loc[j] = ( alpha / 2.0 ) * (uparal_loc_old[j-1] + uparal_loc_old[j+1])
          + ( 1.0 - alpha ) * uparal_loc_old[j] + dt * g ( x );
      }
      if (tid == 0){
        uparal_loc[j0] = 0.0;
      } else {
        x = j0 * dx;
        uparal_loc[j0] = ( alpha / 2.0 ) * (edges[tid*2-1] + uparal_loc_old[j0+1])
        + ( 1.0 - alpha ) * uparal_loc_old[j0] + dt * g ( x );
      }
      if (tid == NUM_THREADS){
        uparal_loc[j1-1] = 0.0;
      } else {
        x = (j1-1) * dx;
        uparal_loc[j1-1] = ( alpha / 2.0 ) * (uparal_loc_old[j1-2] + edges[tid*2+2])
        + ( 1.0 - alpha ) * uparal_loc_old[j1-1] + dt * g ( x );
      }

      std::swap(uparal_loc,uparal_loc_new);
      std::swap(uparal_loc_old,uparal_loc);

      //Update edges
      edges[2*tid]=uparal_loc[0];
      edges[2*tid+1]=uparal_loc[n_loc-1];

      //Loop in time
      for ( int i = 2; i <= m; i++ ) {
        if (tid == 0){
          uparal_loc_new[0] = 0.0;
        } else {
          uparal_loc_new[0] = alpha   * (edges[2*tid-1] + uparal_loc[j+1])
            + 2.0 * ( 1.0 - alpha ) * uparal_loc[j] - uparal_old[j];
        }

        #pragma omp  for
        for ( int j = 1; j < n_loc-1; j++ ) {
          uparal_loc_new[j] = alpha   * (uparal_loc[j-1] + uparal_loc[j+1])
            + 2.0 * ( 1.0 - alpha ) * uparal_loc[j] - uparal_old[j];
        }

        if (tid == NUM_THREADS) {
          uparal_loc_new = 0.0;
        } else {
          uparal_loc_new[n_loc-1] = alpha   * (uparal_loc[j-1] + edges[2*tid+2])
            + 2.0 * ( 1.0 - alpha ) * uparal_loc[j] - uparal_old[j];
        }

        std::swap(uparal_loc,uparal_loc_new);
        std::swap(uparal_loc_old,uparal_loc);

        #pragma omp barrier
        //Update edges
        edges[2*tid]=uparal_loc[0];
        edges[2*tid+1]=uparal_loc[n_loc-1];
      }

      }

    auto stop2 = std::chrono::high_resolution_clock::now();
    auto elapsed2 =  std::chrono::duration<double>(stop2-start2).count();
    std::cout <<"Parallel :   " <<elapsed2 << '\n';







    // Norm
    double norm = 0;
    for (int i = 0; i< m+1; i++){
        for (int j = 0; j< n+1; j++){
            double a = u[i][j]-uparal[i][j];
            norm += a*a;
        }
    }

std::cout << "L2 Norm difference  : " << norm << '\n';

//  Write data file.
outfile.open ( "string_data.dat" );

t = 0.; x = 0.;
for ( int i = 0; i <= m; i++ ) {
    for ( int j = 0; j <= n; j++ ) {
        outfile << "  " << x << "  " << t
        << "  " << u[i][j] << "\n";
        x += dx;
    }
    t += dt; x = 0.;
    outfile << "\n";
}
outfile.close ();

return 0;
# undef m
# undef n
}

// Initial condition
double f ( double x ) {
    double value;

    if ( 0.25 <= x && x <= 0.50 ) {
        value = ( x - 0.25 )*( 0.50 - x );
    }
    else {
        value = 0.0;
    }

    return value;
}

// Derivative of the intial condition
double g ( double x ) {
    double value;

    value = 0.0;
    return value;
}
