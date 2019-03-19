/*!
 * \file  mandel.cxx
 * \brief Compute Mandelbrot set area
 */

// Display to std
#include <iostream>
// Value of PI
#include <cmath>
// OpenMP
#include <omp.h>

# define NPOINTS 1000
# define MAXITER 1000

void testpoint(void);

struct d_complex
{
  double r;
  double i;
};

struct d_complex c;
int numoutside = 0;

int main()
{
  double area, error, eps  = 1.0e-5;

  /* Loop over grid of points in the complex plane which contains
   * the Mandelbrot set, testing each point to see whether
   * it is inside or outside the set.
   */
#pragma omp parallel for default(shared) private(c,eps)
  for (int i=0; i<NPOINTS; i++)
  {
    for (int j=0; j<NPOINTS; j++)
    {
      c.r = -2.0+2.5*(double)(i)/(double)(NPOINTS)+eps;
      c.i = 1.125*(double)(j)/(double)(NPOINTS)+eps;
      testpoint();
    }
  }

  // Calculate area of set and error estimate and output the results
  area=2.0*2.5*1.125*(double)(NPOINTS*NPOINTS-numoutside)/(double)(NPOINTS*NPOINTS);
  error=area/(double)NPOINTS;

  std::cout<<"Area of Mandlebrot set :" << area <<" +/-" << error <<std::endl;
  std::cout<<"Correct answer should be around 1.510659" <<std::endl;

  return EXIT_SUCCESS;
}

  /*!
   * \brief   Test if point is in or out of the set
   * \details Does the iteration z=z*z+c, until |z| > 2 when point is known
   *          to be outside set. If loop count reaches MAXITER,
   *          point is considered to be inside the set
   * \param  void
   * \return void
   */
void testpoint(void)
{
  struct d_complex z;
  double temp;

  z=c;
  for (int iter=0; iter<MAXITER; iter++)
  {
    temp = (z.r*z.r)-(z.i*z.i)+c.r;
    z.i = z.r*z.i*2+c.i;
    z.r = temp;
    if ((z.r*z.r+z.i*z.i)>4.0)
    {
      numoutside++;
      break;
    }
  }
}

