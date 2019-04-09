#include <math.h>

#include "functions.h"

#define PI 4.*atan(1.0)

double f1(double x, double y, int ordre)
{
	switch (ordre)
	{
		case 1:
			return (-2.0 * ( 1.0 + y ) * ( 1.0 - y ));
		case 2:
			return (-2.0 * ( 1.0 + x ) * ( 1.0 - x ));
		default:
			return (( 1.0 - x * x ) * ( 1.0 - y * y ) );  // exact solution
	}
}

double f2(double x, double y, int ordre)
{
	switch (ordre)
	{
		case 1:
			return (-PI*PI*cos(PI*x)*cos(PI*y));
		case 2:
			return (-PI*PI*cos(PI*x)*cos(PI*y));
		default:
			return (cos(PI*x)*cos(PI*y));
	}
}

double f3(double x, double y, int ordre)
{
	switch (ordre)
	{
		case 1:
			return (-PI*PI*cos(PI*x));
		case 2:
			return (-PI*PI*cos(PI*y));
		default:
			return (cos(PI*x) + cos(PI*y));
	}
}

double f4(double x, double y, int ordre)
{
	switch (ordre)
	{
		case 1:
			return 0.;
		case 2:
			return 0.;
		default:
			return 2.;
	}
}

double f5(double x, double y, int ordre)
{
	switch (ordre)
	{
		case 1:
			return 0.;
		case 2:
			return 0.;
		default:
			if (x==0. | y==0.)
				return 2.;
			else
				return 0.;
	}
}
