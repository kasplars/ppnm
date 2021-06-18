#include <math.h>
#include <stdio.h>
#include <complex.h>
#include "auxfile.h"

/*
Created by Kasper Larsen, student number: 201806672.

My exam problem is number 6 since 72 % 22 = 6

PROBLEM STATEMENT:

- 6 - Adaptive integration of complex-valued functins

Implement an adaptive integrator which calculates the integral of a complex-valued function f(z) of a complex variable z along a straight line between two points in the complex plane. 
*/

int main(){
	int numevals = 0;							// for counting the number of evaluations
	double complex f(double complex z){numevals++; return 1/(2*csqrt(z));}	// integrand function
	double complex a = 0, b = 1.0 + 2*I; 					// integration is from point a to point b 
	
	double delta = 0.001; 							// absolute accuracy goal
	double epsilon = 0.001; 						// relative accuracy goal
	
	double complex integral = recadapter(f,a,b,delta,epsilon); 		// calculation of integral using recursive adaptive integration
	
	// output
	printf("Integral of f(z) from %g%+gi to %g%+gi is %g%+gi. Number of evaluations: %i.\n",creal(a),cimag(a),creal(b),cimag(b),
		creal(integral),cimag(integral),numevals);
return 0;
}
