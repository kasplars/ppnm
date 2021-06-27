#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <assert.h>
#include "auxfile.h"

/*
Created by Kasper Larsen, student number: 201806672.

My exam problem is number 6 since 72 % 22 = 6

PROBLEM STATEMENT:

- 6 - Adaptive integration of complex-valued functions

Implement an adaptive integrator which calculates the integral of a complex-valued function f(z) of a complex variable z 
along a straight line between two points in the complex plane. 

*/

int main(){
	// definitions
	int numevals = 0;							// for counting the number of evaluations
	double delta = 0.001; 							// absolute accuracy goal
	double epsilon = 0.001; 						// relative accuracy goal
	double complex a = 0, b = 1.0 + 2*I; 					// integration is from point a to point b 
	double err;								// estimation of the error of the integral value
	double complex f(double complex z){numevals++; return 1/(2*csqrt(z));}	// integrand function
	
	// introduction
	printf("\nA recursive adaptive integrator calculating integrals of complex-valued functions has been implemented in this set of files.\n");
	printf("Recursive adaptive integration with and without the Clenshaw-Curtis variable transformation quadrature is given.\n");
	printf("Tests of the adaptive integrators are given below:\n\n");
	
	// calculation of integral using ordinary recursive adaptive integration
	double complex integral = recadapter(f,a,b,delta,epsilon,&err);
	
	// output
	printf("Results from integration of f(z)=1/(2*sqrt(z)) from %g%+gi to %g%+gi:\n",creal(a),cimag(a),creal(b),cimag(b));
	printf("	calculated integral value:	%15.15g%+gi\n",creal(integral),cimag(integral));
	printf("	true integral value:		sqrt(1+2i)\n");
	printf("	error estimate: 		%15.15g\n",err);
	printf("	number of evaluations: 		%i\n",numevals);
	printf("	relative accuracy:		%5.15g\n",epsilon);
	printf("	absolute accuracy:		%5.15g\n",delta);
	printf("	method:				Ordinary recursive adaptive integration\n\n");
	
	// calculation of integral using the Clenshaw-Curtis variable transformation quadrature:
	numevals = 0;
	integral = clenshaw_curtis(f,a,b,delta,epsilon,&err);
	
	// output
	printf("\nResults from integration of f(z)=1/(2*sqrt(z)) from %g%+gi to %g%+gi:\n",creal(a),cimag(a),creal(b),cimag(b));
	printf("	calculated integral value:	%15.15g%+gi\n",creal(integral),cimag(integral));
	printf("	true integral value:		sqrt(1+2i)\n");
	printf("	error estimate: 		%15.15g\n",err);
	printf("	number of evaluations: 		%i\n",numevals);
	printf("	relative accuracy:		%5.15g\n",epsilon);
	printf("	absolute accuracy:		%5.15g\n",delta);
	printf("	method:				Adaptive integrator with Clenshaw-Curtis\n\n");
	
	
return 0;
}
