#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <assert.h>

double complex integrate(double complex (*f)(double complex), double complex a, double complex b, double complex f2, double complex f3, double delta, double epsilon, int numrecs, double * errorestimate){
	assert(numrecs<100000);			// If numrecs goes beyond the threshold, we should probably tune the paramters.
	double cnorm(double complex z){
		double x = creal(z), y = cimag(z);
		return sqrt(x*x+y*y);
	}
	double complex f1 = f(a+(b-a)/6), f4 = f(a+5.0*(b-a)/6); 
	double complex Q = (b-a)/6*(2.0*f1+f2+f3+2.0*f4), q = (b-a)/4*(f1+f2+f3+f4), err=cnorm(Q-q);	// following eqs. (44), (45) and (46) in the notes.
	if (cnorm(err) < delta+epsilon*cnorm(Q)) {*errorestimate = err; return Q;}			// following eq. (46)
	else return integrate(f,a,(a+b)/2,f1,f2,delta/sqrt(2),epsilon,numrecs+1,errorestimate)+
		    integrate(f,(a+b)/2,b,f3,f4,delta/sqrt(2),epsilon,numrecs+1,errorestimate);
}

double complex recadapter(double complex (*f)(double complex), double complex a, double complex b, double delta, double epsilon, double * errorestimate) {
	double complex f2 = f(a+2.0*(b-a)/6), f3 = f(a+4.0*(b-a)/6); 					// following eqs. 
	int numrecs = 0;
	return integrate(f,a,b,f2,f3,delta,epsilon,numrecs,errorestimate);
}

double complex clenshaw_curtis(double complex (*f)(double complex), double complex a, double complex b, double delta, double epsilon, double * errorestimate){
	// variable transformation x -> (a+b)/2+(a-b)/2*cos(t) makes \int_a^b f(x)dx on the right form.
	double complex g(double complex t){return f((a+b)/2+(a-b)/2*ccos(t))*csin(t)*(b-a)/2;}
	return recadapter(g,0,M_PI,2*delta,2*epsilon,errorestimate);
}

