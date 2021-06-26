#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "functions.h"

double integrate(double (*f)(double), double a, double b, double f2, double f3, double delta, 
	double epsilon, int numrecs, double * errorestimate)
{
	double f1 = f(a+(b-a)/6), f4 = f(a+5.0*(b-a)/6); 
	double Q = (b-a)/6*(2.0*f1+f2+f3+2.0*f4), q = (b-a)/4*(f1+f2+f3+f4), err=fabs(Q-q);
	if (err < delta+epsilon*fabs(Q)) {*errorestimate = err; return Q;}
	else return integrate(f,a,(a+b)/2,f1,f2,delta/sqrt(2),epsilon,numrecs+1,errorestimate)+
		    integrate(f,(a+b)/2,b,f3,f4,delta/sqrt(2),epsilon,numrecs+1,errorestimate);
}

double recadapter(double (*f)(double), double a, double b, double delta, double epsilon, double * errorestimate) {
	if (isinf(a) == -1 && isinf(b) == 1) {return integrate_minftoinf(f,delta,epsilon,errorestimate);}
	if (isinf(a) == -1 && isinf(b) == 0) {return integrate_minftob(f,b,delta,epsilon,errorestimate);}
	if (isinf(a) == 0 && isinf(b) == 1) {return integrate_atoinf(f,a,delta,epsilon,errorestimate);}
	else {
	double f2 = f(a+2.0*(b-a)/6), f3 = f(a+4.0*(b-a)/6); int numrecs = 0;
	return integrate(f,a,b,f2,f3,delta,epsilon,numrecs,errorestimate);
	}
}

double clenshaw_curtis(double (*f)(double),double a,double b,double delta,double epsilon,double * errorestimate){
	// variable transformation x -> (a+b)/2+(a-b)/2*cos(t) makes \int_a^b f(x)dx on the right form.
	double g(double t){return f((a+b)/2+(a-b)/2*cos(t))*sin(t)*(b-a)/2;}
	return recadapter(g,0,M_PI,2*delta,2*epsilon,errorestimate);
}

double integrate_atoinf(double (*f)(double), double a, double delta, double epsilon, double * errorestimate) {
	double ftrans(double t){return f(a+(1.0-t)/t)/pow(t,2);}
	return recadapter(ftrans,0,1,delta,epsilon,errorestimate);
}

double integrate_minftob(double (*f)(double), double b, double delta, double epsilon, double * errorestimate) {
	double ftrans(double t){return f(b-(1.0-t)/t)/pow(t,2);}
	return recadapter(ftrans,0,1,delta,epsilon, errorestimate);
}

double integrate_minftoinf(double (*f)(double), double delta, double epsilon, double * errorestimate) {
	double ftrans(double t){return (f((1.0-t)/t)+f(-(1.0-t)/t))/pow(t,2);}
	return recadapter(ftrans,0,1,delta,epsilon,errorestimate);
}
