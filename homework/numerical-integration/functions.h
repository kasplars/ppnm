#ifndef HAVE_FUNCTIONS_H
#define HAVE_FUNCTIONS_H

double integrate(double (*f)(double), double a, double b, double f2, double f3, double delta, 
	double epsilon, int numrecs, double * errorestimate);
double recadapter(double (*f)(double), double a, double b, double delta, double epsilon, double * errorestimate);
double clenshaw_curtis(double (*f)(double),double a,double b,double delta,double epsilon,double * errorestimate);
double integrate_atoinf(double (*f)(double), double a, double delta, double epsilon, double * errorestimate);
double integrate_minftob(double (*f)(double), double b, double delta, double epsilon, double * errorestimate);
double integrate_minftoinf(double (*f)(double), double delta, double epsilon, double * errorestimate);

#endif
