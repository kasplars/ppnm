#ifndef AUXFILE_H
#define AUXFILE_H

double complex integrate(double complex (*f)(double complex), double complex a, double complex b, double complex f2, double complex f3, double delta, double epsilon, int numrecs, double * errorestimate);
double complex recadapter(double complex (*f)(double complex), double complex a, double complex b, double delta, double epsilon, double * errorestimate);
double complex clenshaw_curtis(double complex (*f)(double complex), double complex a, double complex b, double delta, double epsilon, double * errorestimate);

#endif
