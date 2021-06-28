#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void randomx(int dim, double * a, double * b, double * x);
double vandercorput(int n, int base);
void halton(int n, int dim, double * a, double * b, double * x, int errorestimate);
void plainmc(int dim, double (*f)(int dim, double * x), double * a, double * b, int N, double * res, double * err);
void haltonmc(int dim, double (*f)(int dim, double * x), double * a, double * b, int N, double * res, double * err);

#endif
