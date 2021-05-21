#include "gsl/gsl_vector.h"
#ifndef ANN_H
#define ANN_H

typedef struct {int n; double(*f)(double); double(*f_dif)(double); double(*f_int)(double); gsl_vector * params;} ann;
ann * ann_alloc(int n,double(*f)(double),double(*f_dif)(double), double(*f_int)(double));
void ann_free(ann * network);
double ann_response(ann * network,double x);
double ann_response_dif(ann * network,double x);
double ann_response_int(ann * network,double x,double a);
void ann_train(ann * network,gsl_vector * xs, gsl_vector * ys);

#endif
