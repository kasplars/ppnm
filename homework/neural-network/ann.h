#include "gsl/gsl_vector.h"
#ifndef ANN_H
#define ANN_H

typedef struct {int n; double(*f)(double); double(*phi)(double,double,double,double); double(*f_dif)(double); double(*f_difdif)(double); double(*f_int)(double); gsl_vector * params;} ann;

ann * ann_alloc(int n,double(*f)(double),double(*f_dif)(double), double(*f_int)(double));

ann * ann_alloc_diffeq(int n,double(*phi)(double,double,double,double),double(*f)(double),double(*f_dif)(double),double(*f_difdif)(double));

void ann_free(ann * network);

double ann_response(ann * network,double x);

double ann_response_dif(ann * network,double x);

double ann_response_int(ann * network,double x,double a);

void ann_train(ann * network,gsl_vector * xs, gsl_vector * ys);

void ann_train_diffeq(ann * network, double a, double b, double c, double yc, double dyc);

#endif
