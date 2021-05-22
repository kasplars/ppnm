#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include "library.h"
#include "ann.h"

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

ann * ann_alloc(int n,double(*f)(double),double(*f_dif)(double),double(*f_int)(double)){
	ann * network = malloc(sizeof(ann));
	network->n=n;
	network->f=f;
	network->f_dif=f_dif;
	network->f_int=f_int;
	network->params = gsl_vector_alloc(3*n);
	return network;
}

ann * ann_alloc_diffeq(int n,double(*phi)(double,double,double,double),double(*f)(double),double(*f_dif)(double),double(*f_difdif)(double)){
	ann * network = malloc(sizeof(ann));
	network->n=n;
	network->phi=phi;
	network->f=f;
	network->f_dif=f_dif;
	network->f_difdif=f_difdif;
	network->params = gsl_vector_alloc(3*n);
	return network;
}

void ann_free(ann * network){
	gsl_vector_free(network->params);
	free(network);
}

double ann_response(ann * network,double x){
	double response = 0;
	// params ordering = {ai,bi,wi,ai+1,bi+1,...}
	double a, b, w;
	for (int i = 0; i<network->n; i++) {
		a = gsl_vector_get(network->params,3*i);
		b = gsl_vector_get(network->params,3*i+1);
		w = gsl_vector_get(network->params,3*i+2);
		response += network->f((x-a)/b)*w;
	}
	return response;
}

double ann_response_dif(ann * network,double x){
	double response = 0;
	double a, b, w;
	// params ordering = {ai,bi,wi,ai+1,bi+1,...}
	for (int i = 0; i<network->n; i++) {
		a = gsl_vector_get(network->params,3*i);
		b = gsl_vector_get(network->params,3*i+1);
		w = gsl_vector_get(network->params,3*i+2);
		response += network->f_dif((x-a)/b)*w/b;
	}
	return response;
}

double ann_response_difdif(ann * network,double x){
	double response = 0;
	double a, b, w;
	// params ordering = {ai,bi,wi,ai+1,bi+1,...}
	for (int i = 0; i<network->n; i++) {
		a = gsl_vector_get(network->params,3*i);
		b = gsl_vector_get(network->params,3*i+1);
		w = gsl_vector_get(network->params,3*i+2);
		response += network->f_difdif((x-a)/b)*w/b/b;
	}
	return response;
}

double ann_response_int(ann * network,double x,double initialx){
	double response = 0;
	double a, b, w;
	// params ordering = {ai,bi,wi,ai+1,bi+1,...}
	for (int i = 0; i<network->n; i++) {
		a = gsl_vector_get(network->params,3*i);
		b = gsl_vector_get(network->params,3*i+1);
		w = gsl_vector_get(network->params,3*i+2);
		response += network->f_int((x-a)/b)*w*b - network->f_int((initialx-a)/b)*w*b;
	}
	return response;
}

void ann_train(ann * network,gsl_vector * xs, gsl_vector * ys){
	void costfunc(gsl_vector * p, double * Cp) {
		*Cp = 0;
		gsl_vector_memcpy(network->params,p);
		for (int i = 0; i<xs->size; i++) {
		*Cp += pow(ann_response(network,gsl_vector_get(xs,i))-gsl_vector_get(ys,i),2);
		}	
	}
	gsl_vector * p = gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	int steps = qnewton(costfunc,p,1e-4);
	printf("Number of steps in training the NN: %i\n",steps);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}

void ann_train_diffeq(ann * network, double a, double b, double c, double yc, double dyc){
	void costfunc(gsl_vector * p, double * Cp) {
		double integrand(double x, void * parameters) {
			return pow(network->phi(x,ann_response(network,x),ann_response_dif(network,x),ann_response_difdif(network,x)),2);
		}
		int u = 1000; double err;
		gsl_integration_workspace * w = gsl_integration_workspace_alloc(u);
		gsl_function F;
		F.function = &integrand;
		gsl_vector_memcpy(network->params,p);
		gsl_integration_qags(&F,a,b,0.001,0.01,u,w,Cp,&err);
		*Cp += (b-a)*(pow(ann_response(network,c) - yc,2) + pow(ann_response_dif(network,c) - dyc,2));
		gsl_integration_workspace_free(w);
	}
	gsl_vector * p = gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	int steps = qnewton(costfunc,p,1e-4);
	printf("Number of steps in training the NN: %i\n",steps);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}
