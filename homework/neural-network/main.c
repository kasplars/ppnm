#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <float.h>
#include "library.h"
#include "ann.h"

#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

double gausswavelet(double x) {return x*exp(-pow(x,2));}
double gauss(double x) {return exp(-pow(x,2));}
double wavelet(double x) {return cos(5*x)*exp(-pow(x,2));}

double f(double x) {return x*x*x+1.0;}

void linspace(gsl_vector * xs, double a, double b, int numpoints)  {
	double dx = (b-a)/(numpoints-1);
	for (int i = 0; i<numpoints; i++) {
		gsl_vector_set(xs,i,a+i*dx);
	}
}

int main(){
	int n = 3;
	ann * network=ann_alloc(n,gausswavelet);
	
	int numpoints = 30; double a = -1, b = 1;
	gsl_vector * xs = gsl_vector_alloc(numpoints);
	gsl_vector * ys = gsl_vector_alloc(numpoints);
	
	FILE * data = fopen("data.txt","w");	

	linspace(xs,a,b,numpoints);
	fprintf(data,"#index 0\n");
	for (int i = 0; i<numpoints; i++) {
		gsl_vector_set(ys,i,f(gsl_vector_get(xs,i)));
		fprintf(data,"%10g %10g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i));
	}
	fprintf(data,"\n\n\n");
	
	gsl_vector_set_all(network->params,1.0);
	
	ann_train(network,xs,ys);
	vector_print("p=",network->params);
	
	int plotnum = 100; 
	
	gsl_vector * plotxs = gsl_vector_alloc(plotnum);
	linspace(plotxs,a,b,plotnum);
	double t, trained, real;
	
	fprintf(data,"#index 1\n");
	for (int i = 0; i<plotnum; i++) {
		t = gsl_vector_get(plotxs,i);
		trained = ann_response(network,t);
		real = f(t);
		fprintf(data,"%10g %10g %10g\n",t,trained,real);
	}
	
	fclose(data);
	ann_free(network);
	gsl_vector_free(xs);
	gsl_vector_free(ys);
	gsl_vector_free(plotxs);
return 0;
}
