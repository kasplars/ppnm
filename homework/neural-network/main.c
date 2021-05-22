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
double gausswavelet_dif(double x) {return exp(-pow(x,2))-2*x*x*exp(-pow(x,2));}
double gausswavelet_difdif(double x) {return exp(-pow(x,2))*(4*x*x*x-6*x);}
double gausswavelet_int(double x) {return -exp(-pow(x,2))/2.0;}

double f(double x) {return sin(x);}
double Phi(double x, double y, double yd, double ydd) {return y+ydd;}

void linspace(gsl_vector * xs, double a, double b, int numpoints)  {
	double dx = (b-a)/(numpoints-1);
	for (int i = 0; i<numpoints; i++) {
		gsl_vector_set(xs,i,a+i*dx);
	}
}

int main(){
	// Part A and B
	
	printf("Part A and B\n");
	
	// parameters for program:
	int n = 3; 			// number of deep neurons
	int numpoints = 50;		// number of points to train from 
	double a = -3.14, b = 3.14;	// start and end
	int plotnum = 100; 		// number of points to create from the neural network function
	
	// allocate space for neural network and set default parameters
	ann * network = ann_alloc(n,gausswavelet,gausswavelet_dif,gausswavelet_int);
	for (int i = 0; i<n; i++) {
		gsl_vector_set(network->params,3*i,a+(b-a)*i/(n-1)); // displace uniformly on interval
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
	}
	
	// create training set and make data file
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
	
	// training neural network and printing output
	ann_train(network,xs,ys);
	vector_print("p=",network->params);
	
	// neural network functions are plotted here
	gsl_vector * plotxs = gsl_vector_alloc(plotnum);
	linspace(plotxs,a,b,plotnum);
	double xval, trained, trained_dif, trained_int, real;
	fprintf(data,"#index 1\n");
	for (int i = 0; i<plotnum; i++) {
		xval = gsl_vector_get(plotxs,i);
		trained = ann_response(network,xval);
		trained_dif = ann_response_dif(network,xval);
		trained_int = ann_response_int(network,xval,a);
		real = f(xval);
		fprintf(data,"%10g %10g %10g %10g %10g\n",xval,trained,trained_dif,trained_int,real);
	}
	fprintf(data,"\n\n\n");
	
	// free data
	ann_free(network);
	gsl_vector_free(xs);
	gsl_vector_free(ys);
	gsl_vector_free(plotxs);
	
	printf("The code correctly finds an approximation to the artificially created data set. It also finds\n");
	printf("the approximate derivative and anti-derivative.\n\n\n");
	
	// Part C
	
	printf("Part C\n");
	
	// parameters for program:
	n = 4; 			// number of deep neurons
	a = -3.14, b = 3.14;	// start and end
	plotnum = 100; 		// number of points to create from the neural network function
	double c = 0, yc = 0, dyc = 1;	// boundary conditions for differential equation
	
	// allocate space for neural network and set default parameters
	ann * dnetwork = ann_alloc_diffeq(n,Phi,gausswavelet,gausswavelet_dif,gausswavelet_difdif);
	for (int i = 0; i<n; i++) {
		gsl_vector_set(dnetwork->params,3*i,a+(b-a)*i/(n-1)); // displace uniformly on interval
		gsl_vector_set(dnetwork->params,3*i+1,1);
		gsl_vector_set(dnetwork->params,3*i+2,1);
	}
	
	// training neural network and printing output
	ann_train_diffeq(dnetwork,a,b,c,yc,dyc);
	vector_print("p=",dnetwork->params);
	
	// neural network functions are plotted here
	gsl_vector * plotts = gsl_vector_alloc(plotnum);
	linspace(plotts,a,b,plotnum);
	double t, trainedsolution, realsolution;
	fprintf(data,"#index 2\n");
	for (int i = 0; i<plotnum; i++) {
		t = gsl_vector_get(plotts,i);
		trainedsolution = ann_response(dnetwork,t);
		realsolution = f(t);
		fprintf(data,"%10g %10g %10g\n",t,trainedsolution,realsolution);
	}
	fprintf(data,"\n\n\n");
	
	printf("The code correctly finds an approximation to the solution of the provided dif. eq.\n\n");
	
	
	// free data
	fclose(data);
	ann_free(dnetwork);
	gsl_vector_free(plotts);
	
return 0;
}
