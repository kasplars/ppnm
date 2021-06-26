#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include "functions.h"

double f(double x) {return 1.0/sqrt(x);}
double p(double x) {return 4.0*sqrt(1-pow(x,2));}
double gauss(double x) {return exp(-x*x);}

int main() {
	// Part A and B
	
	printf("Part A and B\n\n");

	double a = 0.0, b = 1.0, delta = 0.01, epsilon = 0.01, errorestimate; int numevals = 0;
	double function(double x){numevals++; return f(x);}
	double integral = recadapter(function,a,b,delta,epsilon,&errorestimate);
	
	printf("Integral of f(x) from %g to %g is %g not using Clenshaw-Curtis. Number of evaluations: %i.\n",a,b,integral,numevals);
	numevals = 0;
	
	integral = clenshaw_curtis(function,a,b,delta,epsilon,&errorestimate);
	
	printf("Integral of f(x) from %g to %g is %g using Clenshaw-Curtis. Number of evaluations: %i.\n",a,b,integral,numevals);
	numevals = 0;
	
	printf("So there is a lot of improvement in numbers of evaluations.\n");
	
	a = 0.0, b = 1.0, delta = 0.000000001, epsilon = 0.000000001;
	double function_new(double x){numevals++; return p(x);}
	double integral1 = clenshaw_curtis(function_new,a,b,delta,epsilon,&errorestimate); int numevals_aux = numevals; numevals = 0;
	double integral2 = recadapter(function_new,a,b,delta,epsilon,&errorestimate);
	printf("The integral of p(x) from %g to %g without and with Clenshaw-Curtis, respectively:\n",a,b);
	printf("%.15g; number of evaluations: %i; error: %.15g\n",integral1,numevals_aux,fabs(M_PI-integral1));
	printf("%.15g; number of evaluations: %i; error: %.15g\n",integral2,numevals,fabs(M_PI-integral2));
	printf("We see that Clenshaw-Curtis takes more evaluations but is more accurate.\n");
	
	int limits = 1e6; double result, error; numevals = 0;
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(limits);
	gsl_function F; 
	double function_newest(double x,void * params){numevals++; return p(x);}
	F.function = &function_newest;
	
	gsl_integration_qags(&F,a,b,delta,epsilon,limits,workspace,&result,&error);
	
	printf("QAGS from GSL library, the result, error and number of function evaluations are: \n%.15g\n%.15g\n%i,\nrespectively.\n",result,error,numevals);
	
	gsl_integration_workspace_free (workspace);
	
	// Part C 
	
	printf("\n\nPart C\n\n");
	
	numevals = 0; delta = 0.000001; epsilon = 0.000001;
	double functionC(double x) {numevals++; return gauss(x);} 
	
	double val = recadapter(functionC,-INFINITY, INFINITY, delta,epsilon,&errorestimate);
	printf("Using own routine:\n\n");
	printf("The gaussian function integrated from -infinity to infity is equal to %.10g with error estimate %.10g. The real value is sqrt(pi). The number of evaluations were %i.\n\n",val,errorestimate,numevals);
	printf("Using GSL routine:\n\n");
	
	numevals = 0;
	gsl_integration_workspace * workspace2 = gsl_integration_workspace_alloc(limits);
	gsl_function G; 
	double functionC2(double x,void * params) {numevals++; return gauss(x);}
	G.function = &functionC2;
	
	gsl_integration_qagi(&G,delta,epsilon,limits,workspace2,&result,&error);
	
	printf("The gaussian function integrated from -infinity to infity is equal to %.10g with error estimate %.10g. The real value is sqrt(pi). The number of evaluations were %i.\n\n",result,error,numevals);
	
	
	
return 0;
}








